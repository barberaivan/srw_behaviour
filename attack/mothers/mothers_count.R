options(scipen = 999)

library(rstan)
library(tidyverse); theme_set(theme_bw())
library(arm)
library(plyr)      # revalue
library(markovchain)
library(bayestestR)
library(mgcv)
library(viridis)
library(posterior)
library(extraDistr) # rtpois

# Functions ---------------------------------------------------------------

softmax <- function(x) exp(x) / sum(exp(x))

normalize <- function(x) x / sum(x)

hdi_lower <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_low"]])
}

hdi_upper <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_high"]])
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

get_map <- function(x) {
  x <- as.numeric(x)
  map <- map_estimate(x)
  return(c(map = map))
}

# separate behaviour column inteo behav_type and exposure
separate_behav <- function(df, col_name = "behav") {
  df <- separate(df, col_name, into = c("behav_raw", "exposure"), 
                 sep = "_", remove = FALSE)
  df$behav_raw <- factor(df$behav_raw, levels = behavs_m_short,
                         labels = behavs_labels)
  df$exposure <- factor(df$exposure, levels = c("S", "UW"), 
                        labels = c("Surface", "Under water"))
  return(df)
}

# Load data -----------------------------------------------------------------

bdata <- read.csv("data/srw_behaviour_data.csv",
                  stringsAsFactors = FALSE, sep = ",")
# data is in a local folder that I won't upload yet.

# Observations every 5 minutes, discrete intervals. 
# behavm = mothers' behaviour at the interval beginning
# behavc = calves' behaviour
#   TS = travel slow,
#   TM = travel medium speed,
#   TF = travel fast,
#   UW = under water,
#   S = surface
#   SA = surface activity, 
#   AUW = under water activity
#   R = rest
#   GAL = resting in galeon position
#   NU = nursing (calves)
# pam = number of kelp gull attacks to mother during the interval
# pac = number of kelp gull attacks to calf
# fam = attack occurrence (binary) to mother
# fac = attack occurrence (binary) to calf 
# follow = id for continuous sequence of observations
# obs = observer id

# Time variable  ----------------------------------------------------------

# create variable minutes since 00:00 en every follow
bdata$min <- sapply(1:nrow(bdata), function(i) {
  sep <- strsplit(bdata$hour[i], ":") %>% unlist() %>%  as.numeric()
  return(sep[1] * 60 + sep[2])
})

# order by year, follow and time
bdata <- bdata[order(bdata$year, bdata$follow, bdata$min), ]

# check there are no jumps higher than 5 min within follows
# fols_temp <- unique(bdata$follow)
# sapply(length(fols_temp), function(f) {
#   return(diff(bdata$min[bdata$follow == fols_temp[f]]) %>% unique)
# }) %>% unique # perfect

# t variable, interval order (time)
bdata$t <- integer(nrow(bdata))
bdata$t[1] <- 1
for(i in 2:nrow(bdata)) {
  result <- ifelse(bdata$follow[i] == bdata$follow[i-1],
                   bdata$t[i-1] + 1,
                   1)
  bdata$t[i] <- result
}

# last observation in the follow
bdata$last <- FALSE
bdata$last[nrow(bdata)] <- TRUE
for(i in 1:(nrow(bdata) - 1)) {
  if(bdata$follow[i] != bdata$follow[i+1]) {
    bdata$last[i] <- TRUE
  }
}

# View(bdata[, c("follow", "t", "last")])
# attack variables --------------------------------------------------------

# Checks for missing data
unique(bdata$pac)
# "PETREL" and some NA
bdata$pac[bdata$pac == "PETREL"] <- 1

unique(as.character(bdata$pam)) # Also NA
which(is.na(bdata$pam))         # In the same rows
which(is.na(bdata$pac))         

which(is.na(bdata$fam))   # no NA in attack presence
which(is.na(bdata$fac))  

# Check follows with NA
fols_na <- bdata$follow[which(is.na(bdata$pac))] %>% unique
#for(f in fols_na) View(bdata[bdata$follow == f, ])

# VR-MS-KM-DT-RS34955ICH ## This follow has too many NA in behaviour because they didn't write UW/S. Delete it.
bdata <- bdata[!bdata$follow == "VR-MS-KM-DT-RS34955ICH", ]

## replace pam and pac with ones in the next follows:
# VR-MS34961JCH 
# VR-MS-KM34965JL
# VR-MS-KM34971BCH
# 1 is the most likely attack number in 1995, and observers recorded attack presence.
bdata[bdata$follow %in% fols_na[-1] & is.na(bdata$pac) & is.na(bdata$pam), 
      c("pam", "pac")] <- 1

bdata$pac <- as.numeric(bdata$pac)
bdata$pam <- as.numeric(bdata$pam)

# mother-calf specific attack variables (not all are used in this model)

bdata$a <- as.numeric(bdata$pam > 0 | bdata$pac > 0)         # attack presence to mother or calf
bdata$ab <- as.numeric(bdata$pam > 0 & bdata$pac > 0)        # attack presence ot both mother and calf
bdata$am_only <- as.numeric(bdata$pam > 0 & bdata$pac == 0)  # attack presence only to mother
bdata$am_marg <- as.numeric(bdata$pam > 0)                   # attack presence to mother independent of what happens to calf
bdata$ac_only <- as.numeric(bdata$pac > 0 & bdata$pam == 0)  # attack presence only to calf
bdata$ac_marg <- as.numeric(bdata$pac > 0)                   # attack presence to calf independent of what happens to calf
bdata$nam <- bdata$pam                                       # number of attacks to mother
bdata$nac <- bdata$pac                                       # number of attacks to calf
# just renamed pam and pac

# check attacks 
aborc <- as.numeric(bdata$ab == 1 | bdata$ac_only == 1)
table(aborc, bdata$am_only)
table(aborc, bdata$ac_marg)
table(bdata$ac_marg, bdata$am_only)

tritri <- bdata[, c("am_only", "ab", "ac_only")]
summary(rowSums(tritri[bdata$a == 1, ])) # these are mutually exclusive.


# Behaviors ---------------------------------------------------------------
unique(bdata$behavm)[order(unique(bdata$behavm))] # it has missing data
table(bdata$behavm)
# reclassifying

# define valid raw behaviors (now I take out the undefined UWs or S)
behavs_m <- c("RUW", "RS", "SR", "GAL", "SR", "SRUW", "SRS",
              "TSUW", "TSS", #"TS", "TSY",
              "TMUW",
              "TMS", "TFS", "TFUW", #"TM", 
              "SA",
              "AUW")
#son na 
# D TUW TS

# define simplified behaviors (m for mother, clust for clustered)
behavs_m_clust <- c("R_S", "R_UW", 
                    "ST_S", "ST_UW", 
                    "FT_S", "FT_UW",
                    "IPA_S", "IPA_UW") 
# (rest, slow travel, fast travel and in place activity) X (surface and under water)
behavs_m_short <- c("R", "ST", "FT", "IPA")
behavs_labels <- c("Rest", "Slow travel", "Fast travel", "In place activity")

uw_cols <- which(behavs_m_clust %in% c("R_UW", "ST_UW", "FT_UW", "IPA_UW"))
s_cols <- which(behavs_m_clust %in% c("R_S", "ST_S", "FT_S", "IPA_S"))

K <- length(behavs_m_clust)

bdata$behavm_clust <- character(nrow(bdata))

# recode behaviours (deleted the undefined between S or UW)
bdata$behavm_clust[bdata$behavm %in% c("RS", "SR", "GAL", "SR", "SRS")] <- "R_S"
bdata$behavm_clust[bdata$behavm %in% c("RUW")] <- "R_UW"

bdata$behavm_clust[bdata$behavm %in% c("TSS")] <- "ST_S"
bdata$behavm_clust[bdata$behavm %in% c("TSUW")] <- "ST_UW"

bdata$behavm_clust[bdata$behavm %in% c("TMS", "TFS")] <- "FT_S"
bdata$behavm_clust[bdata$behavm %in% c("TMUW", "TFUW")] <- "FT_UW"

bdata$behavm_clust[bdata$behavm %in% c("SA")] <- "IPA_S"
bdata$behavm_clust[bdata$behavm %in% c("AUW")] <- "IPA_UW"

bdata$behavm_clust <- factor(bdata$behavm_clust, levels = behavs_m_clust)
# summary(bdata$behavm_clust)

# Find index of useful observations ----------------------------------------

# They must meet all the following conditions:
# (these are different, for attack count!)

# it must not be NA, and has to be >0
# !is.na(am_marg[i]) & !is.na(nam[i]) & nam[i] > 0

# can't be last obs in the follow
# !last[i]
 
# behaviour has to be observed
# !is.na(behavm_clust[i])

# current calf behaviour must be observed
bdata$use <- sapply(1:nrow(bdata), function(i) {
                 with(bdata, 
                  !is.na(am_marg[i]) & !is.na(nam[i]) &
                  nam[i] > 0 &
                  !is.na(behavm_clust[i]) &
                  !last[i]
                 )
               })
bdata$use_num <- as.numeric(bdata$use)
sum(bdata$use_num)

# Assess follow length considering useful data 
# (here only use == TRUE is considered)

fol_agg <- aggregate(use_num ~ follow, bdata, sum)
#plot(table(fol_agg$use_num))

# Erase follows with less than 1 useful observations
fol_shorts <- fol_agg$follow[fol_agg$use_num < 1]
bdata <- bdata[!(bdata$follow %in% fol_shorts), ]
# nrow(bdata) 

# Subset data -------------------------------------------------------------

d <- bdata[, c("nam", "year", "follow", "behavm_clust")]
d <- d[bdata$use, ]
d$follow_id <- as.numeric(factor(d$follow))
N_fol <- length(unique(d$follow)) # 374
d$year_id <- as.numeric(factor(d$year), levels = unique(d$year))

nrow(d)
anyNA(d) # OK

# Data to fit model -------------------------------------------------------

# behaviour matrix
X_behavm <- model.matrix(nam ~ behavm_clust - 1, d)

#bdata_year <- aggregate(nam ~ year, bdata, mean)

sdata <- list(
  N = nrow(d),
  N_fol = max(d$follow_id),
  N_years = max(d$year_id),
  K = ncol(X_behavm),
  #B = B,
  
  nam = d$nam, 
  X_behavm = X_behavm,
  year_id = d$year_id,
  follow_id = d$follow_id,
  #spline_mfit = spline_mfit,
  
  # priors (_p for probability, l for lambda)
  prior_intercepts_sd_l = 5,
  # prior_lin_sd_l = 2,
  # prior_spline_sd_l = 1,
  prior_years_sd_l = 2,
  prior_fol_sd_l = 2
)
str(sdata)
with(d, table(nam, behavm_clust, year_id))
# Too little variation to estimate variation among years

# Model fit ---------------------------------------------------------------

## compilo y ajusto (c for count)
# gmc_code <- stan_model("attack/mothers/mothers_count.stan", verbose = TRUE)
# gmc <- sampling(
#   gmc_code, data = sdata, seed = 148974, refresh = 20,
#   #chains = 1, cores = 1, iter = 10,
#   chains = 6, cores = 6, iter = 2000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   pars = c("means", "gamma",
#            "sd_fol_l", "sd_year_l", "marginal_sd_l",
#            "e_year_l", "e_fol_l")
# )
# gmcs <- summary(gmc)[[1]]
# summary(gmcs[, "n_eff"])
# summary(gmcs[, "Rhat"])
# head(gmcs, n = 20)
# saveRDS(gmc, "models/attack/mothers/mothers_count_samples.rds")
# saveRDS(gmcs, "attack/mothers/mothers_count_summary.rds")

# 1687.77 / 60 = 28.13 min
# model done

gmc <- readRDS("models/attack/mothers/mothers_count_samples.rds")
gmcs <- readRDS("attack/mothers/mothers_count_summary.rds")


# Posterior predictive check ----------------------------------------------

# extract parameters
gammas <- as.matrix(gmc, "gamma")
e_fol_ls <- as.matrix(gmc, "e_fol_l")
e_year_ls <- as.matrix(gmc, "e_year_l")
nsim <- nrow(gammas)

nam_sim <- matrix(NA, length(d$nam), nsim)
for(s in 1:nsim) {
  print(s)
  #s = 1
  e_fol_l <- e_fol_ls[s, ]
  e_year_l <- e_year_ls[s, ]

  for(i in 1:sdata$N) {
    lambda <- exp(
      sdata$X_behavm[i, ] %*% gammas[s, ] +
      e_fol_l[sdata$follow_id[i]] +
      e_year_l[sdata$year_id[i]]
    )
    nam_sim[i, s] <- rtpois(1, lambda = lambda, a = 0)
    # rtpois does not include a, but stan does.
    # In stan the same truncation is T[1, ]
  }
}
saveRDS(nam_sim, "models/attack/mothers/mothers_count_nam_sim.rds")
nam_sim <- readRDS("models/attack/mothers/mothers_count_nam_sim.rds")

nam_sim_tab <- table(as.numeric(nam_sim)) / sum(table(as.numeric(nam_sim)))

nam_obs_tab <- table(d$nam) / sum(table(d$nam))

par(mfrow = c(1, 2))
plot(nam_sim_tab, main = "simulation")
plot(nam_obs_tab, main = "observed", xlim = c(1, max(as.numeric(names(nam_sim_tab)))))
par(mfrow = c(1, 1))
# una belleza

# Behaviour effect ----------------------------------------------------------

# extract parameters
# gamma_hat <- as.matrix(gmc, "gamma")
# sd_marginal <- as.matrix(gmc, "marginal_sd_l")
# nsim <- nrow(gamma_hat)
# 
# means <- gamma_hat
# for(i in 1:nsim) {
#   lambdas <- exp(gamma_hat[i, ] + 0.5 * sd_marginal[i] ^ 2)
#   means[i, ] <- lambdas / (1 - exp(-lambdas))
# }

# get observed means:
count_agg <- aggregate(nam ~ behavm_clust + year + follow, data = d, FUN = mean)
count_agg <- aggregate(nam ~ behavm_clust, data = count_agg, FUN = mean)
count_agg <- separate_behav(count_agg, "behavm_clust")

# arrange d
d <- separate_behav(d, "behavm_clust")


means <- as.matrix(gmc, "means") # it is ordered by behaviour
means_summ <- apply(means, 2, hdmean) %>% t %>% as.data.frame
means_summ$behavm_clust <- behavs_m_clust
means_summ <- separate_behav(means_summ, "behavm_clust")

# plot
ggplot(means_summ,
       aes(x = behav_raw, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = exposure)) +
  # raw data 
  # geom_boxplot(data = d, mapping = aes(x = behav_raw,
  #                                    y = nam, colour = exposure),
  #            position = position_dodge2(width = 0.4),
  #            inherit.aes = F,
  #            alpha = 0.5) +
  # predictions
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  
  # aggregated data (means)
  geom_point(data = count_agg, mapping = aes(x = behav_raw, 
                                             y = nam, colour = exposure),
             position = position_dodge2(width = 0.4),
             inherit.aes = F,
             size = 2, shape = 5) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "exposure") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  #ylab("Pr(attack[t] = 1 | attack[t-1] = 1)") +
  ylab("Mean attack number when mothers were attacked") +  
  xlab("Mother's behaviour")

ggsave("attack/mothers/mothers_count plot.png", 
       width = 17, height = 13, units = "cm")

# Using the behav + year + follow aggregation, the intervals include the observed
# means.


# Export predictions and means --------------------------------------------

rownames(means_summ) <- NULL
export_table <- left_join(means_summ, count_agg[, c("behavm_clust", "nam")], by = "behavm_clust")
names(export_table)[ncol(export_table)] <- "mu_obs"
export_table$mc <- "Mothers"
export_table$var <- "Attack number"

saveRDS(export_table, "attack/mothers/mothers_count_prediction table.rds")

