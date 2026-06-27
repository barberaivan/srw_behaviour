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

softmax <- function(x) exp(x) / sum(exp(x))

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

# Import data -------------------------------------------------------------

bdata <- read.csv("data/srw_behaviour_data.csv",
                  stringsAsFactors = FALSE, sep = ",")

head(bdata)
# bdata_b <- read.csv("/home/ivan/Insync/Whales/Behavior/Behav_data.csv",
#                   stringsAsFactors = FALSE)
# 
# cbind(names(bdata), names(bdata_b))

# Time variable  ----------------------------------------------------------

# Chequeo que los datos estén ordenados por follow y por tiempo.

# variable de minutos desde las 00:00 en cada follow
bdata$min <- sapply(1:nrow(bdata), function(i) {
  sep <- strsplit(bdata$hour[i], ":") %>% unlist() %>%  as.numeric()
  return(sep[1] * 60 + sep[2])
})

head(bdata)

# order by año, follow and time
bdata <- bdata[order(bdata$year, bdata$follow, bdata$min), ]

# chequeo que no haya saltos mayores a 5 min
# fols_temp <- unique(bdata$follow)
# sapply(length(fols_temp), function(f) {
#   return(diff(bdata$min[bdata$follow == fols_temp[f]]) %>% unique)
# }) %>% unique # perfect

# variable t, nro de intervalo

bdata$t <- integer(nrow(bdata))
bdata$t[1] <- 1
for(i in 2:nrow(bdata)) {
  result <- ifelse(bdata$follow[i] == bdata$follow[i-1],
                   bdata$t[i-1] + 1,
                   1)
  bdata$t[i] <- result
}

# order 
bdata <- bdata[order(bdata$year, bdata$follow, bdata$t), ]

# is the last obs in the follow?
bdata$last <- c(sapply(1:(nrow(bdata)-1), function(i) {
  bdata$follow[i] != bdata$follow[i+1]
}), FALSE )

# View(bdata[, c("follow", "t", "last")])
# Attacks -----------------------------------------------------------------

# "PETREL" y algunos NA
unique(bdata$pac)
bdata$pac[bdata$pac == "PETREL"] <- 1

unique(as.character(bdata$pam)) # Hay NA
which(is.na(bdata$pam))         # En las mismas obs.
which(is.na(bdata$pac))         # Probably 3 follows with NA


# delete follows with NA in attack variables.
fols_with_na <- bdata$follow[is.na(bdata$pam) | is.na(bdata$pac)] %>% unique
nrow(bdata)
bdata <- bdata[!(bdata$follow %in% fols_with_na), ]
nrow(bdata)
# se pierden como 300 filas de 1995.

# ningún problema más con ataques.

bdata$pac <- as.numeric(bdata$pac)
bdata$pam <- as.numeric(bdata$pam)

# mother-calf specific attack variables

bdata$a <- as.numeric(bdata$pam > 0 | bdata$pac > 0)         # attack presence to mother or calf
bdata$am_only <- as.numeric(bdata$pam > 0 & bdata$pac == 0)  # attack presence only to mother
bdata$am_marg <- as.numeric(bdata$pam > 0)                   # attack presence to mother independent of what happens to calf
bdata$ac_only <- as.numeric(bdata$pac > 0 & bdata$pam == 0)  # attack presence only to mother
bdata$ac_marg <- as.numeric(bdata$pac > 0)                   # attack presence to mother independent of what happens to calf
bdata$nam <- bdata$pam                                       # number of attacks to mother
bdata$nac <- bdata$pac                                       # number of attacks to calf

# for simplicity
bdata$am <- bdata$am_marg

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
# !is.na(am[i]) & !is.na(nam[i]) & nam[i] > 0

# can't be last obs in the follow
# !last[i]
 
# behaviour has to be observed
# !is.na(behavm_clust[i])

# current mother behaviour must be observed
bdata$use <- sapply(1:nrow(bdata), function(i) {
                 with(bdata, 
                  !is.na(am[i]) & !is.na(nam[i]) &
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

# save deleted follows to include them in the marginal frequency estimation.
# saveRDS(fol_shorts, "behavior5_model_01_short follows deleted (N_use < 1) ids.rds")
# saveRDS(bdata, "behavior5_model_01_data without deleted follows (N_use < 1).rds")

# Subset data -------------------------------------------------------------

d <- bdata[, c("nam", "year", "follow", "behavm_clust")]
d <- d[bdata$use, ]
d$follow_id <- as.numeric(factor(d$follow))
N_fol <- length(unique(d$follow)) # 798
d$year_id <- as.numeric(factor(d$year), levels = unique(d$year))

nrow(d)
anyNA(d) # OK

# Data to fit model -------------------------------------------------------

# behaviour matrix
X_behavm <- model.matrix(nam ~ behavm_clust - 1, d)

# splines
B <- 5 # basis functions with intercept
bdata_year <- aggregate(nam ~ year, bdata, mean)
bs <- smoothCon(s(year, k = B, bs = "cr"), diagonal.penalty = TRUE,
                absorb.cons = TRUE, data = bdata_year)[[1]]
spline_mfit <- cbind(rep(1, nrow(bs$X)), bs$X) 

N_years <- max(d$year_id)

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
# gmc_code <- stan_model("model_gaviots_count_02.stan", verbose = TRUE)
# gmc <- sampling(
#   gmc_code, data = sdata, seed = 148974, refresh = 20,
#   #chains = 1, cores = 1, iter = 10,
#   chains = 6, cores = 6, iter = 2000,
#   control = list(adapt_delta = 0.95, max_treedepth = 15),
#   pars = c("gamma",
#            "sd_fol_l", "sd_year_l", "marginal_sd_l",
#            "e_year_l", "e_fol_l")
# )
# gmcs <- summary(gmc)[[1]]
# summary(gmcs[, "n_eff"])
# summary(gmcs[, "Rhat"])
# head(gmcs, n = 20)
# saveRDS(gmc, "models/attack/mothers_spline years/model_gaviots_count_02_samples.rds")
# saveRDS(gmcs, "attack/mothers_spline years/model_gaviots_count_02_summary.rds")

# 264.665 / 60 = 4.41 min
# model done

gmc <- readRDS("models/attack/mothers_spline years/model_gaviots_count_02_samples.rds")
gmcs <- readRDS("attack/mothers_spline years/model_gaviots_count_02_summary.rds")


# # Posterior predictive check ----------------------------------------------
# 
# # extract parameters
# gammas <- as.matrix(gmc, "gamma")
# e_fol_ls <- as.matrix(gmc, "e_fol_l")
# e_year_ls <- as.matrix(gmc, "e_year_l")
# nsim <- nrow(gammas)
# 
# nam_sim <- matrix(NA, length(d$nam), nsim)
# for(s in 1:nsim) {
#   print(s)
#   #s = 1
#   e_fol_l <- e_fol_ls[s, ] # second column
#   e_year_l <- e_year_ls[s, ]
#   
#   for(i in 1:sdata$N) {
#     lambda <- exp(
#       sdata$X_behavm[i, ] %*% gammas[s, ] + 
#       e_fol_l[sdata$follow_id[i]] + 
#       e_year_l[sdata$year_id[i]]
#     )
#     nam_sim[i, s] <- rtpois(1, lambda = lambda, a = 0)
#     # rtpois does not include a, but stan does.
#     # In stan the same truncation is T[1, ]
#   }
# }
# #anyNA(nam_sim)
# 
# saveRDS(nam_sim, "models/model_gaviots_count_02_nam_sim.rds")
# nam_sim <- readRDS("models/model_gaviots_count_02_nam_sim.rds")
# 
# nam_sim_tab <- table(as.numeric(nam_sim)) / sum(table(as.numeric(nam_sim)))
# 
# nam_obs_tab <- table(d$nam) / sum(table(d$nam))
# 
# par(mfrow = c(1, 2))
# plot(nam_sim_tab, main = "simulation")
# plot(nam_obs_tab, main = "observed", xlim = c(1, max(as.numeric(names(nam_sim_tab)))))
# par(mfrow = c(1, 1))
# # una belleza

# Behaviour effect (bad plot) ---------------------------------------------

# medians <- apply(exp(gammas), 2, hdmean, name = "lambda") %>% t %>% as.data.frame
# medians$behav <- behavs_m_clust
# medians$position <- rep(c("S", "UW"), K / 2)
# 
# ggplot(medians, aes(x = behav, y = lambda_mean, 
#                     ymin = lambda_lower, ymax = lambda_upper,
#                     colour = position)) +
#   geom_point() + 
#   geom_errorbar()
# # Muy poca variabilidad en nro de ataques

# Good plots --------------------------------------------------------------

# extract parameters
gamma_hat <- as.matrix(gmc, "gamma")
sd_marginal <- as.matrix(gmc, "marginal_sd_l")
nsim <- nrow(gamma_hat)

means <- gamma_hat
for(i in 1:nsim) {
  lambdas <- exp(gamma_hat[i, ] + 0.5 * sd_marginal[i] ^ 2)
  means[i, ] <- lambdas / (1 - exp(-lambdas))
}

means

means_summ <- apply(means, 2, hdmean) %>% t %>% as.data.frame
means_summ$behav <- behavs_m_clust

means_summ <- separate(means_summ, "behav", into = c("behav_raw", "position"), 
                       sep = "_", remove = FALSE)
means_summ$behav_raw <- factor(means_summ$behav_raw, levels = c("R", "ST", "FT", "IPA"),
                               labels = c("Rest", "Slow travel", "Fast travel", "In place activity"))
means_summ$position <- factor(means_summ$position, levels = c("S", "UW"),
                              labels = c("Surface", "Under water"))

ggplot(means_summ, 
       aes(x = behav_raw, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = position)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  #ylab("Pr(attack[t] = 1 | attack[t-1] = 1)") +
  ylab("Mean attack number when whales were attacked") +  
  xlab("Mother's behaviour")

ggsave("model_gaviots_count plot.png", 
       width = 17, height = 13, units = "cm")


# Plot with predictions and data ------------------------------------------

names(d)[which(names(d) == "behavm_clust")] <- "behav"
d2 <- separate_behav(d)

names(means_summ)[which(names(means_summ) == "position")] <- "exposure"

# Function to make boxplot without extreme points
# https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  #stats <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), method = 8)
  stats <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), method = 8)
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}

# custom boxplots


ggplot(means_summ, 
         aes(x = behav_raw, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
             colour = exposure)) +
  # data
  # stat_summary(fun.data = calc_stat, geom = "boxplot", alpha = 0,
  #              position = position_dodge2(preserve = "single", padding = 0.1,
  #                                         width = 0.4),
  #              width = 0.7, fatten = 1, lwd = 0.35,
  #              data = d2, mapping = aes(x = behav_raw, colour = exposure,
  #                                       y = nam),
  #              inherit.aes = F) +
  
  geom_jitter(data = d2, mapping = aes(x = behav_raw, colour = exposure,
                                        y = nam),
              position = position_dodge2(preserve = "single", padding = 0.1,
                                         width = 0.4),
              inherit.aes = F, alpha = 0.2) +
  
  # predictions
  geom_point(position = position_dodge2(width = 0.7), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.7)) +#,
                 #size = 3, alpha = 0.5) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Exposure") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  #ylab("Pr(attack[t] = 1 | attack[t-1] = 1)") +
  ylab("Attack number when whales were attacked") +  
  xlab("Mother's behaviour") + 
  scale_y_continuous(breaks = seq(0, 12, by = 2))
  

export_list <- list(
  count_means = means_summ,
  count_data = d2,
  notes = "predictions and data on mothers attack number"
)

saveRDS(export_list, 
        "attack/mothers_spline years/model_gaviots_count_02 predictions for plot.rds")
