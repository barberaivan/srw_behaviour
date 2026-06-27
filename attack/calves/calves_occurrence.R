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
library(logitnorm)

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
  df$behav_raw <- factor(df$behav_raw, levels = behavs_c_short,
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

# Subset year for calves behaviour ----------------------------------------

yc <- unique(bdata[!is.na(bdata$behavc), "year"])
bdata <- bdata[bdata$year %in% yc, ]

# Behaviors ---------------------------------------------------------------
unique(bdata$behavc)[order(unique(bdata$behavc))] # it has missing data
#table(bdata$behavm) %>% plot

# reclassify

# define valid raw behaviors (I take out the undefined UWs or Ss)
behavs_c <- c("RUW", "RS", "GAL", "SR", "SRS", "NU",
              "TSUW", "TSS", 
              "TMUW",
              "TMS", "TFS", "TFUW",  
              "SA", "AS",
              "AUW")

table(bdata$behavc[!bdata$behavc %in% behavs_c]) # abundance of invalid behaviours
# Most are "UW" (under water). 

# define simplified behaviors (m for mother, clust for clustered)
behavs_c_clust <- c("R_S", "R_UW", 
                    "ST_S", "ST_UW", 
                    "FT_S", "FT_UW",
                    "IPA_S", "IPA_UW") 
# (rest, slow travel, fast travel and in place activity) X (surface and under water)
behavs_c_short <- c("R", "ST", "FT", "IPA")
behavs_labels <- c("Rest", "Slow travel", "Fast travel", "In place activity")

uw_cols <- which(behavs_c_clust %in% c("R_UW", "ST_UW", "FT_UW", "IPA_UW"))
s_cols <- which(behavs_c_clust %in% c("R_S", "ST_S", "FT_S", "IPA_S"))

K <- length(behavs_c_clust)

bdata$behavc_clust <- character(nrow(bdata))

# recode behaviours (deleted the undefined between S or UW)
bdata$behavc_clust[bdata$behavc %in% c("RS", "SR", "GAL", "SR", "SRS", "NU")] <- "R_S"
bdata$behavc_clust[bdata$behavc %in% c("RUW")] <- "R_UW"

bdata$behavc_clust[bdata$behavc %in% c("TSS")] <- "ST_S"
bdata$behavc_clust[bdata$behavc %in% c("TSUW")] <- "ST_UW"

bdata$behavc_clust[bdata$behavc %in% c("TMS", "TFS")] <- "FT_S"
bdata$behavc_clust[bdata$behavc %in% c("TMUW", "TFUW")] <- "FT_UW"

bdata$behavc_clust[bdata$behavc %in% c("SA", "AS")] <- "IPA_S"
bdata$behavc_clust[bdata$behavc %in% c("AUW")] <- "IPA_UW"

bdata$behavc_clust <- factor(bdata$behavc_clust, levels = behavs_c_clust)
# summary(bdata$behavc_clust)


# Find index of useful observations ----------------------------------------

# They must meet all the following conditions:

# it must not be NA
# !is.na(ac_marg[i]) & !is.na(nam[i])

# it must not be a follow beggining
# bdata$t[i] != 1 

# previous value must not be NA
# !is.na(ac_marg[i-1]) & !is.na(nam[i-1])

# can't be last obs in the follow
# !last[i]

# current mother behaviour must be observed
bdata$use <- c(FALSE, 
               sapply(2:nrow(bdata), function(i) {
                 with(bdata, 
                      !is.na(ac_marg[i]) & 
                      t[i] != 1 & 
                      !is.na(ac_marg[i-1]) &
                      !is.na(behavc_clust[i]) &
                      !last[i]
                 )
               }))
bdata$use_num <- as.numeric(bdata$use)
sum(bdata$use_num)

# Assess follow length considering useful data 
# (here only use == TRUE is considered)
fol_agg <- aggregate(use_num ~ follow, bdata, sum)
#plot(table(fol_agg$use_num))

# Erase follows with less than 2 useful observations
# (this way, we can estimate the follow random effect, and there are very few
# follows with just one useful obs)
fol_shorts <- fol_agg$follow[fol_agg$use_num < 2]
bdata <- bdata[!(bdata$follow %in% fol_shorts), ]
# nrow(bdata) 

# Follow and year id -------------------------------------------------------

# follow num
bdata$follow_id <- as.numeric(factor(bdata$follow, levels = unique(bdata$follow)))
# number of follows
N_fol <- length(unique(bdata$follow)) # 581

# Year id 
bdata$year_id <- as.numeric(factor(bdata$year, levels = unique(bdata$year)))


# Subset data -------------------------------------------------------------

d <- bdata[, c("ac_marg", "year_id", "follow_id", "behavc_clust")]
d$ac_marg_prev <- NA
for(i in 2:nrow(d)) d$ac_marg_prev[i] <- bdata$ac_marg[i-1]
d <- d[bdata$use, ]

nrow(d)
anyNA(d) # OK

d2 <- d
names(d2)[which(names(d2) == "behavc_clust")] <- "behav"
d2 <- separate_behav(d2)
saveRDS(d2, "models/attack/calves/calves_occurrence_data_subset.rds")

# Data to fit model -------------------------------------------------------

# behaviour matrix
X_behavc <- model.matrix(ac_marg ~ behavc_clust - 1, d)
N_years <- max(d$year_id)

sdata <- list(
  N = nrow(d),
  N_fol = max(d$follow_id),
  N_years = max(d$year_id),
  K = ncol(X_behavc),

  ac_marg = d$ac_marg, 
  X_behavc = X_behavc,
  ac_marg_prev = d$ac_marg_prev,
  year_id = d$year_id,
  follow_id = d$follow_id,

  # priors (_p for probability, l for lambda)
  prior_intercepts_sd_p = 2,
  prior_years_sd_p = 1,
  prior_fol_sd_p = 1
)
str(sdata)

# Model fit ---------------------------------------------------------------

# compilo y ajusto (o for occurrence)
# gmo_code <- stan_model("attack/calves/calves_occurrence.stan", verbose = TRUE)
# gmo <- sampling(
#   gmo_code, data = sdata, seed = 5468, refresh = 10,
#   #chains = 1, cores = 1, iter = 10,
#   chains = 10, cores = 10, iter = 2500, warmup = 1000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   pars = c("alpha", "beta",
#            "sd_fol_p", "sd_year_p", "marginal_sd_p",
#            "e_year_p", "e_fol_p")
# )
# # 1199.86 / 60 = 20 min
# # tiraba neff bajo antes, por eso tanta iter
# gmos <- summary(gmo)[[1]]
# summary(gmos[, "n_eff"])
# summary(gmos[, "Rhat"])
# head(gmos, n = 20)
# saveRDS(gmo, "models/attack/calves/calves_occurrence_samples.rds")
# saveRDS(gmos, "attack/calves/calves_occurrence_summary.rds")

# 1290.28 / 60 = 21 min
# model done, perfect.

gmo <- readRDS("models/attack/calves/calves_occurrence_samples.rds")
gmos <- readRDS("attack/calves/calves_occurrence_summary.rds")

# # Posterior predictive check ----------------------------------------------

# extract parameters
alphas <- as.matrix(gmo, "alpha")
betas <- as.matrix(gmo, "beta")
e_fol_ps <- as.matrix(gmo, "e_fol_p")
e_year_ps <- as.matrix(gmo, "e_year_p")

nsim <- nrow(alphas)

# For attack ocurrence just compute p and summarize between follows.
# (to be able to notice overdispersion)

p_sim <- matrix(NA, length(d$ac_marg), nsim)
# for(s in 1:nsim) {
#   print(s)
#   #s = 1
#   alpha <- alphas[s, ]
#   beta <- betas[s, ]
#   e_fol_p <- e_fol_ps[s, ]
#   e_year_p <- e_year_ps[s, ]
# 
#   for(i in 1:sdata$N) {
#     p_sim[i, s] = plogis(
#       sdata$X_behavc[i, ] %*% alpha * (1 - sdata$ac_marg_prev[i]) +
#       sdata$X_behavc[i, ] %*% beta * sdata$ac_marg_prev[i] +
#       e_fol_p[sdata$follow_id[i]] +
#       e_year_p[sdata$year_id[i]]
#     )
#   }
# }
# saveRDS(p_sim, "models/attack/calves/calves_occurrence_p_sim.rds")

p_obs <- aggregate(ac_marg ~ follow_id, d, mean)[, 2]
fol_length <- aggregate(ac_marg ~ follow_id, d, length)[, 2]
# p_sim_agg <- aggregate(p_sim ~ follow_id, d, mean)[, -1] %>% as.matrix()

# saveRDS(p_sim_agg, "models/attack/calves/calves_occurrence_p_sim_agg.rds")
p_sim_agg <- readRDS("models/Attack analysis/calves_8 behavs/model_gaviots_occurrence_02_p_sim_agg.rds")

longs <- which(fol_length > 5)
p_sim_agg_dens <- density(as.numeric(p_sim_agg[longs, ]), from = 0, to = 1)

plot(p_sim_agg_dens, xlim = c(0, 1), main = "attack occurrence ppcheck",
     xlab = "p")
lines(density(p_obs[longs], from = 0, to = 1), col = 2)
# puede estar un poco sobredisperso, pero no hay mucho que hacer.

# zeros props
ac_marg_sim <- rbinom(prod(dim(p_sim)), size = 1, p = p_sim)
ac_marg_sim_tab <- table(ac_marg_sim) / prod(dim(p_sim))
ac_marg_sim_tab; table(d$ac_marg) / nrow(d)
# A la prop general le pega reeee bien.

# Behaviour effects --------------------------------------------------------------

# extract parameters
alphas <- as.matrix(gmo, "alpha")
betas <- as.matrix(gmo, "beta")
sd_marginal <- as.matrix(gmo, "marginal_sd_p")
sd_fol <- as.matrix(gmo, "sd_fol_p")
e_year <- as.matrix(gmo, "e_year_p")

nsim <- nrow(alphas)

pdata <- expand.grid(behav = behavs_c_clust, 
                     ac_marg_prev = c(0, 1))
pdata <- separate_behav(pdata, "behav")

npred <- nrow(pdata)
pb_mat <- model.matrix(rep(1, npred) ~ behav - 1, pdata)
# 
# 
# years_cols <- rep(1:N_years, each = K)

# predictions_samples <- matrix(NA, npred, nsim)
# pb <- txtProgressBar(min = 0, max = npred, style = 3)
# for(i in 1:npred) {
#   #i = 1
#   X <- pb_mat[i, ]
#   
#   lp <- alphas %*% X * (1 - pdata$ac_marg_prev[i]) +
#         betas %*% X * pdata$ac_marg_prev[i]
# 
#   # Compute probability marginal to random effects
#   p_local <- numeric(nsim)
#   for(s in 1:nsim) {
#     #s = 1
#     p_local[s] <- momentsLogitnorm(mu = lp[s], sigma = sd_marginal[s])["mean"]
#   }
# 
#   predictions_samples[i, ] <- p_local
# 
#   setTxtProgressBar(pb, i)
# }
# saveRDS(predictions_samples, "models/attack/calves/calves_occurrence_prediction samples.rds")
predictions_samples <- readRDS("models/attack/calves/calves_occurrence_prediction samples.rds")

# summarize
predictions <- apply(predictions_samples, 1, hdmean, name = "p") %>% t %>% as.data.frame

pdata <- separate(pdata, "behav")

pdata2 <- cbind(pdata, predictions)
pdata2$attack_prev <- NA
pdata2$attack_prev[pdata2$ac_marg_prev == 1] <- "Keep attacking"
pdata2$attack_prev[pdata2$ac_marg_prev == 0] <- "Start attacking"
pdata2$attack_prev <- factor(pdata2$attack_prev, 
                             levels = c("Start attacking", "Keep attacking"))

ggplot(pdata2, 
       aes(x = behav_raw, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  facet_wrap(vars(attack_prev), ncol = 1) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  ylab("Kelp gull attack probability") +  
  xlab("Calves' behaviour")
# ggsave("attack/calves/calves_occurrence_plot.png", 
#        width = 18, height = 21, units = "cm")


# Summarizing attack occurrence (steady states) ------------------------------

head(pdata)
str(predictions_samples)

# separate samples for previous attack == 0 and == 1:
prev0_samples <- predictions_samples[pdata$ac_marg_prev == 0, ]
prev1_samples <- predictions_samples[pdata$ac_marg_prev == 1, ]

# create new data.frame and matrix for steady states:
pdata_steady <- pdata[pdata$ac_marg_prev == 0, -which(names(pdata) == "ac_marg_prev")]
samples_steady <- prev0_samples

# compute steady states:
# pb <- txtProgressBar(min = 0, max = nsim, style = 3)
# for(s in 1:nsim) {
#   for(i in 1:nrow(pdata_steady)) {
#     #i = 1; s = 1
#     probs1 <- c(prev0_samples[i, s], prev1_samples[i, s])
#     probs0 <- 1 - probs1
#     tmat <- cbind(probs0, probs1)
#     colnames(tmat) <- c("no attack", "attack")
#     mc <- new("markovchain", states = colnames(tmat), transitionMatrix = tmat,
#               name = "mc")  
#     ss <- steadyStates(mc)
#     samples_steady[i, s] <- ss[1, "attack"]
#   }
#   setTxtProgressBar(pb, s)
# }
# saveRDS(samples_steady, "models/attack/calves/calves_steady samples.rds")
samples_steady <- readRDS("models/attack/calves/calves_steady samples.rds")

steady_summ <- apply(samples_steady, 1, hdmean, name = "mu") %>% t %>% as.data.frame
pred_steady <- cbind(pdata_steady, steady_summ)

# do steady states match the observed frequencies?
# first, consider aggregating over previous attack, too.
bdata$prev_at <- NA
for(i in 1:nrow(bdata)) {
  if(bdata$use[i]) {
    bdata$prev_at[i] <- bdata$am_marg[i-1]
  }
}

# aggregate data
agg1 <- aggregate(ac_marg ~ behavc_clust + year + follow + prev_at, bdata[bdata$use, ], mean)
agg2 <- aggregate(ac_marg ~ behavc_clust, agg1, mean)
agg2 <- separate_behav(agg2, "behavc_clust")

ggplot(pred_steady, 
       aes(x = behav_raw, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 4,
             data = agg2, mapping = aes(x = behav_raw, y = ac_marg, 
                                        colour = exposure),
             inherit.aes = FALSE) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  # facet_wrap(vars(attack_prev), ncol = 1) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  ylab("Kelp gull attack probability") +  
  xlab("mothers' behaviour") +
  ylim(0, 0.55)

# Using steady states the means are more uncertain, but very similar as computing
# the weighted average. In both cases, there are large discrepancies, mainly
# for IPA_UW


# Sumarizing attack occurrence (weighted means) ------------------------------

head(pdata)
str(predictions_samples)

# separate samples for previous attack == 0 and == 1:
prev0_samples <- predictions_samples[pdata$ac_marg_prev == 0, ]
prev1_samples <- predictions_samples[pdata$ac_marg_prev == 1, ]

# create new data.frame and matrix for weighted means:
pdata_w <- pdata[pdata$ac_marg_prev == 0, -which(names(pdata) == "ac_marg_prev")]

# compute weights:
a_w <- d$ac_marg_prev %>% mean
samples_w <- prev0_samples * (1 - a_w) + prev1_samples * a_w
 
samples_w[1:5, 1:3]

w_summ <- apply(samples_w, 1, hdmean, name = "y") %>% t %>% as.data.frame
pred_w <- cbind(pdata_w, w_summ)

ggplot(pred_w, 
       aes(x = behav_raw, y = y_mean, ymin = y_lower, ymax = y_upper,
           colour = position)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 4,
             data = agg3, mapping = aes(x = behav_raw, y = ac_marg, 
                                        colour = position),
             inherit.aes = FALSE) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  # facet_wrap(vars(attack_prev), ncol = 1) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  ylab("Kelp gull attack probability") +  
  xlab("Calves' behaviour") +
  ylim(0, 0.55)
# ggsave("attack/calves/calves_occurrence steady_plot.png", 
#        width = 17, height = 12, units = "cm")


# Attack mean number --------------------------------------------------------

# conditional attack number (given they are attacked)
gmc <- readRDS("models/attack/calves/calves_count_samples.rds")
atnum <- as.matrix(gmc, "means") %>% t
atnum_summ <- apply(atnum, 1, hdmean, name = "y") %>% t %>% as.data.frame
atnum_plot <- cbind(pdata_steady, atnum_summ)
set.seed(2)
take_id <- sample(1:nsim, size = ncol(atnum), replace = F)

# marginal attack number.
# Exptectation for Hurdle model taken from Cameron & Trivedi (2013) p. 137 and from 
# https://data.princeton.edu/wws509/notes/countmoments
atnum_marg <- samples_steady[, take_id] * atnum
atnum_marg_summ <- apply(atnum_marg, 1, hdmean, name = "y") %>% t %>% as.data.frame
atnum_marg_plot <- cbind(pdata_steady, atnum_marg_summ)

# merge and plot
pred_steady$var <- "Attack probability (steady states)"
atnum_plot$var <- "Attack number when attacked"
atnum_marg_plot$var <- "Attack number (unconditional)"

atplot <- rbind(pred_steady, atnum_plot, atnum_marg_plot)
atplot$var <- factor(atplot$var, levels = c(pred_steady$var[1],
                                            atnum_plot$var[1],
                                            atnum_marg_plot$var[1]))

ggplot(atplot, 
       aes(x = behav_raw, y = y_mean, ymin = y_lower, ymax = y_upper,
           colour = position)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  facet_wrap(vars(var), ncol = 1, scales = "free_y") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 15)) +
  ylab("Kelp gull attack probability") +  
  xlab("Calves' behaviour")

ggsave("attack/calves/calves_complete plot.png", 
       width = 17, height = 20, units = "cm")


# Export predictions and means --------------------------------------------

rownames(pred_steady) <- NULL
export_table <- left_join(pred_steady, agg2[, c("behav_raw", "exposure", "ac_marg")],
                          by = c("behav_raw", "exposure"))
names(export_table)[ncol(export_table)] <- "mu_obs"
export_table$mc <- "Calves"
export_table$var <- "Attack probability"

saveRDS(export_table, "attack/calves/calves_occurrence_prediction table.rds")
