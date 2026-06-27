# Model with spline for years

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

# it must not be NA
# !is.na(am[i]) & !is.na(nam[i])

# it must not be a follow beggining
# bdata$t[i] != 1 

# previous value must not be NA
# !is.na(am[i-1]) & !is.na(nam[i-1])

# can't be last obs in the follow
# !last[i]

# current mother behaviour must be observed
bdata$use <- c(FALSE, 
               sapply(2:nrow(bdata), function(i) {
                 with(bdata, 
                  !is.na(am[i]) & !is.na(nam[i]) &
                  t[i] != 1 & 
                  !is.na(am[i-1]) & !is.na(nam[i-1]) &
                  !is.na(behavm_clust[i]) &
                  !last[i]
                 )
               }))
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


# Follow and year id -------------------------------------------------------

# follow num
bdata$follow_id <- as.numeric(factor(bdata$follow))
# number of follows
N_fol <- length(unique(bdata$follow)) # 1787

# Year id 
bdata$year_id <- as.numeric(factor(bdata$year), levels = unique(bdata$year))


# Subset data -------------------------------------------------------------

d <- bdata[, c("am", "nam", "year_id", "follow_id", "behavm_clust")]
d$am_prev <- NA
for(i in 2:nrow(d)) d$am_prev[i] <- bdata$am[i-1]
d <- d[bdata$use, ]

nrow(d)
anyNA(d) # OK


# Data to fit model -------------------------------------------------------

# behaviour matrix
X_behavm <- model.matrix(am ~ behavm_clust - 1, d)

# splines
B <- 5 # basis functions with intercept
bdata_year <- aggregate(am ~ year, bdata, mean)
bs <- smoothCon(s(year, k = B, bs = "cr"), diagonal.penalty = TRUE,
                absorb.cons = TRUE, data = bdata_year)[[1]]
spline_mfit <- cbind(rep(1, nrow(bs$X)), bs$X) 

# follow data, identify which ones have no attacks
d_fol <- aggregate(am ~ follow_id, d, sum)
str(d_fol)

follows_with <- which(d_fol$am > 0)
follows_without <- which(d_fol$am == 0)
#length(follows_with) + length(follows_without) == nrow(d_fol)

# rows with attack to evaluate nam
rows_with <- which(d$am > 0)
N_years <- max(d$year_id)

sdata <- list(
  N = nrow(d),
  N_fol = max(d$follow_id),
  N_years = max(d$year_id),
  K = ncol(X_behavm),
  B = B,
  
  am = d$am, 
  X_behavm = X_behavm,
  am_prev = d$am_prev,
  year_id = d$year_id,
  follow_id = d$follow_id,
  spline_mfit = spline_mfit,
  
  # priors (_p for probability, l for lambda)
  prior_intercepts_sd_p = 2,
  prior_lin_sd_p = 1.5,
  prior_spline_sd_p = 0.7,
  prior_years_sd_p = 1,
  prior_fol_sd_p = 1

)


# Model fit ---------------------------------------------------------------

# compilo y ajusto (o for occurrence)
# gmo_code <- stan_model("model_gaviots_occurrence_02.stan", verbose = TRUE)
# gmo <- sampling(
#   gmo_code, data = sdata, seed = 5468, refresh = 1,
#   #chains = 1, cores = 1, iter = 10,
#   chains = 6, cores = 6, iter = 2000,
#   control = list(adapt_delta = 0.95, max_treedepth = 15),
#   pars = c("sd_alpha", "sd_beta",
#            "sd_fol_p", "sd_year_p", "marginal_sd_p",
#            "alpha", "beta",
#            "e_year_p",
#            "e_fol_p")
# )
# gmos <- summary(gmo)[[1]]
# summary(gmos[, "n_eff"])
# summary(gmos[, "Rhat"])
# head(gmos, n = 20)
# saveRDS(gmo, "models/attack/mothers_spline years/model_gaviots_occurrence_02_samples.rds")
# saveRDS(gmos, "attack/mothers_spline years/model_gaviots_occurrence_02_summary.rds")

# 1290.28 / 60 = 21 min
# model done, perfect.

gmo <- readRDS("models/attack/mothers_spline years/model_gaviots_occurrence_02_samples.rds")
gmos <- readRDS("attack/mothers_spline years/model_gaviots_occurrence_02_summary.rds")

# # Posterior predictive check ----------------------------------------------
# 
# # extract parameters
# alphas <- as.matrix(gmo, "alpha")
# betas <- as.matrix(gmo, "beta")
# e_fol_ps <- as.matrix(gmo, "e_fol_p")
# e_year_ps <- as.matrix(gmo, "e_year_p")
# 
# nsim <- nrow(alphas)
# 
# # For attack ocurrence just compute p and summarize between follows.
# # (to be able to notice overdispersion)
# 
# p_sim <- matrix(NA, length(d$am), nsim)
# 
# for(s in 1:nsim) {
#   print(s)
#   #s = 1
#   alpha <- matrix(alphas[s, ], sdata$K, sdata$N_years)
#   beta <- matrix(betas[s, ], sdata$K, sdata$N_years)
#   #matrix(colnames(alphas), K, N_years)
#   e_fol_p <- e_fol_ps[s, ]
#   e_year_p <- e_year_ps[s, ]
#   
#   for(i in 1:sdata$N) {
#     p_sim[i, s] = plogis(
#       sdata$X_behavm[i, ] %*% alpha[, sdata$year_id[i]] * (1 - sdata$am_prev[i]) + 
#       sdata$X_behavm[i, ] %*% beta[, sdata$year_id[i]] * sdata$am_prev[i] +        
#       e_fol_p[sdata$follow_id[i]] + 
#       e_year_p[sdata$year_id[i]]
#     )
#   }
# }
# saveRDS(p_sim, "models/model_gaviots_occurrence_02_p_sim.rds")
# 
# p_obs <- aggregate(am ~ follow_id, d, mean)[, 2]
# fol_length <- aggregate(am ~ follow_id, d, length)[, 2]
# p_sim_agg <- aggregate(p_sim ~ follow_id, d, mean)[, -1] %>% as.matrix()
# 
# saveRDS(p_sim_agg, "models/model_gaviots_occurrence_02_p_sim_agg.rds")
# #p_sim_agg <- readRDS("models/model_gaviots_occurrence_02_p_sim_agg.rds")
# 
# longs <- which(fol_length >= 1)
# p_sim_agg_dens <- density(as.numeric(p_sim_agg[longs, ]), from = 0, to = 1)
# 
# plot(p_sim_agg_dens, xlim = c(0, 1), main = "attack occurrence ppcheck",
#      xlab = "p")
# lines(density(p_obs[longs], from = 0, to = 1), col = 2)
# 
# # zeros props
# am_sim <- rbinom(prod(dim(p_sim)), size = 1, p = p_sim)
# am_sim_tab <- table(am_sim) / prod(dim(p_sim))
# am_sim_tab; table(d$am) / nrow(d)
# # A la prop general le pega reeee bien.


# Miro parámetros crudamente -------------------------------------------------

# datap0 <- data.frame(behav = rep(behavs_m_clust, N_years),
#                      year = rep(unique(bdata$year), each = K))
# datap <- rbind(datap0, datap0)
# datap$prev <- rep(c("no attack", "attack"), each = nrow(datap0))
# 
# alpha_summ <- apply(plogis(alphas), 2, hdmean, name = "p") %>% t %>% as.data.frame
# beta_summ <- apply(plogis(betas), 2, hdmean, name = "p") %>% t %>% as.data.frame
# 
# predictions <- cbind(datap, rbind(alpha_summ, beta_summ))
# predictions$behav <- factor(predictions$behav, levels = behavs_m_clust)
# 
# ggplot(predictions, aes(x = year, y = p_mean, ymin = p_lower, ymax = p_upper,
#                         colour = behav, fill = behav)) + 
#   geom_line() + 
#   geom_ribbon(alpha = 0.3, color = NA) +
#   facet_wrap(vars(prev))
# 
# ggplot(predictions, aes(x = behav, y = p_mean, ymin = p_lower, ymax = p_upper,
#                         colour = prev, fill = prev)) + 
#   geom_point() + 
#   geom_errorbar() +
#   facet_wrap(vars(year))
# # parece que no ha cambiado a lo largo del tiempo.


# Good plots --------------------------------------------------------------

# extract parameters
alphas <- as.matrix(gmo, "alpha")
betas <- as.matrix(gmo, "beta")
sd_marginal <- as.matrix(gmo, "marginal_sd_p")
sd_fol <- as.matrix(gmo, "sd_fol_p")
e_year <- as.matrix(gmo, "e_year_p")

nsim <- nrow(alphas)

pdata <- expand.grid(behav = behavs_m_clust, 
                     year = unique(bdata$year),
                     am_prev = c(0, 1))
pdata$year_id <- as.numeric(as.factor(pdata$year))


pdata <- separate(pdata, "behav", into = c("behav_raw", "position"), 
                      sep = "_", remove = FALSE)
pdata$behav_raw <- factor(pdata$behav_raw, levels = behavs_m_short,
                          labels = behavs_labels)
pdata$position <- factor(pdata$position, levels = c("S", "UW"),
                         labels = c("Surface", "Under water"))

npred <- nrow(pdata)
pb_mat <- model.matrix(rep(1, npred) ~ behav - 1, pdata)

years_cols <- rep(1:N_years, each = K)

predictions_samples <- matrix(NA, npred, nsim)
pb <- txtProgressBar(min = 0, max = npred, style = 3)
for(i in 1:npred) {
  #i = 1
  y <- pdata$year_id[i]
  X <- pb_mat[i, ]


  alpha_local <- alphas[, years_cols == y]
  beta_local <- betas[, years_cols == y]


  lp <- alpha_local %*% X * (1 - pdata$am_prev[i]) +
        beta_local %*% X * pdata$am_prev[i] +
        e_year[, y]

  # Compute probability at logit scale, marginal to follow random effects
  p_local <- numeric(nsim)
  for(s in 1:nsim) {
    #s = 1
    p_local[s] <- momentsLogitnorm(mu = lp[s], sigma = sd_fol[s])["mean"]
  }

  predictions_samples[i, ] <- p_local

  setTxtProgressBar(pb, i)
}

saveRDS(predictions_samples, "models/attack/mothers_spline years/model_gaviots_occurrence_02 prediction samples.rds")
predictions_samples <- readRDS("models/attack/mothers_spline years/model_gaviots_occurrence_02 prediction samples.rds")


# summarize
predictions <- apply(predictions_samples, 1, hdmean, name = "p") %>% t %>% as.data.frame

pdata2 <- cbind(pdata, predictions)

ggplot(pdata2[pdata2$am_prev == 1, ], 
       aes(x = behav_raw, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = position)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  facet_wrap(vars(year), ncol = 3) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  #ylab("Pr(attack[t] = 1 | attack[t-1] = 1)") +
  ylab("Probability of keep attacking") +  
  xlab("Mother's behaviour")
# ggsave("model_gaviots_occurrence_02 plot keep attacking.png", 
#        width = 18, height = 21, units = "cm")

ggplot(pdata2[pdata2$am_prev == 0, ], 
       aes(x = behav_raw, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = position)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Position") +
  facet_wrap(vars(year), ncol = 3) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  #ylab("Pr(attack[t] = 1 | attack[t-1] = 1)") +
  ylab("Probability of start attacking") +  
  xlab("Mother's behaviour")
# ggsave("model_gaviots_occurrence_02 plot start attacking.png", 
#        width = 18, height = 21, units = "cm")


# Steady states -----------------------------------------------------------

head(pdata)
str(predictions_samples)

# separate samples for previous attack == 0 and == 1:
prev0_samples <- predictions_samples[pdata$am_prev == 0, ]
prev1_samples <- predictions_samples[pdata$am_prev == 1, ]

# create new data.frame and matrix for steady states:
pdata_steady <- pdata[pdata$am_prev == 0, -which(names(pdata) == "am_prev")]
samples_steady <- prev0_samples

# compute steady states:
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for(s in 1:nsim) {
  for(i in 1:nrow(pdata_steady)) {
    # i = 1; s = 1
    probs1 <- c(prev0_samples[i, s], prev1_samples[i, s])
    probs0 <- 1 - probs1
    tmat <- cbind(probs0, probs1)
    colnames(tmat) <- c("no attack", "attack")
    mc <- new("markovchain", states = colnames(tmat), transitionMatrix = tmat,
              name = "mc")
    ss <- steadyStates(mc)
    samples_steady[i, s] <- ss[1, "attack"]
  }
  setTxtProgressBar(pb, s)
}
saveRDS(samples_steady, "models/attack/mothers_spline years/model_gaviots_occurrence_02_steady samples.rds")
samples_steady <- readRDS("models/attack/mothers_spline years/model_gaviots_occurrence_02_steady samples.rds")

samples_steady[1:5, 1:3]

steady_summ <- apply(samples_steady, 1, hdmean, name = "mu") %>% t %>% as.data.frame
pred_steady <- cbind(pdata_steady, steady_summ)
names(pred_steady)[which(names(pred_steady) == "position")] <- "exposure"

# do steady states match the observed frequencies?
# first, consider aggregating over previous attack, too.
bdata$prev_at <- NA
for(i in 1:nrow(bdata)) {
  if(bdata$use[i]) {
    bdata$prev_at[i] <- bdata$am_marg[i-1]
  }
}

# aggregate data
agg1 <- aggregate(am_marg ~ behavm_clust + year + follow, bdata[bdata$use, ], mean)
agg2 <- aggregate(am_marg ~ behavm_clust + year, agg1, mean)
agg2 <- separate_behav(agg2, "behavm_clust")

ggplot(pred_steady, 
       aes(x = behav_raw, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 20,
             data = agg2, mapping = aes(x = behav_raw, y = am_marg, 
                                        colour = exposure),
             inherit.aes = FALSE) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Exposure") +
  facet_wrap(vars(year), ncol = 3) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  ylab("Kelp gull attack probability") +  
  xlab("mothers' behaviour") +
  ylim(0, 0.55)
# faily good fit.
ggsave("attack/mothers_spline years/model_gaviots_occurrence_02 plot steady by year.png",
       width = 17, height = 21, units = "cm")

# Facet by behaviour
ggplot(pred_steady, 
       aes(x = year, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_point(position = position_dodge2(width = 0.4), size = 3.5, shape = 20,
             data = agg2, mapping = aes(x = year, y = am_marg, 
                                        colour = exposure),
             inherit.aes = FALSE) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Exposure") +
  facet_wrap(vars(behav_raw), nrow = 4) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  ylab("Kelp gull attack probability") +  
  xlab("year") +
  ylim(0, 0.55)

ggsave("attack/mothers_spline years/model_gaviots_occurrence_02 plot steady by year-year x.png",
       width = 17, height = 21, units = "cm")



# Prediction averaging years (steady states) -------------------------------

predmeans <- aggregate(samples_steady ~ behav_raw + exposure,
                       data = pred_steady, FUN = mean)
names(predmeans)[1:6]
pdata_agg <- predmeans[, 1:2]
pred_agg <- apply(as.matrix(predmeans[, 3:(nsim + 2)]), 1, hdmean, name = "p") %>% 
  t %>% as.data.frame

pdata_avg <- cbind(pdata_agg, pred_agg)
agg3 <- aggregate(am_marg ~ behav_raw + exposure, agg2, mean)

ggplot(pdata_avg,
       aes(x = behav_raw, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  geom_point(agg3, mapping = aes(x = behav_raw, y = am_marg, colour = exposure),
             inherit.aes = F,
             position = position_dodge2(width = 0.4), size = 3.5, shape = 20) + 
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = "Exposure") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0)) +
  ylab("Attack probability") +  
  xlab("Mother's behaviour")

ggsave("attack/mothers_spline years/model_gaviots_occurrence_02 plot avg steady.png", 
       width = 17, height = 16, units = "cm")


export_list_occurrence <- list(
  pred_steady_years = pred_steady,
  data_avg_years = agg2,
  pred_steady_avg = pdata_avg,
  data_avg = agg3,
  notes = "_years has year included, and _avg is averaged over years"
)

saveRDS(export_list_occurrence,
  "attack/mothers_spline years/model_gaviots_occurrence_02 predictions and data for plots.rds")
