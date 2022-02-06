options(scipen = 999)

library(tidyverse); theme_set(theme_bw())

library(rstan)           # Stan!
library(markovchain)     # compute steady state distribution for markov chains
library(bayestestR)      # highest density intervals
library(mgcv)            # spline basis set up
library(viridis)         # nice colors
library(posterior)       # summaries and diagnostics for posterior samples
library(ggh4x)           # facet_nested
#library(arm)             # logit # this generates problems with traceplot function

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

# Load data -----------------------------------------------------------------

setwd("/home/ivan/Insync/Whales/srw_behaviour") # not necesary here, but necessary to run from raw R.
bdata <- read.csv("srw_behaviour_data.csv",
                  stringsAsFactors = FALSE, sep = ",")

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

# Attacks -----------------------------------------------------------------

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

# mother-calf specific attack variables

bdata$a <- as.numeric(bdata$pam > 0 | bdata$pac > 0)         # attack presence to mother or calf
bdata$am_only <- as.numeric(bdata$pam > 0 & bdata$pac == 0)  # attack presence only to mother
bdata$am_marg <- as.numeric(bdata$pam > 0)                   # attack presence to mother independent of what happens to calf
bdata$ac_only <- as.numeric(bdata$pac > 0 & bdata$pam == 0)  # attack presence only to mother
bdata$ac_marg <- as.numeric(bdata$pac > 0)                   # attack presence to mother independent of what happens to calf
bdata$nam <- bdata$pam                                       # number of attacks to mother
bdata$nac <- bdata$pac                                       # number of attacks to calf
# just renamed pam and pac

# Behaviors ---------------------------------------------------------------

unique(bdata$behavm)[order(unique(bdata$behavm))] # it has missing data
#table(bdata$behavm) %>% plot

# reclassify

# define valid raw behaviors (I take out the undefined UWs or Ss)
behavs_m <- c("RUW", "RS", "SR", "GAL", "SR", "SRUW", "SRS",
              "TSUW", "TSS", #"TS", "TSY",
              "TMUW",
              "TMS", "TFS", "TFUW", #"TM", 
              "SA",
              "AUW")

table(bdata$behavm[!bdata$behavm %in% behavs_m]) # abundance of invalid behaviours
# Most are "UW" (under water). I could also use others, such as TS, but 
# too few observations fall in this kind of category, so it won't change the 
# model fit, and including them would make the likelihood harder to code.

# define simplified behaviors (m for mother, clust for clustered)
behavs_m_clust <- c("R_S", "R_UW", "ST_S", "ST_UW", "HE_S", "HE_UW") # (rest, slow travel and High Energy) X (S or UW)
behavs_m_short <- c("R", "ST", "HE")

K <- length(behavs_m_clust)

bdata$behavm_clust <- character(nrow(bdata))

# recode behaviours (deleted the undefined between S or UW)
bdata$behavm_clust[bdata$behavm %in% c("RS", "SR", "GAL", "SR", "SRS")] <- "R_S"
bdata$behavm_clust[bdata$behavm %in% c("RUW")] <- "R_UW"

bdata$behavm_clust[bdata$behavm %in% c("TSS")] <- "ST_S"
bdata$behavm_clust[bdata$behavm %in% c("TSUW")] <- "ST_UW"

bdata$behavm_clust[bdata$behavm %in% c("TMS", "TFS", "SA")] <- "HE_S"
bdata$behavm_clust[bdata$behavm %in% c("TMUW", "TFUW", "AUW")] <- "HE_UW"

bdata$behavm_clust <- factor(bdata$behavm_clust, levels = behavs_m_clust)
# summary(bdata$behavm_clust)

# What are the high energy behaviours?
# (table(bdata$behavm[bdata$behavm_clust == "HE_UW"]) / 
#   sum( table(bdata$behavm[bdata$behavm_clust == "HE_UW"]) )) %>% plot
# AUW       TFUW       TMUW 
# 0.04784191 0.13572543 0.81643266 
# Mostly travel under water

# (table(bdata$behavm[bdata$behavm_clust == "HE_S"]) / 
#     sum( table(bdata$behavm[bdata$behavm_clust == "HE_S"]) )) %>% plot
# SA        TFS        TMS 
# 0.76685834 0.01411396 0.21902771 
# Mostly surface activity

# Define uw column
# UW means that the observer was sure the whale was under water, but she/he 
# could not distinguish between R, ST or HE
bdata$uw <- 0
unique(bdata$behavm)[!unique(bdata$behavm) %in% behavs_m]
# UW and TUW will be UW, the others are completely NA
bdata$uw[bdata$behavm %in% c("UW", "TUW")] <- 1
sum(bdata$uw) / (nrow(bdata) - length(unique(bdata$follow))) 
# around 8 % is UW

# View(bdata[, c("follow", "t", "behavm_clust", "uw")])

# Find index of observations with previous and current behavior observed -----

# Computing the likelihood of these observations in Stan is much easier and 
# efficient than for others, because we can use matrix multiplication, avoiding
# a few loops.
# These observations are also the most informative, because behaviour transitions
# larger than one interval have a flatter likelihood.

# (consider UW unobserved for now.)

y_num <- as.numeric(bdata$behavm_clust)

# They must meet the following conditions:
# it must not be NA
# !is.na(Y[i]) & 

# it must not be a follow beggining
# bdata$t[i] != 1 

# previous behavior must not be NA

bdata$cons <- c(FALSE, 
                sapply(2:length(y_num), function(i) {
                  !is.na(y_num[i]) & (bdata$t[i] != 1) & !is.na(y_num[i-1])
                }))
bdata$cons_num <- as.numeric(bdata$cons)


# Filter follows by length ------------------------------------------------

# Assess follow length considering the most useful data 
# (here only cons == TRUE is considered)

fol_agg <- aggregate(cons_num ~ follow, bdata, sum)
#plot(table(fol_agg$cons_num))

# Delete follows with less than P consecutive observations
fol_shorts <- fol_agg$follow[fol_agg$cons_num < 4]
bdata <- bdata[!(bdata$follow %in% fol_shorts), ]

# Are there NAs in t = 1? 
(starting_behav <- summary(bdata$behavm_clust[bdata$t == 1]))
# yes, just a few
summary(bdata$behavm_clust[bdata$t == 1]) / length(bdata$behavm_clust[bdata$t == 1])

# Follow starts and ends ----------------------------------------------------

# (Stan doesn't allow ragged data structures)

# number of follows
(N_fol <- length(unique(bdata$follow))) # 1650

fol_agg <- aggregate(t ~ year + follow, FUN = length, data = bdata)
fol_agg <- fol_agg[order(fol_agg$year, fol_agg$follow), ]
bdata <- bdata[order(bdata$year, bdata$follow, bdata$t), ]

(mean(fol_agg$t *5)) # mean follow duration = 55.62 min
# plot(order(unique(bdata$follow))~ order(unique(fol_agg$follow))) # OK
# plot(table(fol_agg$t))

# identify each follow with an integer
bdata$follow_num <- factor(bdata$follow, levels = unique(bdata$follow)) %>% as.numeric

# number of intervals by follow
N_int <- fol_agg$t
# plot(table(N_int))

# follows end index
int_end <- cumsum(N_int)
# max(int_end) == nrow(bdata) # perfect

# follows beggining index
int_start <- int_end - N_int + 1

N <- sum(N_int)

# cbind(int_start, int_end)   # perfect


# Inputs for design matrix ---------------------------------------------------

# code behavior in matrix form, keeping NAs
y_mat <- model.matrix.lm(pam ~ behavm_clust - 1, data = bdata, na.action = "na.pass")
# turn NA into 0
y_mat[is.na(bdata$behavm_clust), ] <- 0 


# define index of current observed behaviors with previous observed behavior (cons == TRUE)
cons_id <- which(bdata$cons)
N_cons <- length(cons_id)

# Check if is there any NA at cons_id or before
sum(as.numeric(sapply(cons_id, function(i) {
  anyNA(bdata$behavm_clust[c(i, i-1)])
})))

# subset rows of y_mat at t-1 from cons_id
# This will be used as desing matrix, because the previous behaviour is 
# a predictor for the current.
y_mat_cons <- y_mat[cons_id - 1, ]
# anyNA(y_mat_cons)
# the same subsetting will be done in Stan for the z vector

# Numeric response variable -----------------------------------------------

y_vec <- as.numeric(bdata$behavm_clust)
# Stan doesn't accept NA, so we rename them
y_vec[is.na(y_vec)] <- K + 1


# Year data ---------------------------------------------------------------

bdata$year_factor <- factor(bdata$year, levels = unique(bdata$year))
bdata$year_num <- as.numeric(bdata$year_factor)
N_years <- length(levels(bdata$year_factor))

# Subset year_num to consecutive observations to compute matrix multiplication
# year by year
year_num_cons <- bdata$year_num[cons_id]

# Define starts and ends of each year to subset the design matrix
year_length <- table(year_num_cons) %>% unname
year_end <- cumsum(year_length)
year_start <- year_end - year_length + 1

# define period (1955 or 2000s) to vary z params
bdata$period_id <- 2
bdata$period_id[bdata$year == 1995] <- 1


# Basis functions for year effect ----------------------------------------

bdata_year <- aggregate(fam ~ year, bdata, mean)

B <- 5 # number of basis functions including the intercept
b <- B - 1 # without intercept
bs <- smoothCon(s(year, k = B, bs = "cr"), diagonal.penalty = TRUE,
                absorb.cons = TRUE, data = bdata_year,
                scale.penalty = TRUE)[[1]]
spline_mfit <- bs$X

# Dimentions check ---------------------------------------------------------

dim(y_mat_cons); length(y_vec[cons_id]); length(y_vec)

# Identify jumps in behavior data ------------------------------------------

# Jumps are sequences of unobserved behaviours (including UW) bounded by observed
# behaviours within a follow. 
# The jump start is the first NA observation, and
# the jump end is the first observed data after the NA sequence.
# We will need to compute the probabilities vector from the first NA (start)
# up to the jump end, which is observed (for t in start:end).
# The last probability vector is the likelihood for the jump-end observation.
# We take into account that UW missing observations can't be surface behaviours 
# (R_S, ST_S or HE_S).

jumps <- matrix(NA, 1, 2)
colnames(jumps) <- c("start", "end")

i = 1
while(i < nrow(bdata)) {
  
  # if it's observed, it might be a jump start
  if (!is.na(bdata$behavm_clust[i])) {
    
    start_candidate <- with(bdata,
      !is.na(behavm_clust[i]) &   # it's observed
      is.na(behavm_clust[i+1]) &  # following obs is NA
      follow[i] == follow[i+1]    # they belong to the same follow
    )

    if(start_candidate) {
      # has this jump landed somewhere?

      # the closest position to end:
      end_proposal <- i+2

      # it's an end of jump if it's observed at end_proposal and if
      # they belong to the same follow as start_candidate
      is_end <- !is.na(bdata$behavm_clust[end_proposal]) &
                bdata$follow[i] == bdata$follow[end_proposal]
      
      # If the end is not found, advance the end_proposal until the follow
      # is over.
      while(!is_end & (bdata$follow[i] == bdata$follow[end_proposal + 1])) {
        end_proposal <- end_proposal + 1
        is_end <- !is.na(bdata$behavm_clust[end_proposal])
      }

      if(is_end) {
        jumps <- rbind(jumps, c(i+1, end_proposal))
        
        
        # Advance i to the end_proposal - 1, because i will grow 1 3 lines below
        i <- end_proposal - 1
      }
    }
  }
  i <- i + 1  ## the problem is likely here: when it's end, it jumps to the next, without evaluating
  # whether it's also a start
}

# delete first row, which is NA
jumps <- jumps[-1, ]

# missing sequence length distribution
table(jumps[, 2] - jumps[, 1])
table(jumps[, 2] - jumps[, 1]) / nrow(jumps)
# most have 1 or 2 NAs

# order jumps data

jumps <- as.data.frame(jumps)
jumps$year <- bdata$year_num[jumps$start]
jumps$follow_id <- bdata$follow_num[jumps$start]
N_jumps <- nrow(jumps)

# problem at follow 8
# MS-KM-DT-RS34957HL

# occurs when an obs is both start and end of jump 



# Identify UW endings of follows ------------------------------------------

# These are follows that ended in unobserved behaviours containing at least
# one UW in the ending sequence.
# Examples:
#   R_UW, R_UW, NA, UW, NA
#   R_UW, R_UW, UW, UW, UW
#   R_UW, R_UW, NA, NA, UW
# In these cases, if t starts at 1 and it's the first follow in the data 
# end_start = 3 and end_end = 5. Subsequent end_start and end_end match their
# corresponding row in the dataset.
# Only the log_likelihood of UW will be added to the log_posterior.

uw_cols <- which(behavs_m_clust %in% c("R_UW", "ST_UW", "HE_UW"))
s_cols <- which(behavs_m_clust %in% c("R_S", "ST_S", "HE_S"))

uw_ends <- matrix(NA, 1, 2)
temp <- matrix(NA, 1, 2)
colnames(uw_ends) <- colnames(temp) <- c("start", "end")

i = 2
while(i <= nrow(bdata)) {
  
  # if it's na, it's a candidate for ending
  if (is.na(bdata$behavm_clust[i])) {
    
    f <- bdata$follow[i]
    time_i <- bdata$t[i]
    remaining <- with(bdata, behavm_clust[follow == f & t > time_i])
    last <- i + length(remaining)
    
    is_uw_ending <- with(bdata,
      time_i != 1 &                 # it's not a follow start
      !is.na(behavm_clust[i-1]) &   # previous behaviour was observed
      all(is.na(remaining)) #&       # remaining observations in the follow are missing
        
      #(uw[i:last] %>% sum) > 0      # there is at least one UW
      # this last condition is wrong, because I want to estimate the NAs, not only the UW
    )
    
    if(is_uw_ending) {
      temp[1, "start"] <- i
      i <- last # andvance i to the follow end, so the loop can move on later
      temp[1, "end"] <- i
      uw_ends <- rbind(uw_ends, temp)
    }
  }
  
  i <- i + 1 # if an end was found, it moves from last to the following follow
}

# delete first row, which is NA
uw_ends <- uw_ends[-1, ]

# uw ending length distribution
table(uw_ends[, 2] - uw_ends[, 1] + 1)
table(uw_ends[, 2] - uw_ends[, 1] + 1) / nrow(uw_ends)

uw_ends <- as.data.frame(uw_ends)
uw_ends$year <- bdata$year_num[uw_ends$start]
uw_ends$follow_id <- bdata$follow_num[uw_ends$start]

N_uw_ends <- nrow(uw_ends)

# bdata$behavm[uw_ends$start]
# bdata$behavm[uw_ends$end]

# Check jumps and uw_ends are correctly identified -------------------------

bdata$missing_obs <- character(nrow(bdata))
bdata$missing_obs[jumps$end] <- "jump_end"
bdata$missing_obs[jumps$start] <- "jump_start"
bdata$missing_obs[uw_ends$start] <- "uw_end_start"
bdata$missing_obs[uw_ends$end] <- "uw_end_end"

#View(bdata[, c("follow", "follow_num", "t", "behavm_clust", "uw", "missing_obs")])
# divine
# View(bdata[is.na(bdata$behavm_clust) & bdata$uw == 0, 
#            c("follow", "t", "behavm_clust", "uw", "missing_obs")])

sum(bdata$uw) / sum(is.na(bdata$behavm_clust)) # 93 % of NA are UW
sum(bdata$uw) / nrow(bdata) # 6.3 % of obs are UW


# Prediction data ----------------------------------------------------------
# computed outside stan to avois C stack size issues.

z_unique <- seq(0, 1, by = 0.025)
prev_states_text <- rep(behavs_m_clust, each = length(z_unique))
prev_states_temp <- rep(1:K, each = length(z_unique))
new_predmat_states <- matrix(0, length(prev_states_temp), K)
for(i in 1:nrow(new_predmat_states)) {
  new_predmat_states[i, prev_states_temp[i]] <- 1
}

z_rep <- rep(z_unique, K)
X_pred <- cbind(new_predmat_states, new_predmat_states * z_unique)

# This prediction matrix is used every year, so nrow(X_pred) is the number of 
# predictions by year 
year_length_pred <- nrow(X_pred) 
# dim(X_pred)


# Stan data ----------------------------------------------------------------

# index matrix to parameterize transition matrix with zeroes in the diagonal
ind_matrix <- matrix(NA, K-1, K)
for(j in 1:K) ind_matrix[, j] <- matrix(1:K, ncol = 1)[-j, ]

sdata <- list(
  # Ns and indexes
  N = N, N_cons = length(cons_id), K = K, 
  b = b, # number of basis without intercept
  N_fol = N_fol, int_start = int_start, int_end = int_end,
  N_years = N_years, 
  year_start = year_start, year_end = year_end, 
  year_length = year_length, 
  period_id = bdata$period_id,
  
  N_jumps = N_jumps,
  jump_start = jumps$start,
  jump_end = jumps$end,
  jump_year = jumps$year,
  
  N_uw_ends = N_uw_ends,
  end_start = uw_ends$start,
  end_end = uw_ends$end,
  end_year = uw_ends$year,
  
  uw = bdata$uw,
  
  cons_id = cons_id,
  cons_id_prev = cons_id - 1,
  
  ind_matrix = ind_matrix,
  s_cols = s_cols,

  # Variables for fitting
  y_mat = y_mat,
  y_vec = y_vec,
  # attack
  a = bdata$a,              # attack presence to mother or calf
  am_only = bdata$am_only,  # attack presence only to mother
  am_marg = bdata$am_marg,  # attack presence to mother independent of what happens to calf
  ac_only = bdata$ac_only,  # attack presence only to mother
  ac_marg = bdata$ac_marg,  # attack presence to mother independent of what happens to calf
  nam = bdata$nam,          # number of attacks to mother
  nac = bdata$nac,          # number of attacks to calf
  # spline
  spline_mfit = spline_mfit, # without intercept
  
  # priors
  
  # transition parameters intercepts sds
  prior_alpha_intercept_sd = 2,
  prior_beta_intercept_sd = 2,    
  
  # tmat spline sd prior (half-normal)
  prior_sd_spline_tmat_sd = 0.4, 
  # prior for linear term of the spline
  prior_sd_spline_lin_tmat_sd = 0.6,
  # tmat random effects sd prior (half-normal)
  prior_sd_years_sd = 1.5 / 2, 
  # These priors have to be relatively narrow to achieve convergence
  
  # z params
  prior_zpar_sd = 3,
  
  # Limits for increase and decay (mothers):
  # (still not fitted with restriction)    --------------VOLVER A CORRERLO CON ESTO
  #                                                      Y CON LA SPLINE EN ZPAR
  L_dec = 0.1745958,
  L_inc = 0.3187079 # better 
)
L_dec = 0.1745958#, reaches 0.1 after 1 h without attacks
L_inc = 0.4376587#, reaches 0.9 after 4 consecutive intervals with attack.
# L_inc = 0.3187079, reaches 0.9 after 6 consecutive intervals with attack.

# makes sense to let z increase little only with attacks to mothers because
# in years with much attack to calves it can reach higher values. 

str(sdata)

# N_fol # 1526
# range(fol_agg$t) # 5 to 53


# Model fit ---------------------------------------------------------------

# Prior simulations for the model in
# <behaviour_mothers_analysis_prior_simulations.R>

## compile and sample
# bm_code <- stan_model(file = "behaviour_mothers_model.stan", verbose = TRUE)
# warmup <- 1000
# N_samples <- 1200
# bm <- sampling(
#   bm_code, seed = 5454, refresh = 1,
#   chains = 10, cores = 10, iter = warmup + N_samples, warmup = warmup,
#   #chains = 1, cores = 1, iter = 4, # used to test the stan code and data
#   data = sdata,
#   control = list(adapt_delta = 0.95, max_treedepth = 15),
#   pars = c("sd_spline_tmat", "sd_years",
#            "decay",
#            "increase_b", "increase_b3",
#            "increase_m",
#            "increase_c", "increase_c3",
#            "kappa", "delta_m", "delta_c", "lambda_c", "decay_logit",
#            "alpha", "beta",                        # splines + random effects
#            "alpha_intercept", "beta_intercept",    # intercepts
#            "alpha_raneff",                         # random effects
#            "alpha_spline_coef", "beta_spline_coef",# splineS coefficients
#            "z_0")#, "z") # estimate without saving z because it uses too much memory
#                          # that fills up the 16 GB RAM when I try to extract it.
#                          # Compute z afterwards in R, sample-wise.
# )
# 21735.3 / 3600 = 6.037583 h running
# bms <- summary(bm)[[1]]
# head(bms)
# summary(bms[, "n_eff"])
# summary(bms[, "Rhat"])
## saveRDS(bm, "files/behaviour_mothers_model_samples_with_z.R")
# saveRDS(bm, "files/behaviour_mothers_model_samples.R")


# Load model
bm <- readRDS("files/behaviour_mothers_model_samples.R")

traceplot(bm, pars = "increase_m") + ylim(0, 1) + geom_hline(yintercept = L_inc) +
  geom_hline(yintercept = 0.319)
# ~0.319 sería para que llegue a 0.9 en 6 ataques seguidos
# en el 95 gran parte de la posterior queda fuera de rango.

traceplot(bm, pars = "decay") + ylim(0, 1) + geom_hline(yintercept = L_dec)
# con decay no hay problema.

traceplot(bm, pars = "increase_c") + ylim(0, 1)
traceplot(bm, pars = "increase_c3") + ylim(0, 1)

incm <- as.matrix(bm, pars = "increase_m")
str(incm)
sum(incm[, 1] < L_inc) / nrow(incm) # 39.375 % below lower limit
sum(incm[, 1] < 0.319) / nrow(incm) # 14.63333 % below 0.319

# veamos los z en el model 30. habra zs que se acercan a 1?
file.choose()
bm30 <- readRDS("/home/ivan/Whales_offline/behavior4_model_30_samples.R")
z30 <- as.matrix(bm30, pars = "z")
str(z30)
z30s <- apply(z30, 2, hdmean) %>% t %>% as.data.frame
z30m <- ecdf(z30s$mu_mean)
which(z30m == 0.9)
str(z30m)
z30m(1) - z30m(0.95) # almost 10 % reaches close to 1
z3095 <- ecdf(z30s$mu_mean[bdata$year == "1995"])
z3095(0.95) # 3 % is above 0.95

plot(density(z30s$mu_mean, from = 0, to = 1))
# sí, que los hay, los hay. 
plot(density(z30s$mu_mean[bdata$year == "1995"], from = 0, to = 1))
plot(density(z30s$mu_mean[bdata$year != "1995"], from = 0, to = 1, na.rm = TRUE))
# pero no en 1995.

traceplot(bm30, pars = "increase_m") + ylim(0, 1)+ geom_hline(yintercept = L_inc) +
  geom_hline(yintercept = 0.319)

incm30 <- as.matrix(bm30, pars = "increase_m")
str(incm30)
sum(incm30[, 1] < L_inc) / nrow(incm30) # 0.1858333 below lower limit
sum(incm30[, 1] < 0.319) / nrow(incm30) # 0.06016667 below 0.319

# Convergence diagnostics and effective sample size -----------------------

# check_hmc_diagnostics(bm)
# get_bfmi(bm)
# get_low_bfmi_chains(bm)
# # Everything is OK.
# 
# # summaries for model using posterior package
# 
# # turn samples into df
# bmdf <- as_draws_df(bm)
# 
# # exclude some parameters
# out1 <- grep("_intercept", colnames(bmdf))
# out2 <- grep("_raneff", colnames(bmdf))
# out3 <- which(colSums(bmdf) == 0) %>% unname
#  
# bmdf_subset <- bmdf[, -c(out1, out2, out3)]
#  
# bm_summary <- summarise_draws(bmdf_subset)
# # it takes long, save it:
# saveRDS(bm_summary, "files/behaviour_mothers_model_summary.R")
# 
# # Which are the worst-sampled parameters? 
# summary(bm_summary[, c("rhat", "ess_bulk", "ess_tail")])
# head(bm_summary[order(bm_summary$ess_bulk), ], n = 20)
# # min n_eff is 1818.


# Estimation of unobserved behaviours: JUMPS ------------------------------

bdata$behavm_hat <- NA
year_cols <- rep(1:N_years, K*K)

# To estimate unobserved behaviours we need z, but we didn't save it in the
# stanfit because it was a huge matrix that overloaded the 16 GB RAM. So here
# we create a function to extract z as a function of the follow (f) and the posterior
# sample (i).

# Extract posterior samples

alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
z0_samples <- as.matrix(bm, pars = "z_0")
zpar_samples <- as.matrix(bm, pars = c("kappa", "delta_m", "delta_c", 
                                       "lambda_c", "decay"))
N_samples <- nrow(alpha_samples) # 10 # to test, use 10

# objects used to compute z:
zpar_period_id <- list(
  first = c("kappa[1]", "delta_m[1]", "delta_c", "lambda_c", "decay[1]"),
  last = c("kappa[2]", "delta_m[2]", "delta_c", "lambda_c", "decay[2]")
)

# attack variables in matrix form
attack_matrix <- as.matrix(bdata[, c("a", "am_only", "ac_only", "nac", "ac_marg")])

# function to compute z
zcalc <- function(f, i) { # follow, iteration (numeric)
  #f = 1; i = 1
  z0 <- z0_samples[i, f]
  zpars <- zpar_samples[i, zpar_period_id[[ sdata$period_id[f] ]] ]
  names(zpars) <- c("kappa", "delta_m", "delta_c", "lambda_c", "decay")
  
  z <- numeric(N_int[f] + 1) # add one which is time = 0 (remove later)
  z[1] <- z0
  
  at <- rbind(rep(0, 5), attack_matrix[int_start[f] : int_end[f], ])
  
  # Loop
  for(t in (2:length(z))) {
    z[t] <- 
      z[t-1] +
      # increase
      at[t, "a"] * (1 - z[t-1]) * plogis(
        zpars["kappa"] + # intercept 1 attack both
        (-1) * zpars["delta_m"] * at[t, "am_only"] + # intercept 1 attack mother
        (-1) * zpars["delta_c"] * at[t, "ac_only"] + # intercept 1 attack calf
        zpars["lambda_c"] * (at[t, "nac"] - 1) * at[t, "ac_marg"] # effect of number of attacks to calf 
      ) +
      # decay
      (-1) * (1 - at[t, "a"]) * zpars["decay"] * z[t-1]
  }
  return(z[-1])
}


# list containing all the possible sequences in every jump with their marginal 
# posterior probability
jump_seq_hat <- vector(mode = "list", length = N_jumps)


for(j in 1:N_jumps) { 
  #j = 53
  print(paste("jump", j, "/", N_jumps))
  
  N_miss <- jumps$end[j] - jumps$start[j]
  
  N_trans <- N_miss + 1 # number of transitions to compute likelihood
  width <- N_miss + 2   # number of observations including previous and last observed (bounds)
  
  # list cotaining all possible observed values
  interval_list <- vector(mode = "list", length = N_miss)
  # get bdata$uw to know whether there were UW
  missing_uw <- bdata$uw[jumps$start[j] : (jumps$end[j] - 1)]
  for(m in 1:N_miss) {
    interval_list[[m]] <- 1:K
    if(missing_uw[m] == 1) interval_list[[m]] <- uw_cols
  }
  names(interval_list) <- paste(rep("miss", N_miss), 1:N_miss, sep = "_")
  
  seqs <- expand.grid(interval_list)
  N_seqs <- nrow(seqs)
  
  # define starting and ending behaviors (observed)
  start_behav <- y_vec[jumps$start[j] - 1] # remember that start is the first missing obs
  end_behav <- y_vec[jumps$end[j]]         # end is the observed behaviour
  
  # include starting and ending behaviors in the sequence matrix
  seqs_wide <- cbind(start = rep(start_behav, N_seqs),
                     seqs, 
                     end = rep(end_behav, N_seqs))
  
  # matrix to fill with the likelihood of every sequence for every posterior
  # sample
  seqs_probs <- matrix(NA, N_seqs, N_samples)
  
  # define year of the jump
  y <- jumps$year[j]
  
  # get follow
  f <- jumps$follow_id[j]
  
  # loop over posterior samples and sequences
  
  for(i in 1:N_samples) { # to N_samples
  #for(i in 1:10) { # to N_samples
    #print(paste("jump", j, "/", N_jumps, "- iteration", i))
    #i = 628
    #print(paste("jump ", j, ", sample ", i, sep = ""))
    
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    
    z_temp <- zcalc(f, i)
    # z_temp has length N_int[f], but we subset jumps according to their position
    # in the whole dataset, so follows don't start at 1 but at int_start[f].
    # to translate the jump$start[j] into the local z position (starting at 1),
    # we get
    start_local <- jumps$start[j] - int_start[f] + 1
    end_local <- jumps$end[j] - int_start[f] + 1
    # We need z from jump$start[j]-1 to jump$end[j]-1 because z[t-1] defines the 
    # transition matrix from t-1 to t.
    # But remember that the z vector is local, so we use start_local and end_local.
    z <- z_temp[(start_local - 1) : (end_local - 1)]
    
    # Create transition matrices
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      #int = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) %>% t
    }
    
    # Transition matrices must be modified if UW is in the missing sequence:
    # if a tmat leads to an UW behaviour, the surface columns (s_cols) must be set to 0, 
    # and the (modified) rows have to be normalized. The loop goes up to 
    # N_trans - 1 because the last transition leads always to a fully observed behaviour, 
    # so it can't be UW.
    for(tran in 1:(N_trans - 1)) { 
      if(missing_uw[tran] == 1) {
        tmats[, s_cols, tran] <- 0
        tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) %>% t
      }
    }
    
    # Loop over sequences
    for(s in 1:N_seqs) {
      #s = 1
      trans_probs <- numeric(N_trans) # vector to hold transition probabilities
      # Compute probabilities in each transition
      for(t in 1:N_trans) {
        trans_probs[t] <- tmats[seqs_wide[s, t], seqs_wide[s, t+1], t]
      }
      
      seqs_probs[s, i] <- trans_probs %>% log %>% sum %>% exp
    }# 
    
  } # end loop over posterior samples
  # it takes a lot of time by jump
  
  # compute probability of every sequence marginal to the posterior samples
  #str(seqs_probs)
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  #plot(seqs_probs_norm)
  #sum(seqs_probs_norm)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) %>% as.data.frame
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] %>% as.numeric
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
    
  # Fill most probable behaviours
  bdata$behavm_hat[jumps$start[j] : (jumps$end[j] - 1)] <- most_prob_seq_text 
  
  # Fill seqs result
  jump_seq_hat[[j]] <- seqs_result
}

# a long jump:
# jump_seq_hat[[58]]
# plot(jump_seq_hat[[58]]$p)
# check whether it worked
# follows with jumps
fwj <- bdata$follow[jumps$start] %>% unique
View(bdata[bdata$follow %in% fwj,
           c("follow", "follow_num", "t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "files/behaviour data with imputed jumps_R object_temporal.R")


# Estimation of unobserved behaviours: ENDS -------------------------------

# Missing behaviours for which there are no posterior observed behaviours.

# list containing all the possible sequences in every ending with their marginal 
# posterior probability
ending_seq_hat <- vector(mode = "list", length = N_uw_ends)

# Compute most likely sequences
for(j in 1:N_uw_ends) { 
  #j = 32
  print(paste("ending", j, "/", N_uw_ends))
  
  N_miss <- uw_ends$start[j]:uw_ends$end[j] %>% length
  N_trans <- N_miss # in this case, same transitions as missing intervals 
  width <- N_miss + 1 # number of observations including start
  
  # list cotaining all possible observed values
  interval_list <- vector(mode = "list", length = N_miss)
  # get bdata$uw to know whether there were UW
  missing_uw <- bdata$uw[uw_ends$start[j] : (uw_ends$end[j])]
  for(m in 1:N_miss) {
    interval_list[[m]] <- 1:K
    if(missing_uw[m] == 1) interval_list[[m]] <- uw_cols
  }
  names(interval_list) <- paste(rep("miss", N_miss), 1:N_miss, sep = "_")
  
  seqs <- expand.grid(interval_list)
  N_seqs <- nrow(seqs)
  
  # define starting behaviour (observed)
  start_behav <- y_vec[uw_ends$start[j] - 1] # remember that uw_ends$start is the first missing obs
  
  # include starting behaviour in the sequence matrix
  seqs_wide <- cbind(start = rep(start_behav, N_seqs),
                     seqs)
  
  if(start_behav == K + 1) stop(paste("starting behavior is NA at ending ", j, sep = ""))
  
  # matrix to fill with the likelihood of every sequence for every posterior
  # sample
  seqs_probs <- matrix(NA, N_seqs, N_samples)
  
  # define year of the ending
  y <- uw_ends$year[j]

  # get follow id
  f <- uw_ends$follow_id[j]
  
  # loop over posterior samples and sequences
  for(i in 1:N_samples) {
  #for(i in 1:10) {  
    #i = 1
    #print(paste("take off ", j, ", sample ", i, sep = ""))
    
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    
    z_temp <- zcalc(f, i)
    # z_temp has length N_int[f], but we subset jumps according to their position
    # in the whole dataset, so follows don't start at 1 but at int_start[f].
    # to translate the jump$start[j] into the local z position (starting at 1),
    # we get
    start_local <- uw_ends$start[j] - int_start[f] + 1
    end_local <- uw_ends$end[j] - int_start[f] + 1
    # We need z from jump$start[j]-1 to jump$end[j]-1 because z[t-1] defines the 
    # transition matrix from t-1 to t.
    # But remember that the z vector is local, so we use start_local and end_local.
    z <- z_temp[(start_local - 1) : (end_local - 1)]
    
    # make transition matrices array
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      tran = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) %>% t
      
      # Transitions to UW have zero probability of surface behaviours.
      if(missing_uw[tran] == 1) {
        tmats[, s_cols, tran] <- 0
        tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) %>% t
      }
    }
    
    # Loop over sequences
    for(s in 1:N_seqs) {
      #s = 1
      trans_probs <- numeric(N_trans) # vector to hold transition probabilities
      # Compute probabilities in each transition
      for(t in 1:N_trans) {
        trans_probs[t] <- tmats[seqs_wide[s, t], seqs_wide[s, t+1], t]
      }
      
      seqs_probs[s, i] <- trans_probs %>% log %>% sum %>% exp
    }# 
  } # end loop over posterior samples
  # it takes a lot of time by jump
  
  # compute probability of every ending sequence marginal to the posterior samples
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  #plot(seqs_probs_norm)
  #sum(seqs_probs_norm)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) %>% as.data.frame
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] %>% as.numeric
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
  
  # Fill most probable behaviours
  bdata$behavm_hat[uw_ends$start[j] : (uw_ends$end[j])] <- most_prob_seq_text 
  
  # Fill seqs result
  ending_seq_hat[[j]] <- seqs_result
  
}

# check whether it worked
# follows with uw_endings
fwuwe <- bdata$follow[uw_ends$start] %>% unique
View(bdata[bdata$follow %in% fwuwe,
           c("follow", "follow_num", "t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "files/behaviour data with imputed jumps and endings_R object_temporal.R")


# Estimation of unobserved behaviours: BEGINNING --------------------------

# I'll call them startlost

starts_na <- int_start[y_vec[int_start] == K + 1]

startlosts <- data.frame(start = starts_na, end = numeric(length(starts_na)),
                         behav_end = numeric(length(starts_na)))
# in this case, end position is the first observed behavior in the follow,
# similar to the jumps encoding.

# starting from the biginning, find where is the closest observed behavior, 
# register it and set its position as end.

N_startlosts <- nrow(startlosts)

for(t in 1:N_startlosts) {
  #t = 1
  start <- startlosts$start[t]
  
  post <- K + 1 # K + 1 is NA
  i <- 0
  while(post == K + 1) {
    post <- y_vec[start + (i + 1)]
    i <- i + 1
  }
  
  startlosts$end[t] <- start + i  # position of first observed behaviour
  startlosts$behav_end[t] <- post # first observed behaviour (numeric)
}

# add year
startlosts$year <- bdata$year_num[startlosts$start]
# add follow id 
startlosts$follow_id <- bdata$follow_num[startlosts$start]

bdata$missing_obs[startlosts$start] <- "startlost_start"
bdata$missing_obs[startlosts$end] <- "startlost_end"

# check visually:
# fwsl <- bdata$follow[startlosts$start]
# View(bdata[bdata$follow %in% fwsl, c("follow", "t", "behavm_clust", "missing_obs")])
# rownames don't match the positions in the table, but indexind is done by position,
# not by rownames, so it's not a problem.

# list containing all the possible sequences in every startlost with their marginal 
# probability
startlost_seq_hat <- vector(mode = "list", length = N_startlosts)

# Compute most probable sequences
for(j in 1:N_startlosts) {
  #j = 21
  print(paste("startlost", j, "/", N_startlosts))
  
  N_miss <- startlosts$end[j] - startlosts$start[j]
  N_trans <- N_miss # in this case, as many transitions as missing intervals 
  width <- N_miss + 1 # number of observations including end
  
  # list cotaining all possible observed values
  interval_list <- vector(mode = "list", length = N_miss)
  # get bdata$uw to know whether there were UW
  missing_uw <- bdata$uw[startlosts$start[j] : (startlosts$end[j] - 1)]
  for(m in 1:N_miss) {
    interval_list[[m]] <- 1:K
    if(missing_uw[m] == 1) interval_list[[m]] <- uw_cols
  }
  names(interval_list) <- paste(rep("miss", N_miss), 1:N_miss, sep = "_")
  
  seqs <- expand.grid(interval_list)
  N_seqs <- nrow(seqs)
  
  # define ending behavior
  end_behav <- y_vec[startlosts$end[j]]
  if(end_behav == K + 1) stop(paste("ending behavior is NA at startlost ", j, sep = ""))
  
  # include ending behavior in the sequence matrix
  seqs_wide <- cbind(seqs, end = rep(end_behav, N_seqs))
  
  # matrix to fill with the likelihood of every sequence for every posterior
  # sample
  seqs_probs <- matrix(NA, N_seqs, N_samples)
  
  # define year of the startlost
  y <- startlosts$year[j]

  # get follow id
  f <- startlosts$follow_id[j]
  
  # loop over posterior samples and sequences
  for(i in 1:N_samples) { # to N_samples
  #for(i in 1:100) {
    #i = 1
    #print(paste("startlost ", j, ", sample ", i, sep = ""))
    
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    
    z_temp <- zcalc(f, i)
    # z_temp has length N_int[f], but we subset jumps according to their position
    # in the whole dataset, so follows don't start at 1 but at int_start[f].
    # to translate the jump$start[j] into the local z position (starting at 1),
    # we get
    start_local <- startlosts$start[j] - int_start[f] + 1
    end_local <- startlosts$end[j] - int_start[f] + 1
    # We need z from jump$start[j] to jump$end[j]-1 because z[t-1] defines the 
    # transition matrix from t-1 to t.
    # But remember that the z vector is local, so we use start_local and end_local.
    z <- z_temp[(start_local) : (end_local - 1)]
    # (z0 is the real beggining, but we get it in another way)
    
    # get transition matrices array
    # these transitions go from the 1st missing behaviour (t = 1 in the follow)
    # up to the first observed behaviour.
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      #tran = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) %>% t
    }
    
    # Transitions to UW have zero probability of surface behaviours.
    # Notice this does not affect the last transition, because it leads to an
    # observed behaviour. So, if there is only 1 missing start, it's omitted.
    if(N_trans > 1) {
      for(tran in 1:(N_trans - 1)) {
        if(missing_uw[tran + 1] == 1) {
          tmats[, s_cols, tran] <- 0
          tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) %>% t
        }
      }
    }
    # Compute probability of first unobserved behavior unconditional to previous
    # behaviour (steady states distribution for z0)
    z0 <- z0_samples[i, f]
    
    tmat_0 <- apply(alpha + beta * z0, 1, softmax) %>% t
    mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat_0,
              name = "mc")  
    ss <- steadyStates(mc)
    # if the first behaviour is UW, correct probabilities to account for that.
    if(missing_uw[1] == 1) {
      ss[s_cols] <- 0
      ss <- normalize(ss)
    }
    
    # Loop over sequences
    for(s in 1:N_seqs) {
      #s = 1
      trans_probs <- numeric(N_trans) # vector to hold transition probabilities
      # Compute probabilities in each transition
      for(t in 1:N_trans) {
        # s = 1
        # t = 2
        trans_probs[t] <- tmats[seqs_wide[s, t], seqs_wide[s, t+1], t]
      }
      
      p_first <- ss[seqs[s, 1]]
      
      trans_probs_wide <- c(p_first, trans_probs)
      seqs_probs[s, i] <- trans_probs_wide %>% log %>% sum %>% exp
    }# 
  } # end loop over posterior samples
  
  # compute probability of every starting sequence marginal to the posterior samples
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) %>% as.data.frame
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] %>% as.numeric
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
  
  # Fill most probable behaviours
  bdata$behavm_hat[startlosts$start[j] : (startlosts$end[j] - 1)] <- most_prob_seq_text 
  
  # Fill seqs result
  startlost_seq_hat[[j]] <- seqs_result
  
}

# check whether it worked
# follows with lost starts
fwsl <- bdata$follow[startlosts$start] %>% unique
View(bdata[bdata$follow %in% fwsl,
           c("follow", "follow_num","t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "files/behaviour data with imputed jumps, endings and starts_R object_temporal.R")



# Merge and export observed and estimated behaviours ------------------------

# This is necessary to compute the observed behaviours proportion by year.
bdata$behavm_full <- bdata$behavm_clust
bdata$behavm_full[is.na(bdata$behavm_clust)] <- bdata$behavm_hat[is.na(bdata$behavm_clust)]

# check there are not missing behaviours
which(is.na(bdata$behavm_full)) # perfect

# View(bdata[bdata$follow %in% bdata$follow[is.na(bdata$behavm_full)],
#            c("follow", "t", "behavm_clust", "uw", "behavm_hat", "missing_obs")])

# save data with imputed behaviours
write.csv(bdata, "files/behaviour data with imputed behaviours.csv")
saveRDS(bdata, "files/behaviour data with imputed behaviours_R object.R")


hat_list <- list(jumps = jump_seq_hat, endings = ending_seq_hat, starts = startlost_seq_hat)
saveRDS(hat_list, "files/behaviour estimated missing sequences_R list.R")
# Estimation of unobserved behaviours finished.



hat_list <- readRDS("files/behaviour estimated missing sequences_R list.R")
str(hat_list)
hat_list$starts[[1]]

bdata2 <- readRDS("files/behaviour data with imputed behaviours_R object.R")

par(mfrow = c(1, 2))
plot(table(bdata2$behavm_clust[!is.na(bdata2$behavm_clust)]) / length(bdata2$behavm_clust[!is.na(bdata2$behavm_clust)]))
plot(table(bdata2$behavm_full) / length(bdata2$behavm_full))
par(mfrow = c(1, 1))

# Prediction under attacks scenarios --------------------------------------

# Compute behavior probabilities under a continuous attack scenario starting
# from undisturbed state and then the opposite.
# Starting behavior probabilities are steady states.

# samples for alpha and beta
alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
year_cols <- rep(1:N_years, K*K)
N_pred <- year_length_pred

# create prediction list
N_samples <- nrow(alpha_samples)
N_int_pred <- 31
years_unique <- unique(bdata$year)

pred_attack <- array(NA, dim = c(N_int_pred, K + 1, N_years, N_samples), 
                     dimnames = list(time = 0:(N_int_pred -1),
                                     behav = c("z", behavs_m_clust),
                                     year = years_unique,
                                     sample = 1:N_samples))

pred_no_attack <- pred_attack

attack_on <- rep(1, N_int_pred)
attack_off <- rep(0, N_int_pred)

increase <- as.matrix(bm, pars = "increase_b")
decay <- as.matrix(bm, pars = "decay")

#str(increase) # 2 columnas.
period_year <- c(1, rep(2, N_years - 1))

# here I will align z with the behavior probabilities in the same interval.
# So, now it's z[t] the one determining the transition probs from t-1 to t.
# The model was written differently to get a differente alignment to
# perform matrix multiplication.
# Here z depends on the previos attacks (which are constant in these two 
# scenarios).

### DON'T RUN, HEAVY LOOP ###

# pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
# for(i in 1:N_samples) {
#   #i = 1
#   
#   # loop over years
#   for(y in 1:N_years) {
#     #y = 1
#     alpha <- matrix(alpha_samples[i, which(year_cols == y)], K, K)
#     beta <- matrix(beta_samples[i, which(year_cols == y)], K, K)
#     
#     inc <- increase[i, period_year[y]]
#     dec <- decay[i, period_year[y]]
#     
#     # Attack scenario ---------------------------
#     
#     z <- 0 # starts calm
#     pred_attack[1, "z", y, i] <- z 
#     tmat <- apply(alpha + beta * z, 1, softmax) %>% t
#     mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat,
#               name = "mc")  
#     ss <- steadyStates(mc)
#     pred_attack[1, 2:(K+1), y, i] <- ss
#     
#     
#     for(t in 2:nrow(pred_attack)) {
#       #t = 2
#       pred_attack[t, "z", y, i] <- pred_attack[t-1, "z", y, i] + 
#         inc * (1 - pred_attack[t-1, "z", y, i])
#       tmat <- apply(alpha + beta * pred_attack[t, "z", y, i], 1, softmax) %>% t
#       
#       pred_attack[t, 2:(K+1), y, i] <- pred_attack[t-1, 2:(K+1), y, i] %*% tmat
#     }
#     
#     #pred_attack[, , i] 
#     
#     # No attack scenario --------------------------
#     
#     z <- 1 # starts upset
#     pred_no_attack[1, "z", y, i] <- z
#     tmat <- apply(alpha + beta * z, 1, softmax) %>% t
#     mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat,
#               name = "mc")  
#     ss <- steadyStates(mc)
#     pred_no_attack[1, 2:(K+1), y, i] <- ss
#     
#     for(t in 2:nrow(pred_no_attack)) {
#       #t = 2
#       pred_no_attack[t, "z", y, i] <- pred_no_attack[t-1, "z", y, i] -
#         dec * pred_no_attack[t-1, "z", y, i]
#       tmat <- apply(alpha + beta * pred_no_attack[t, "z", y, i], 1, softmax) %>% t
#       
#       pred_no_attack[t, 2:(K+1), y, i] <- pred_no_attack[t-1, 2:(K+1), y, i] %*% tmat
#     }
#     
#     #pred_no_attack[, , i] 
#   }
#   setTxtProgressBar(pb, i)
# }
# 
# pred_samples <- list(attack = pred_attack,
#                      no_attack = pred_no_attack)
# saveRDS(pred_samples, "model_30 predictions attack scenarios samples object.R")

pred_samples <- readRDS("model_30 predictions attack scenarios samples object.R")

at_temp <- apply(pred_samples$attack, 1:3, hdmean, name = "prob") %>% as.data.frame.table
at_no_temp <- apply(pred_samples$no_attack, 1:3, hdmean, name = "prob") %>% as.data.frame.table

at_temp$scenario <- "Attack"
at_no_temp$scenario <- "After attack"

predictions <- rbind(
  pivot_wider(at_temp, names_from = "Var1", values_from = "Freq"),
  pivot_wider(at_no_temp, names_from = "Var1", values_from = "Freq")
)

predictions$scenario <- factor(predictions$scenario, 
                               levels = c("Attack",
                                          "After attack"))

predictions$time <- as.numeric(as.character(predictions$time))
predictions$min <- predictions$time * 5
predictions$year <- factor(predictions$year, levels = years_unique)

str(predictions)
# Separate behav and uw

predictions <- separate(predictions, "behav", into = c("behav_raw", "position"), 
                        sep = "_", remove = FALSE)
predictions$behav_raw <- factor(predictions$behav_raw, levels = c("R", "ST", "HE", "z"),
                                labels = c("Rest", "Slow travel", "High energy", "z"))
predictions$position <- factor(predictions$position, levels = c("S", "UW"),
                               labels = c("Surface", "Under water"))

# guardo tabla
saveRDS(predictions, "model_30 predictions attack scenarios table object.R")
predictions <- readRDS("model_30 predictions attack scenarios table object.R")

# Plots
colors <- c("#440154", "#35b779")

filt1 <- with(predictions, 
  min <= 70 & 
  behav != "z"
)

ggplot(predictions[filt1, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = position)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_grid(rows = vars(year), cols = vars(behav_raw)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  ylab("Behavior probability") +
  xlab("Time (min)") +
  ggtitle("Behaviour in response to kelp gull micropredation") #+

ggsave("model_30 attack scenario predictions all years.png", height = 50, width = 20,
       units = "cm")
# nice


filt2 <- with(predictions, 
  min <= 70 & 
  behav != "z" &
  year %in% c("1995", "2004", "2009", "2015") 
)

ggplot(predictions[filt2, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = position)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_grid(rows = vars(year), cols = vars(behav_raw)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  ylab("Behavior probability") +
  xlab("Time (min)") +
  ggtitle("Behaviour in response to kelp gull micropredation") #+

ggsave("model_30 attack scenario predictions few years.png", height = 20, width = 20,
       units = "cm")


# Disturbed, undisturbed and observed by year -------------------------------

predictions2 <- predictions[predictions$time == 0, ]
# The starting condition in Attack is undisturbed, and in After attacks it's 
# disturbed
predictions2$scenario <- factor(predictions2$scenario, 
                                levels = c(levels(predictions$scenario), "Observed"),
                                labels = c("Undisturbed", "Disturbed", "Observed"))

predictions2$position <- factor(predictions2$position, 
                                levels = levels(predictions$position),
                                labels = c("S", "UW"))

predictions2$year_num <- as.numeric(as.character(predictions2$year))
predictions2$period <- factor("new", levels = c("old", "new"))
predictions2$period[predictions2$year_num == 1995] <- "old"

ghost1 <- predictions2[with(predictions2, period == "old" & behav != "z"), ]
ghost1$year_num <- 1994.5
ghost2 <- predictions2[with(predictions2, period == "old" & behav != "z"), ]
ghost2$year_num <- 1995.5
ghosts <- rbind(ghost1, ghost2)

# average data by year

# first, create matrix with full behaviour column, which includes estimated ones.
y_mat_full <- model.matrix.lm(pam ~ behavm_full - 1, data = bdata, na.action = "na.pass")

y_agg <- aggregate(y_mat_full ~ bdata$year, FUN = mean, na.rm = TRUE)
colnames(y_agg) <- c("year", behavs_m_clust)
y_agg$year_num <- as.numeric(as.character(y_agg$year))
y_agg <- pivot_longer(y_agg, which(colnames(y_agg) %in% behavs_m_clust), 
                      values_to = "prob_mean", names_to = "behav")
y_agg$scenario <- "Observed"
y_agg <- separate(y_agg, "behav", into = c("behav_raw", "position"), 
                  sep = "_", remove = FALSE)
y_agg$behav_raw <- factor(y_agg$behav_raw, levels = c("R", "ST", "HE"),
                          labels = c("Rest", "Slow travel", "High energy"))
y_agg$position <- factor(y_agg$position, levels = c("S", "UW"))
y_agg$prob_lower <- NA
y_agg$prob_upper <- NA
y_agg$period <- factor("new", levels = c("old", "new"))
y_agg[y_agg$year == 1995, "period"] <- "old"

# filter
filt3 <- with(predictions2, behav != "z")
colors <- c("#440154", "#35b779", "black")


# subset useful columns and merge.
p <- predictions2[filt3, c("year_num", "period", "behav_raw", "position", 
                           "scenario", "prob_mean", "prob_lower", "prob_upper")]
w <- y_agg[, c("year_num", "period", "behav_raw", "position", 
               "scenario", "prob_mean", "prob_lower", "prob_upper")]

ppp <- rbind(p, w)
ppp$scenario

ggplot(ppp, 
       aes(x = year_num, 
           y = prob_mean, ymin = prob_lower, ymax = prob_upper, 
           colour = scenario, shape = scenario)) + 
  geom_point(position = position_dodge2(width = 0.3), size = 3) +
  geom_linerange(position = position_dodge2(width = 0.3)) +
  geom_point(data = ghosts, mapping = aes(x = year_num, y = prob_lower), alpha = 0) +
  #geom_point(data = y_agg, mapping = aes(x = year_num, y = prob_mean)) + 
  #facet_grid(rows = vars(behav_raw, position), cols = vars(period), scales = "free", space = "free") +
  facet_nested(rows = vars(behav_raw, position), cols = vars(period), scales = "free", space = "free") +
  scale_color_manual(values = colors, name = "Condition") +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(1, 1, 20)) +
  guides(shape = FALSE) +
  ylab("Behavior probability") +
  xlab("Year") +
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1995, 2005, 2008, 2012, 2016)) +
  ggtitle("Temporal variation in behaviour") #+

ggsave("model_30 behaviours ts.png", height = 20, width = 20, unit = "cm")



# Kelp gull effects by year -----------------------------------------------

undisturbed <- pred_samples$attack[1, 2:(K+1), , ]
disturbed <- pred_samples$no_attack[1, 2:(K+1), , ]

y_agg2 <- aggregate(sdata$y_mat ~ bdata$year, FUN = mean, na.rm = TRUE)[, -1] %>% as.matrix %>% t

eff_pot_samples <- disturbed - undisturbed
eff_obs_samples <- eff_pot_samples
for(i in 1:dim(eff_obs_samples)[3]) {
  eff_obs_samples[, , i] <- y_agg2 - undisturbed[, , i]
}

eff_pot <- apply(eff_pot_samples, 1:2, hdmean, name = "eff") %>% as.data.frame.table
eff_obs <- apply(eff_obs_samples, 1:2, hdmean, name = "eff") %>% as.data.frame.table
eff_pot$effect_type <- "Potential"
eff_obs$effect_type <- "Observed"
effs <- rbind(eff_pot, eff_obs)

effs <- pivot_wider(effs, names_from = "Var1", values_from = "Freq")
effs <- separate(effs, "behav", into = c("behav_raw", "position"), 
                 sep = "_", remove = FALSE)
effs$behav_raw <- factor(effs$behav_raw, levels = c("R", "ST", "HE"),
                         labels = c("Rest", "Slow travel", "High energy"))
effs$position <- factor(effs$position, levels = c("S", "UW"))
effs$effect_type <- factor(effs$effect_type, levels = c("Potential", "Observed"))

effs$year <- as.numeric(as.character(effs$year))
effs$period <- factor("new", levels = c("old", "new"))
effs$period[effs$year == 1995] <- "old"

effghost1 <- effs[effs$year == 1995, ]
effghost1$year <- 1994.5
effghost2 <- effs[effs$year == 1995, ]
effghost2$year <- 1995.5
effghosts <- rbind(effghost1, effghost2)

colors_eff <- c("#440154", "#35b779")

ggplot(effs, 
       aes(x = year, 
           y = eff_mean, ymin = eff_lower, ymax = eff_upper, 
           colour = effect_type)) + 
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_point(position = position_dodge2(width = 0.3), size = 3, shape = 1) +
  geom_linerange(position = position_dodge2(width = 0.3)) +
  geom_point(data = effghosts, mapping = aes(x = year, y = eff_lower), alpha = 0) +
  facet_nested(rows = vars(behav_raw, position), cols = vars(period), 
               scales = "free", space = "free") +
  scale_color_manual(values = colors_eff, name = "Effect type") +
  ylab("Behavior probability difference (x - undisturbed)") +
  xlab("Year") +
  ggtitle("Kelp gull effects on whales behaviour") + 
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1995, 2005, 2008, 2012, 2016)) +
  ggtitle("Temporal variation in behaviour") #+

ggsave("model_30 behaviours effects ts.png", height = 20, width = 20, unit = "cm")


# Compute transition matrix by year ---------------------------------------

# samples for alpha and beta
alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
year_cols <- rep(1:N_years, K*K)

N_samples <- nrow(alpha_samples)
N_pred <- year_length_pred

# predictions array
pred_arr <- array(NA, dim = c(N_pred, K, N_years, N_samples),
                  dimnames = list(z_prev = paste(z_unique, prev_states_text, sep = "/"),
                                  behav_post = behavs_m_clust,
                                  year = 1:N_years,
                                  sample = 1:N_samples))

# Loop over posterior samples
pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
for(i in 1:N_samples) {
  #print(i)
  # loop over years
  for(y in 1:N_years) {
    #y = 1
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    coef <- rbind(alpha, beta)
    
    pred_arr[, , y, i] <- apply(X_pred %*% coef, 1, softmax) %>% t
  }
  setTxtProgressBar(pb, i)
}

# summarize and paste rows

tmat_data <- apply(pred_arr, 1:3, hdmean) %>% as.data.frame.table
tmat_data <- pivot_wider(tmat_data, names_from = "Var1", values_from = "Freq")

head(tmat_data)
tmat_data <- separate(tmat_data, z_prev, into = c("z", "behav_prev"), sep = "/")
tmat_data$z <- as.numeric(tmat_data$z)
tmat_data$year_real <- factor(sapply(1:nrow(tmat_data), function(i) {
  unique(bdata$year)[tmat_data$year[i]]
}), levels = unique(bdata$year))

tmat_data$behav_post <- factor(tmat_data$behav_post, levels = behavs_m_clust)
tmat_data$behav_prev <- factor(tmat_data$behav_prev, levels = behavs_m_clust)

tmat_data <- separate(tmat_data, "behav_post", into = c("behav_raw_post", "position_post"), 
                      sep = "_", remove = FALSE)
tmat_data <- separate(tmat_data, "behav_prev", into = c("behav_raw_prev", "position_prev"), 
                      sep = "_", remove = FALSE)


tmat_data$behav_raw_prev <- factor(tmat_data$behav_raw_prev, levels = c("R", "ST", "HE"),
                                   labels = c("Rest", "Slow travel", "High energy"))
tmat_data$position_prev <- factor(tmat_data$position_prev, levels = c("S", "UW"))

tmat_data$behav_raw_post <- factor(tmat_data$behav_raw_post, levels = c("R", "ST", "HE"),
                                   labels = c("Rest", "Slow travel", "High energy"))
tmat_data$position_post <- factor(tmat_data$position_post, levels = c("S", "UW"))

# tmat_data$period <- "2004 to 2018"
# tmat_data$period[tmat_data$year_real == "1995"] <- "1995"
# tmat_data$period <- factor(tmat_data$period, levels = c("1995", "2004 to 2018"))

saveRDS(tmat_data, "model_30_predictions tmat_data object.R")

# All years
# ggplot(tmat_data, aes(x = z, y = mu_mean, ymin = mu_lower, ymax = mu_upper, 
#                       colour = period, fill = year_real, group = year_real,
#                       alpha = period, size = period)) +
#   geom_line() +
#   scale_alpha_manual(values = c(1, 0.6)) +
#   scale_size_manual(values = c(1, 0.6)) +
#   scale_color_manual(values = c("red", "gray20")) +
#   #geom_ribbon(alpha = 0.2, color = NA) +
#   facet_grid(cols = vars(behav_post), rows = vars(behav_prev)) +
#   ggtitle("Model 30", subtitle = "spline + raneff") + 
#   #ylim(0, 0.85) +
#   ylab("Transition probabilities")
# # ggsave("model_30 tmat_01.png", height = 20, width = 20,
# #        units = "cm")

years_plot <- c("1995", "2018") 
# years_plot <- c("1995", "2010", "2018") # quantiles
# years_plot <- c("1995", "2004", "2008", "2018")
vir <- "D"
ggplot(tmat_data[tmat_data$year_real %in% years_plot, ], 
       aes(x = z, y = mu_mean, ymin = mu_lower, ymax = mu_upper, 
           colour = year_real, fill = year_real, group = year_real)) +
  geom_line() +
  scale_color_viridis(option = vir, discrete = TRUE, name = "Year", end = 0.7) +
  scale_fill_viridis(option = vir, discrete = TRUE, name = "Year", end = 0.7) +
  geom_ribbon(alpha = 0.2, color = NA) +
  facet_nested(cols = vars(behav_raw_post, position_post), 
               rows = vars(behav_raw_prev, position_prev)) +
  theme(panel.grid.minor = element_blank()) +
  ylab("Transition probability (row to column)")

ggsave("model_30 tmat_02.png", height = 20, width = 30,
       units = "cm")


