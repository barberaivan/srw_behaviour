options(scipen = 999)

library(rstan)           # Stan!
library(markovchain)     # compute steady state distribution for markov chains
library(bayestestR)      # highest density intervals
library(mgcv)            # spline basis set up
library(viridis)         # nice colors
library(posterior)       # summaries and diagnostics for posterior samples
library(ggh4x)           # facet_nested

# tidyverse not working on RStudio because of problems with dbplyr
# library(tidyverse); theme_set(theme_bw())
library(ggplot2); theme_set(theme_bw())
library(tidyr)

# Functions ---------------------------------------------------------------

# Inverse multinomial link function
softmax <- function(x) exp(x) / sum(exp(x))

# (for probabilities)
normalize <- function(x) x / sum(x)

# Functions to summaryze
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

# compute delta index to compare distributions (1 - overlap)
delta_probs <- function(x) sum(abs(x)) / 2

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

# setwd() not needed: open srw_behaviour.Rproj so the working dir is the repo root # not necesary here, but necessary to run from raw R.
bdata <- read.csv(
  "data/srw_behaviour_data_1995-2025.csv",
  stringsAsFactors = FALSE, sep = ","
)
# data is in a local folder that I won't upload yet.

# Data description

# year = year of the observation
# month = month of the observation
# date = date as YYYY-mm-dd
# obs = observer id
# hour = hour of the interval, as hh:mm, character
# group = group of observations from the same animal-day-observer, 
#   used to define a follow
# golfo = gulf, as 
#   L, CH or Doradillo
# follow = id for continuous sequence of observations
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

# Time variable  ----------------------------------------------------------

# create variable minutes since 00:00 en every follow
bdata$min <- sapply(1:nrow(bdata), function(i) {
  sep <- strsplit(bdata$hour[i], ":") |> unlist() |> as.numeric()
  return(sep[1] * 60 + sep[2])
})

# order by year, follow and time
bdata <- bdata[order(bdata$year, bdata$follow, bdata$min), ]

# # check there are no jumps higher than 5 min within follows
# fols_temp <- unique(bdata$follow)
# sapply(length(fols_temp), function(f) {
#   return(diff(bdata$min[bdata$follow == fols_temp[f]]) |> unique())
# }) |> unique() # perfect

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

bdata$pam <- as.numeric(bdata$pam)
bdata$pac <- as.numeric(bdata$pac)

# Checks for missing data
unique(bdata$pac)
# "PETREL" and some NA
bdata$pac[bdata$pac == "PETREL"] <- 1

unique(as.character(bdata$pam)) # Also NA
which(is.na(bdata$pam))         # In the same rows
which(is.na(bdata$pac))         

# If pam/pac are NA, but fam fac have data, copy the same. If present, one 
# attack is the most likely.
bdata$pam[which(is.na(bdata$pam))] <- bdata$fam[which(is.na(bdata$pam))]
bdata$pac[which(is.na(bdata$pac))] <- bdata$fac[which(is.na(bdata$pac))]

# recreate fam anc fac by binarizing pam-pac
bdata$fam <- as.numeric(bdata$pam > 0)
bdata$fac <- as.numeric(bdata$pac > 0)

which(is.na(bdata$fam)) # only two obs with NA, in the same follow.
which(is.na(bdata$fam))

# Check follows with NA
rows_na <- which(is.na(bdata$pam) | is.na(bdata$pac))
fols_na <- bdata$follow[rows_na] |> unique()

# Another follow with NAs
fols_out <- c(fols_na, "VR-MS-KM-DT-RS34955ICH")
# VR-MS-KM-DT-RS34955ICH ## This follow has too many NA in behaviour because they didn't write UW/S. Delete it.
bdata <- bdata[!(bdata$follow %in% fols_out), ]

# mother-calf specific attack variables (not all are used in this model)
bdata$a <- as.numeric(bdata$pam > 0 | bdata$pac > 0)         # attack presence to mother or calf
bdata$ab <- as.numeric(bdata$pam > 0 & bdata$pac > 0)        # attack presence to both mother and calf
bdata$am_only <- as.numeric(bdata$pam > 0 & bdata$pac == 0)  # attack presence only to mother
bdata$am_marg <- as.numeric(bdata$pam > 0)                   # attack presence to mother independent of what happens to calf
bdata$ac_only <- as.numeric(bdata$pac > 0 & bdata$pam == 0)  # attack presence only to calf
bdata$ac_marg <- as.numeric(bdata$pac > 0)                   # attack presence to calf independent of what happens to mother
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
#table(bdata$behavm) |> plot
#(table(bdata$behavm) / sum(table(bdata$behavm)) ) |> plot

# define valid raw behaviors (I take out the undefined UWs or Ss)
behavs_m <- c("RUW", "RW", "RS", "SR", "GAL", "GP", "SR", "SRS",
              "TSUW", "TSS", #"TS", "TSY",
              "TMUW",
              "TMS", "TFS", "TFUW", #"TM", 
              "SA", "AS",
              "AUW", "UWA")

table(bdata$behavm[!bdata$behavm %in% behavs_m]) # abundance of invalid behaviours
# Most are "UW" (under water). I could also use others, such as TS, but 
# too few observations fall in this kind of category, so it won't change the 
# model fit, and including them would make the likelihood harder to code.

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
bdata$behavm_clust[bdata$behavm %in% c("RS", "SR", "GAL", "GP", "SR", "SRS")] <- "R_S"
bdata$behavm_clust[bdata$behavm %in% c("RUW", "RW")] <- "R_UW"

bdata$behavm_clust[bdata$behavm %in% c("TSS")] <- "ST_S"
bdata$behavm_clust[bdata$behavm %in% c("TSUW")] <- "ST_UW"

bdata$behavm_clust[bdata$behavm %in% c("TMS", "TFS")] <- "FT_S"
bdata$behavm_clust[bdata$behavm %in% c("TMUW", "TFUW")] <- "FT_UW"

bdata$behavm_clust[bdata$behavm %in% c("SA", "AS")] <- "IPA_S"
bdata$behavm_clust[bdata$behavm %in% c("AUW", "UWA")] <- "IPA_UW"

bdata$behavm_clust <- factor(bdata$behavm_clust, levels = behavs_m_clust)
# summary(bdata$behavm_clust)

# What are the high energy behaviours?
# (table(bdata$behavm[bdata$behavm_clust == "HE_UW"]) / 
#   sum( table(bdata$behavm[bdata$behavm_clust == "HE_UW"]) )) |> plot
# AUW       TFUW       TMUW 
# 0.04784191 0.13572543 0.81643266 
# Mostly travel under water

# (table(bdata$behavm[bdata$behavm_clust == "HE_S"]) / 
#     sum( table(bdata$behavm[bdata$behavm_clust == "HE_S"]) )) |> plot
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
# around 9 % is UW

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
(N_fol <- length(unique(bdata$follow))) # 2033

fol_agg <- aggregate(t ~ year + follow, FUN = length, data = bdata)
fol_agg <- fol_agg[order(fol_agg$year, fol_agg$follow), ]
bdata <- bdata[order(bdata$year, bdata$follow, bdata$t), ]

(mean(fol_agg$t *5)) # mean follow duration = 56.77 min
# plot(order(unique(bdata$follow))~ order(unique(fol_agg$follow))) # OK
# plot(table(fol_agg$t))

# identify each follow with an integer
bdata$follow_num <- factor(bdata$follow, levels = unique(bdata$follow)) |> as.numeric()

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

# add column indicating last observation by follow (useful to summarize attacks)
bdata$last_not <- TRUE
bdata$last_not[bdata$int_end] <- FALSE

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
year_length <- table(year_num_cons) |> unname()
year_end <- cumsum(year_length)
year_start <- year_end - year_length + 1


# Basis functions for year effects (splines) -------------------------------

# for transition matrix parameters:
bdata_year <- aggregate(fam ~ year, bdata, mean)
B <- 5 # number of basis functions including the intercept
b <- B - 1 # without intercept
bs <- smoothCon(s(year, k = B, bs = "cr"), diagonal.penalty = TRUE,
                absorb.cons = TRUE, data = bdata_year,
                scale.penalty = TRUE)[[1]]
spline_mfit <- bs$X
matplot(bdata_year$year, spline_mfit)

# for z dynamics parameters
Bz <- 5 # 3 # number of basis functions including the intercept
bsz <- smoothCon(s(year, k = Bz, bs = "cr"), diagonal.penalty = TRUE,
                 absorb.cons = TRUE, data = bdata_year,
                 scale.penalty = TRUE)[[1]]
spline_zfit <- cbind(bsz$X, rep(1, N_years))
matplot(bdata_year$year, spline_zfit)# , type = "l"

# columns are ordered from most wiggly basis function to constant.
# spline_zfit contains an intercept column, while 
# spline_mfit does not. It's because the intercepts for transition matrix 
# parameters are more variable than the spline. In that case, the spline affects
# full columns of the tmat, while the spline for z affects the whole z parameters.

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
        
      #(uw[i:last] |> sum) > 0      # there is at least one UW
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

# In stan I will only evaluate missing endings with UW present 
# (these are truly UW endings), but uw_ends also contains endings full of NA.
# So I subset them.

uw_ends$has_uw <- FALSE
for(i in 1:nrow(uw_ends)) {
  #i = 1
  sequ <- uw_ends$start[i] : uw_ends$end[i]
  uw_ends$has_uw[i] <- sum(bdata$uw[sequ]) > 0
}

uw_ends_subset <- uw_ends[uw_ends$has_uw == TRUE, ]


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

sum(bdata$uw) / sum(is.na(bdata$behavm_clust)) # 94 % of NA are UW
sum(bdata$uw) / nrow(bdata) # 6.35 % of obs are UW

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
  b = b, # number of basis functions for tmat spline without intercept
  Bz = Bz, # number of basis functions for zparams spline WITH intercept
  N_fol = N_fol, int_start = int_start, int_end = int_end,
  N_years = N_years, 
  year_start = year_start, year_end = year_end, 
  year_length = year_length, 
  year = bdata$year_num,
  
  N_jumps = N_jumps,
  jump_start = jumps$start,
  jump_end = jumps$end,
  jump_year = jumps$year,
  
  N_uw_ends = nrow(uw_ends_subset),
  end_start = uw_ends_subset$start,
  end_end = uw_ends_subset$end,
  end_year = uw_ends_subset$year,
  
  uw = bdata$uw,
  
  cons_id = cons_id,
  cons_id_prev = cons_id - 1,
  
  ind_matrix = ind_matrix,
  s_cols = s_cols,
  half_K = K / 2,

  # Variables for fitting
  y_mat = y_mat,
  y_vec = y_vec,
  # attack
  a = bdata$a,              # attack presence to mother or calf
  ab = bdata$ab,            # attack presence to both mother and calf
  ac_only = bdata$ac_only,  # attack presence only to mother
  ac_marg = bdata$ac_marg,  # attack presence to mother independent of what happens to calf
  nac = bdata$nac,          # number of attacks to calf
  # splines
  spline_mfit = spline_mfit, # without intercept
  spline_zfit = spline_zfit, # with intercep
  
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
  
  # z params priors
  prior_lambdac_sd = 2,
  
  # the following are priors for the coefficients of raw splines, and their order
  # is like columns in spline_zfit:
  # quadratic, linear, intercept.
  # dec refers to decay, incm to increase_m (they live in logit scale);
  # incbc control the deltas, which live in logit but the spline is parameterized
  # at its log, because delta must be >0.
  # deltas = exp(spline_raw_delta),
  # so they have different priors, which were decided in
  # /home/ivan/Insync/Whales/Behavior 4/restriccines para spline de increase_m ideas.R
   
  prior_dec_sds = rep(1.5, Bz),  #c(1, 1.5, 1.5),
  prior_incm_sds = rep(1.5, Bz), #c(1, 1.5, 1.5),
  prior_incbc_sds = rep(1, Bz),  #c(0.5, 1, 1),
  
  # Limits for increase and decay (mothers):
  L_dec = 0.1745958,
  L_inc = 0.3187079 # better 
)
L_dec = 0.1745958 # reaches 0.1 after 1 h without attacks.
#L_inc = 0.4376587 # reaches 0.9 after 4 consecutive intervals with attack.
L_inc = 0.3187079 # reaches 0.9 after 6 consecutive intervals with attack.

# makes sense to let z increase little only with attacks to mothers because
# in years with many attacks to calves it can reach higher values. 

str(sdata)

# N_fol # 2033
# range(fol_agg$t) # 5 to 53


# Model fit ---------------------------------------------------------------

# Prior simulations for the model in
# <behaviour_mothers_analysis_prior_simulations.R>

# compile and sample
bm_code <- stan_model(file = "behaviour/behaviour_mothers_model.stan", verbose = TRUE)
warmup <- 500
N_samples <- 1500
bm <- sampling(
  bm_code, seed = 5454, refresh = 1,
  chains = 8, cores = 8, iter = warmup + N_samples, warmup = warmup,
  #chains = 1, cores = 1, iter = 4, # used to test the stan code and data
  data = sdata,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  pars = c("sd_spline_tmat", "sd_years",
           "decay",
           "increase_m",
           "increase_b", "increase_b3",
           "increase_c", "increase_c3",
           "kappa", "delta_b", "delta_c", "lambda_c",
           "alpha", "beta",                        # splines + random effects
           "alpha_intercept", "beta_intercept",    # intercepts
           "alpha_raneff",                         # random effects
           "alpha_spline_coef", "beta_spline_coef",# tmat splines coefficients
           "zpar_spline_coef",
           "z_0") # sample without saving z because it uses too much memory
                  # that fills up the 16 GB RAM when I try to extract it.
                  # Compute z afterwards in R, for a matrix uses much less memory.
)
bms <- summary(bm)[[1]]
head(bms)
summary(bms[, "n_eff"])
summary(bms[, "Rhat"])
saveRDS(bm, "models/behaviour/behaviour_mothers_model_samples.rds")
saveRDS(sdata, "models/behaviour/behaviour_mothers_model_data.rds")
# done

# it took 21696.3/3600 = 6.02 h

# Load model
bm <- readRDS("models/behaviour/behaviour_mothers_model_samples.rds")

# Convergence diagnostics and effective sample size -----------------------

check_hmc_diagnostics(bm)
get_bfmi(bm)
get_low_bfmi_chains(bm)
# Everything is OK.

# summaries for model using posterior package

# turn samples into df
bmdf <- as_draws_df(bm)

# exclude some parameters
out1 <- grep("_raneff", colnames(bmdf))
out2 <- which(colSums(bmdf) == 0) |> unname()

bmdf_subset <- bmdf[, -c(out1, out2)]

bm_summary <- summarise_draws(bmdf_subset)
# it takes long, save it:
saveRDS(bm_summary, "behaviour/files/behaviour_mothers_model_summary.rds")

# Which are the worst-sampled parameters?
summary(bm_summary[, c("rhat", "ess_bulk", "ess_tail")])
head(bm_summary[order(bm_summary$ess_bulk), ], n = 20)
# min ess is 1408 with 1000 * 8 samples.

(get_elapsed_time(bm) |> rowSums() |> max()) / 3600 # 6.026767 h running

# Computing z -------------------------------------------------------------

# The latent variable z was not saved in the stan object because it occupied
# too much memory. However, computing it and saving it as a matrix is not so 
# heavy. This will be useful to estimate missing behaviours and might be used
# to compute likelihoods and stuff like that.

z0_samples <- as.matrix(bm, pars = "z_0") # [iters, follows]
N_samples <- nrow(z0_samples)
# get z parameters
kappa = as.matrix(bm, pars = c("kappa"))
delta_b = as.matrix(bm, pars = c("delta_b"))
delta_c = as.matrix(bm, pars = c("delta_c"))
lambda_c = as.matrix(bm, pars = c("lambda_c"))
decay = as.matrix(bm, pars = c("decay"))

z_samples <- matrix(NA, sdata$N, N_samples) # [obs, iters]

# Compute z_samples at first interval in each follow
#N_samples <- 10#nrow(z0_samples)

# IMPORTANT: R does not support (+) or (-) at the line beginning.

for(i in 1:N_samples) {
  print(i)
  for(f in 1:N_fol) {
    z_samples[int_start[f], i] <- z0_samples[i, f] +
      # increase
      bdata$a[int_start[f]] * (1 - z0_samples[i, f]) * plogis(
    
        # intercept: reference is attack presence to mother
        kappa[i, sdata$year[int_start[f]]] +                           
      
        # if both are attacked, the intercept increases:
        delta_b[i, sdata$year[int_start[f]]] * bdata$ab[int_start[f]] -
    
        # if only calf is attacked, the intercept decreases:
        delta_c[i, sdata$year[int_start[f]]] * bdata$ac_only[int_start[f]] +                           
    
        # if there is more than 1 attack to calf, the inc proportion rises
        lambda_c[i] * (bdata$nac[int_start[f]] - 1) * bdata$ac_marg[int_start[f]]  
    ) -
    # decay
    (1 - bdata$a[int_start[f]]) * decay[i, sdata$year[int_start[f]]] * z0_samples[i, f]
  }
}

# remaining z_samples's
for(i in 1:N_samples) {
  print(i)
  
  for(f in 1:N_fol) {
    
    for(t in (int_start[f] + 1):int_end[f]) {
      
      z_samples[t, i] <- z_samples[t-1, i] +
        # increase
        bdata$a[t] * (1 - z_samples[t-1, i]) * plogis(
      
          # intercept: reference is attack presence to mother
          kappa[i, sdata$year[t]] +                          
      
          # if both are attacked, the intercept increases:
          delta_b[i, sdata$year[t]] * bdata$ab[t] -
      
          # if only calf is attacked, the intercept decreases:
          delta_c[i, sdata$year[t]] * bdata$ac_only[t] +                         
      
          # if there is more than 1 attack to calf, the inc proportion rises
          lambda_c[i] * (bdata$nac[t] - 1) * bdata$ac_marg[t]  
        ) -
        # decay
        (1 - bdata$a[t]) * decay[i, sdata$year[t]] * z_samples[t-1, i]
    }
  }
}
dimnames(z_samples) <- list("t" = 1:nrow(bdata), "iter" = 1:N_samples)
saveRDS(z_samples, "models/behaviour/behaviour_mothers_model_samples_z.rds")
z_samples <- readRDS("models/behaviour/behaviour_mothers_model_samples_z.rds")
View(cbind(bdata[, c("follow", "t", "a")], z_samples[, 1:5]))


# Disturbance score in 1995 ----------------------------------------------
# remove the last observation by follow, as attacks were not recorded there.

z95_samples <- z_samples[bdata$year == 1995 & bdata$last_not, ]
z95_means <- colMeans(z95_samples) 
hdmean(z95_means, name = "z") |> round(digits = 3)
# z_lower     z_mean    z_upper 
#   0.075     0.100      0.132 
(bdata$year == 1995) |> sum # 1603 datos.

bdata$fam[bdata$year == 1995 & bdata$last_not] |> mean # 0.108 
# muy similar a lo que dice Rowntree:
 
# "We recorded behavioral activity at the beginning of 3,499
# 5-min intervals that were distributed among the observation sites 
# as follows: Cliff Hut 1,824; Fracaso 759; Golfo Nuevo 916. [...]
# Gull attacks occurred in 17%
# of the 3,499 5-min intervals, but the frequency was much higher at Fracaso
# (36%) than at the nearby Cliff Hut (12%) or at Golfo Nuevo (12%)"

# Los datos sobre pares madre-cría son menos. Tenían 2783, que dividieron en 
# no attack y attack. En la primera categoría eliminaron los 12 primeros intervalos
# de cada follow, que ya de por sí no tenían ataques. 
# De esta forma, Los N con ataque eran
enes <- c(248, 221, 140, 110, 484, 68, 64, 52, 42, 34, 26, 19, 18, 16, 16, 13, 10, 8, 7, 49)
cumsum(enes)[enes == 34] # aprox hasta ahí son 50 min, que llega a 1463,
# siendo menos que lo que dicen: 1246 (habrá habido NAs?)
# Dicen que en total tenían 280 (undisturbed) + 1246 (disturbed up to 50 min)
# 280 + 1246 = 1526

# Sobre el avg disturbance score:
# Assuming a disturbance half-life of 7.5 *in, the exponential decay model
# estimates an average disturbance score of 0.24 over all follows of lone mother-
# calf pairs. This figure implies that a mother-calf pair typically spent the equiv-
# alent of 24% of its daylight hours (or 2.9 h of a 12-h day) in a state of gull-
# induced disturbance. Lower (5.1 min) and higher (9.9 min) estimates of the
# half-life of disturbance give average disturbance scores of 0.19 and 0.27, re-
# spectively.

# Si simulo con la misma tasa de ataque, obtengo lo mismo? Debería obtener algo
# más alto suponiendo ataques independientes. 

# Usando la attack rate de modres-crias en todo el data set (0.19)

# y = exp(-decay * t)
# 0.5 = exp(-decay * 7.5) # half-life is 7.5 min
# log(0.5) = -decay * 7.5
# -log(0.5) / 7.5 = decay
decay_fit <- -log(0.5) / 7.5

dist_score <- function(prob = 0.19, decay = decay_fit, nsim = 1000, replicates = 1000) {
  at <- matrix(rbinom(nsim * replicates, size = 1, prob = prob), nsim, replicates)
  score <- matrix(0, nsim, replicates)
  
  score[1, ] <- at[1, ]
  for(t in 2:nsim) {
    # for whales not attacked at t:
    score[t, at[t, ] == 0] <- score[t-1, at[t, ] == 0] * exp(-decay * 5) # el coef de decay es similar al estimado por mí para z!
    # for whales attacked at t:
    score[t, at[t, ] == 1] <- 1
  }
  return(mean(score))
}

dist_score(prob = 0.19, decay = decay_fit) # 0.387
dist_score(prob = 0.12, decay = decay_fit) # 0.27
# bad match. Será por usar ataque independiente? Eso debería sobreestimar el score.
# Uso lo observado y comparon con la misma tasa de ataque

# uso ataque observado en 1995
at_95 <- bdata$fa[bdata$year == 1995] |> as.numeric()

score_hand <- numeric(length(at_use))
score_hand[1] <- at_use[1]
for(t in 2:length(at_use)) {
  score_hand[t] <- at_use[t-1] * exp(-decay * 5) * (1 - at_use[t]) +
                   at_use[t]
}
mean(score_hand); range(score_hand)
curve(exp(-decay * x), from = 0, to = 60)
# con este criterio, usando datos del 95, obtenemos un score medio de 
# 0.17.


# Estimation of unobserved behaviours: JUMPS ------------------------------

bdata$behavm_hat <- NA
year_cols <- rep(1:N_years, K*K)

# Extract posterior samples

alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")

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
 
  # loop over posterior samples and sequences
  
  for(i in 1:N_samples) { 
  #for(i in 1:10) { 
    #i = 628

    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)

    z <- z_samples[(jumps$start[j] - 1) : (jumps$end[j] - 1), i]
    # take z from start-1 to end-1 because z[t-1] defines the transition matrix
    # from t-1 to t
    
    # Create transition matrices
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      #int = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) |> t()
    }
    
    # Transition matrices must be modified if UW is in the missing sequence:
    # if a tmat leads to an UW behaviour, the surface columns (s_cols) must be set to 0, 
    # and the (modified) rows have to be normalized. The loop goes up to 
    # N_trans - 1 because the last transition leads always to a fully observed behaviour, 
    # so it can't be UW.
    for(tran in 1:(N_trans - 1)) { 
      if(missing_uw[tran] == 1) {
        tmats[, s_cols, tran] <- 0
        tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) |> t()
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
      
      seqs_probs[s, i] <- trans_probs |> log() |> sum() |> exp()
    }# 
    
  } # end loop over posterior samples
  # it takes a lot of time by jump
  
  # compute probability of every sequence marginal to the posterior samples
  #str(seqs_probs)
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  #plot(seqs_probs_norm)
  #sum(seqs_probs_norm)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) |> as.data.frame()
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] |> as.numeric()
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
    
  # Fill most probable behaviours
  bdata$behavm_hat[jumps$start[j] : (jumps$end[j] - 1)] <- most_prob_seq_text 
  
  # Fill seqs result
  jump_seq_hat[[j]] <- seqs_result
}

# a long jump:
# jump_seq_hat[[58]]
# check whether it worked
# follows with jumps
fwj <- bdata$follow[jumps$start] |> unique()
View(bdata[bdata$follow %in% fwj,
           c("follow", "follow_num", "t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "models/behaviour/behaviour_mothers_data with imputed jumps.rds")


# Estimation of unobserved behaviours: ENDS -------------------------------

# Missing behaviours for which there are no posterior observed behaviours.

# list containing all the possible sequences in every ending with their marginal 
# posterior probability
ending_seq_hat <- vector(mode = "list", length = N_uw_ends)

# Compute most likely sequences
for(j in 1:N_uw_ends) { 
  #j = 32
  print(paste("ending", j, "/", N_uw_ends))
  
  N_miss <- uw_ends$start[j]:uw_ends$end[j] |> length()
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

  # loop over posterior samples and sequences
  for(i in 1:N_samples) {
  #for(i in 1:10) {  
    #i = 1
  
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    
    z <- z_samples[(uw_ends$start[j] - 1) : (uw_ends$end[j] - 1), i]
    # take z from start-1 to end-1 because z[t-1] defines the transition matrix
    # from t-1 to t

    # make transition matrices array
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      # tran = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) |> t()
      
      # Transitions to UW have zero probability of surface behaviours.
      if(missing_uw[tran] == 1) {
        tmats[, s_cols, tran] <- 0
        tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) |> t()
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
      
      seqs_probs[s, i] <- trans_probs |> log() |> sum() |> exp()
    }# 
  } # end loop over posterior samples
  # it takes a lot of time by jump
  
  # compute probability of every ending sequence marginal to the posterior samples
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  #plot(seqs_probs_norm)
  #sum(seqs_probs_norm)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) |> as.data.frame()
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] |> as.numeric()
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
  
  # Fill most probable behaviours
  bdata$behavm_hat[uw_ends$start[j] : (uw_ends$end[j])] <- most_prob_seq_text 
  
  # Fill seqs result
  ending_seq_hat[[j]] <- seqs_result
  
}

# check whether it worked
# follows with uw_endings
fwuwe <- bdata$follow[uw_ends$start] |> unique()
View(bdata[bdata$follow %in% fwuwe,
           c("follow", "follow_num", "t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "models/behaviour/behaviour_mothers_data with imputed jumps and endings.rds")


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
  for(i in 1:N_samples) { 
  #for(i in 1:10) {
    #i = 1
    
    alpha <- matrix(alpha_samples[i, which(year_cols == y)],
                    K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)],
                   K, K)
    
    z <- z_samples[(startlosts$start[j]) : (startlosts$end[j] - 1), i]
    # We need z from startlosts$start[j] to startlosts$end[j]-1 because z[t-1] 
    # defines the transition matrix from t-1 to t.
    
    # get transition matrices array
    # these transitions go from the 1st missing behaviour (t = 1 in the follow)
    # up to the first observed behaviour.
    tmats <- array(NA, c(K, K, N_trans))
    for(tran in 1:N_trans) {
      #tran = 1
      tmats[, , tran] <- apply(alpha + beta * z[tran], 1, softmax) |> t()
    }
    
    # Transitions to UW have zero probability of surface behaviours.
    # Notice this does not affect the last transition, because it leads to an
    # observed behaviour. So, if there is only 1 missing start, it's omitted.
    if(N_trans > 1) {
      for(tran in 1:(N_trans - 1)) {
        if(missing_uw[tran + 1] == 1) {
          tmats[, s_cols, tran] <- 0
          tmats[, , tran] <- apply(tmats[, , tran], 1, normalize) |> t()
        }
      }
    }
    
    # Compute probability of first unobserved behavior unconditional to previous
    # behaviour (steady states distribution for z0)
    z0 <- z0_samples[i, f]
    
    tmat_0 <- apply(alpha + beta * z0, 1, softmax) |> t()
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
      seqs_probs[s, i] <- trans_probs_wide |> log() |> sum() |> exp()
    }# 
  } # end loop over posterior samples
  
  # compute probability of every starting sequence marginal to the posterior samples
  seqs_probs_means <- apply(seqs_probs, 1, mean, na.rm = TRUE)
  seqs_probs_norm <- seqs_probs_means / sum(seqs_probs_means)
  
  seqs_result <- cbind(seqs, p = seqs_probs_norm) |> as.data.frame()
  seqs_result <- seqs_result[order(seqs_result$p, decreasing = TRUE), ]
  rownames(seqs_result) <- 1:N_seqs
  
  most_prob_seq <- seqs_result[1, 1:N_miss] |> as.numeric()
  most_prob_seq_text <- behavs_m_clust[most_prob_seq]
  
  # Fill most probable behaviours
  bdata$behavm_hat[startlosts$start[j] : (startlosts$end[j] - 1)] <- most_prob_seq_text 
  
  # Fill seqs result
  startlost_seq_hat[[j]] <- seqs_result
  
}

# check whether it worked
# follows with lost starts
fwsl <- bdata$follow[startlosts$start] |> unique()
View(bdata[bdata$follow %in% fwsl,
           c("follow", "follow_num","t", "behavm_clust", "behavm_hat", "missing_obs")])
# save
saveRDS(bdata, "models/behaviour/behaviour_mothers_data with imputed jumps, endings and starts.rds")



# Merge and export observed and estimated behaviours ------------------------

# This is necessary to compute the observed behaviours proportion by year.
bdata$behavm_full <- bdata$behavm_clust
bdata$behavm_full[is.na(bdata$behavm_clust)] <- bdata$behavm_hat[is.na(bdata$behavm_clust)]

# check there are not missing behaviours
which(is.na(bdata$behavm_full)) # perfect

# View(bdata[bdata$follow %in% bdata$follow[is.na(bdata$behavm_full)],
#            c("follow", "t", "behavm_clust", "uw", "behavm_hat", "missing_obs")])

# save data with imputed behaviours
write.csv(bdata, "models/behaviour/behaviour_mothers_data with imputed behaviours.csv")
saveRDS(bdata, "models/behaviour/behaviour_mothers_data with imputed behaviours.rds")


hat_list <- list(jumps = jump_seq_hat, endings = ending_seq_hat, starts = startlost_seq_hat)
saveRDS(hat_list, "models/behaviour/behaviour_mothers_estimated missing sequences.rds")
# Estimation of unobserved behaviours finished.

hat_list <- readRDS("models/behaviour/behaviour_mothers_estimated missing sequences.rds")
str(hat_list)
hat_list$starts[[1]]

bdata2 <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")

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

increase <- as.matrix(bm, pars = "increase_m")
decay <- as.matrix(bm, pars = "decay")

#str(increase) # 15 columnas.
# period_year <- c(1, rep(2, N_years - 1))

# here I will align z with the behavior probabilities in the same interval.
# So, now it's z[t] the one determining the transition probs from t-1 to t.
# The model was written differently to get a differente alignment to
# perform matrix multiplication.
# Here z depends on the previos attacks (which are constant in these two 
# scenarios).

### DON'T RUN, HEAVY LOOP ###

pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
for(i in 1:N_samples) {
  #i = 1

  # loop over years
  for(y in 1:N_years) {
    #y = 1
    alpha <- matrix(alpha_samples[i, which(year_cols == y)], K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)], K, K)

    inc <- increase[i, y]
    dec <- decay[i, y]

    # Attack scenario ---------------------------

    z <- 0 # starts calm
    pred_attack[1, "z", y, i] <- z
    tmat <- apply(alpha + beta * z, 1, softmax) |> t()
    mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat,
              name = "mc")
    ss <- steadyStates(mc)
    pred_attack[1, 2:(K+1), y, i] <- ss


    for(t in 2:nrow(pred_attack)) {
      #t = 2
      pred_attack[t, "z", y, i] <- pred_attack[t-1, "z", y, i] +
        inc * (1 - pred_attack[t-1, "z", y, i])
      tmat <- apply(alpha + beta * pred_attack[t, "z", y, i], 1, softmax) |> t()

      pred_attack[t, 2:(K+1), y, i] <- pred_attack[t-1, 2:(K+1), y, i] %*% tmat
    }

    #pred_attack[, , i]

    # No attack scenario --------------------------

    z <- 1 # starts upset
    pred_no_attack[1, "z", y, i] <- z
    tmat <- apply(alpha + beta * z, 1, softmax) |> t()
    mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat,
              name = "mc")
    ss <- steadyStates(mc)
    pred_no_attack[1, 2:(K+1), y, i] <- ss

    for(t in 2:nrow(pred_no_attack)) {
      #t = 2
      pred_no_attack[t, "z", y, i] <- pred_no_attack[t-1, "z", y, i] -
        dec * pred_no_attack[t-1, "z", y, i]
      tmat <- apply(alpha + beta * pred_no_attack[t, "z", y, i], 1, softmax) |> t()

      pred_no_attack[t, 2:(K+1), y, i] <- pred_no_attack[t-1, 2:(K+1), y, i] %*% tmat
    }

    #pred_no_attack[, , y, i]
  }
  setTxtProgressBar(pb, i)
}

pred_samples <- list(attack = pred_attack,
                     no_attack = pred_no_attack)
saveRDS(pred_samples, "models/behaviour/behaviour_mothers_predictions_attack scenario samples.rds")
### end heavy loop

pred_samples <- readRDS("models/behaviour/behaviour_mothers_predictions_attack scenario samples.rds")

at_temp <- apply(pred_samples$attack, 1:3, hdmean, name = "prob") |> as.data.frame.table()
at_no_temp <- apply(pred_samples$no_attack, 1:3, hdmean, name = "prob") |> as.data.frame.table()

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

predictions <- separate_behav(predictions)

# guardo tabla
saveRDS(predictions, "behaviour/files/behaviour_mothers_predictions_attack scenarios table.rds")
predictions <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table.rds")

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
           linetype = exposure)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  #facet_grid(rows = vars(year), cols = vars(behav_raw)) +
  facet_nested(cols = vars(behav_raw, exposure), rows = vars(year)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behavior probability") +
  xlab("Time (min)") #+

ggsave("behaviour/figures/behaviour_mothers_predictions_attack scenarios all years.png", 
       height = 35, width = 20, units = "cm")


filt2 <- with(predictions, 
  min <= 70 & 
  behav != "z" &
  year %in% c("1995", "2004", "2009", "2015") 
)

ggplot(predictions[filt2, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = exposure)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_nested(cols = vars(behav_raw, exposure), rows = vars(year)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behavior probability") +
  xlab("Time (min)") #+
  #ggtitle("Behaviour in response to kelp gull micropredation") #+

ggsave("behaviour/figures/behaviour_mothers_predictions_attack scenarios few years.png", 
       height = 15, width = 23, units = "cm")


# Explore delay between z and probabilities
predictions$exposure <- as.character(predictions$exposure)
predictions$exposure[predictions$behav_raw == "z"] <- "both"
filt_z <- with(predictions, 
                min <= 70 & 
                year %in% c("1995") & 
                scenario == "After attack" & 
                behav_raw %in% c("Rest", "z") & 
                exposure != "Under water"
)

predictions$var_type <- "behav"
predictions$var_type[predictions$behav == "z"] <- "z"

ggplot(predictions[filt_z, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  #facet_wrap(cols = vars(behav_raw, exposure), rows = vars(year)) +
  facet_wrap(vars(behav_raw), nrow = 2) +
  #scale_color_manual(values = colors) + 
  #scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behavior probability") +
  xlab("Time (min)")


# The same plot but averaging years (2013 to 2018)
years_calves <- as.character(c(2013, 2015:2018))

str(pred_samples$attack[, , years_calves, ])

tmp1 <- apply(pred_samples$attack[, , years_calves, ], c(1, 2, 4), mean)
at_temp_avg <- apply(tmp1, 1:2, hdmean, name = "prob") |> as.data.frame.table()
at_temp_avg$scenario <- "Attack"

tmp2 <- apply(pred_samples$no_attack[, , years_calves, ], c(1, 2, 4), mean)
at_no_temp_avg <- apply(tmp2, 1:2, hdmean, name = "prob") |> as.data.frame.table()
at_no_temp_avg$scenario <- "After attack"

ptemp <- rbind(at_temp_avg, at_no_temp_avg)

predictions_avg <- pivot_wider(ptemp, names_from = "Var1", values_from = "Freq")

predictions_avg$time <- as.numeric(as.character(predictions_avg$time))
predictions_avg$min <- predictions_avg$time * 5

predictions_avg <- predictions_avg |> separate_behav()

# save
saveRDS(predictions_avg, "behaviour/files/behaviour_mothers_predictions_attack scenarios table_avg.rds")
predictions_avg <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table_avg.rds")

# Plots
colors <- c("#440154", "#35b779")

filt1 <- with(predictions_avg, 
              min <= 70 & 
                behav != "z"
)

ggplot(predictions_avg[filt1, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = exposure)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  #facet_grid(rows = vars(year), cols = vars(behav_raw)) +
  facet_nested(cols = vars(behav_raw, exposure)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behavior probability") +
  xlab("Time (min)") 

# later merge this plot with the calves' one, choosing 3 years


# Disturbed, undisturbed and observed by year -------------------------------


predictions2 <- predictions[predictions$time == 0, ]
# The starting condition in Attack is undisturbed, and in After attacks it's 
# disturbed
predictions2$scenario <- factor(predictions2$scenario, 
                                levels = c(levels(predictions$scenario), "Observed"),
                                labels = c("Undisturbed", "Disturbed", "Observed"))

predictions2$exposure <- factor(predictions2$exposure, 
                                levels = c("Surface", "Under water"))

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
#bdata <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")

# ATENTION, USE behavm_full HERE
bdata <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")
y_mat_full <- model.matrix.lm(pam ~ behavm_full - 1, data = bdata, na.action = "na.pass")
#y_mat_full <- model.matrix.lm(pam ~ behavm_clust - 1, data = bdata, na.action = "na.pass")

y_agg <- aggregate(y_mat_full ~ bdata$year, FUN = mean, na.rm = TRUE)
colnames(y_agg) <- c("year", behavs_m_clust)
y_agg$year_num <- as.numeric(as.character(y_agg$year))
y_agg <- pivot_longer(y_agg, which(colnames(y_agg) %in% behavs_m_clust), 
                      values_to = "prob_mean", names_to = "behav")
y_agg$scenario <- "Observed"
y_agg <- separate_behav(y_agg)
y_agg$prob_lower <- NA
y_agg$prob_upper <- NA
y_agg$period <- factor("new", levels = c("old", "new"))
y_agg[y_agg$year == 1995, "period"] <- "old"

# filter
filt3 <- with(predictions2, behav != "z")
colors <- c("#440154", "#35b779", "black")


# subset useful columns and merge.
p <- predictions2[filt3, c("year_num", "period", "behav_raw", "exposure", 
                           "scenario", "prob_mean", "prob_lower", "prob_upper")]
w <- y_agg[, c("year_num", "period", "behav_raw", "exposure", 
               "scenario", "prob_mean", "prob_lower", "prob_upper")]

ppp <- rbind(p, w)
filt4 <- !(ppp$behav_raw == "In place activity" & ppp$exposure == "Under water")


saveRDS(list(plot_data = ppp, ghosts = ghosts),
        "behaviour/files/behaviour_mothers_predictions disturbed and undisturbed by year.rds")


ggplot(ppp,#[filt4, ], 
       aes(x = year_num, 
           y = prob_mean, ymin = prob_lower, ymax = prob_upper, 
           colour = scenario, shape = scenario)) + 
  geom_point(exposure = position_dodge2(width = 0.3), size = 3) +
  geom_linerange(exposure = position_dodge2(width = 0.3), show.legend = FALSE) +
  geom_point(data = ghosts, mapping = aes(x = year_num, y = prob_lower), alpha = 0) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(period), 
               #scale = "free", space = "free") +
               scales = "free", space = "free_x"
  )+
  scale_color_manual(values = colors, name = "Condition") +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(1, 1, 20)) +
  guides(shape = "none") +
  ylab("Behavior probability") +
  xlab("Year") +
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1995, 2005, 2010, 2015, 2020, 2025)) 

ggsave("behaviour/figures/behaviour_mothers_ts_free_y.png", height = 25, width = 20, unit = "cm")


# Kelp gull effects by year -----------------------------------------------

undisturbed <- pred_samples$attack[1, 2:(K+1), , ]
disturbed <- pred_samples$no_attack[1, 2:(K+1), , ]

y_agg2 <- aggregate(y_mat_full ~ bdata$year, FUN = mean, na.rm = TRUE)[, -1] |> as.matrix() |> t()

eff_pot_samples <- disturbed - undisturbed
eff_obs_samples <- eff_pot_samples
for(i in 1:dim(eff_obs_samples)[3]) {
  eff_obs_samples[, , i] <- y_agg2 - undisturbed[, , i]
}

eff_pot <- apply(eff_pot_samples, 1:2, hdmean, name = "eff") |> as.data.frame.table()
eff_obs <- apply(eff_obs_samples, 1:2, hdmean, name = "eff") |> as.data.frame.table()
eff_pot$effect_type <- "Potential"
eff_obs$effect_type <- "Observed"
effs <- rbind(eff_pot, eff_obs)

effs <- pivot_wider(effs, names_from = "Var1", values_from = "Freq")
effs <- separate_behav(effs)
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
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(period), 
               scales = "free", space = "free") +
  scale_color_manual(values = colors_eff, name = "Effect type") +
  ylab("Behavior probability difference (x - undisturbed)") +
  xlab("Year") +
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1995, 2005, 2008, 2012, 2016)) 

ggsave("behaviour/figures/behaviour_mothers_ts effects.png", height = 20, width = 20, unit = "cm")


# Total effects (summarizing behaviours) -----------------------------------

# Let's call a_obs the observed attack scenario, and a_0, the attack absence 
# scenario, which results in undisturbed behaviour distribution (behav_undisturbed). 
# To measure the net effect of attacks (a_obs) on behaviour, we compare a
# behaviour distribution under a_obs with the distribution for a_0 (undisturbed).
# One approach is to assume that behav_obs, the observed proportion of behaviours
# by year, is the distribution under a_obs. However, the difference 
# |behav_obs - behav_undisturbed| 
# is caused by the difference a_obs - a_0 and by the model error, because the 
# model does not predict perfectly the observed data. This would probably inflate
# the attack effects. A better approach, lacking model error is compared two 
# predicted distributions: behav_obspred - behav_undisturbed.
# behav_obspred is the behaviour distribution by year predicted by the model using
# the estimated z values in each interval. We initialize every follow in the first
# observed behaviour, but simulate the others in a recursive way. 
# We'll compute both effects, but first, we need to simulate behav_obspred.

# Predicted behaviours for estimated z values (behav_obspred) -------------

# Compute the predicted distribution under observed attack (z) scenario.
# The predicted distribution will be that one, starting from observed behaviours
# in time t = 1 and simulating the following given the estimated z values. 

set.seed(342)
str(z_samples)
pred_dist_samples <- matrix(NA, nrow(z_samples), ncol(z_samples))
# get observed initial value
pred_dist_samples[int_start, ] <- as.numeric(bdata$behavm_full[int_start])


for(s in 1:ncol(z_samples)) {   # posterior samples
  print(s)
  for(f in 1:N_fol) {           # follows
    for(i in (int_start[f] + 1):int_end[f]) {    # intervals
      # get parameters
      # f = 1; s = 1; i = int_start[f] + 1 ## test
      y <- bdata$year_num[i]
      row <- pred_dist_samples[i-1, s] # get previous behaviour
      alpha <- matrix(alpha_samples[s, which(year_cols == y)], K, K)
      beta <- matrix(beta_samples[s, which(year_cols == y)], K, K)
      probs <- softmax(alpha[row, ] + beta[row, ] * z_samples[i-1, s])
      # z[t] affects transition from y[t] to y[t+1])
      pred_dist_samples[i, s] <- sample(1:K, size = 1, prob = probs)
    }
  }
} 
saveRDS(pred_dist_samples, "models/behaviour/behaviour_mothers_model_samples_predict_observed_attack.rds")
pred_dist_samples <- readRDS("models/behaviour/behaviour_mothers_model_samples_predict_observed_attack.rds")
# summarize by year
# make the same structure as prediction_samples (compare later)
pred_dist_samples_year <- array(NA, dim = c(N_years, K, N_samples),
                                dimnames = list(year = years_unique, 
                                                behav = behavs_m_clust,
                                                sample = 1:N_samples))
for(s in 1:N_samples) {
  print(s)
  for(y in 1:N_years) {
    # "data" summary
    #y = 1; s = 1
    nums <- pred_dist_samples[bdata$year_num == y, s]
    ttt <- model.matrix(~ factor(nums, levels = as.character(1:K)) - 1)
    pred_dist_samples_year[y, , s] <- colMeans(ttt) |> unname()
  }
}
saveRDS(pred_dist_samples_year, "models/behaviour/behaviour_mothers_model_samples_predict_observed_attack_years.rds")
pred_dist_samples_year <- readRDS("models/behaviour/behaviour_mothers_model_samples_predict_observed_attack_years.rds")


pred_dist <- apply(pred_dist_samples_year, 1:2, hdmean, name = "pred") 
# ojo con el hdmean para cantidades discretas. (muchos ceros)
pred_dist <- as.data.frame.table(pred_dist)
pred_dist_wide <- pivot_wider(pred_dist, names_from = "Var1", values_from = "Freq")

w2 <- y_agg[, c("year_num", "behav", "prob_mean")]
names(w2)[1:2] <- names(pred_dist_wide)[1:2]
w2$year; 
pred_dist_wide$year <- as.numeric(as.character(pred_dist_wide$year))



saveRDS(pred_dist_wide, "behaviour/files/behaviour_mothers_model_pred dist wide.rds")
pred_dist_wide <- readRDS("behaviour/files/behaviour_mothers_model_pred dist wide.rds")

ggplot(pred_dist_wide, aes(x = year, y = pred_mean, ymin = pred_lower, ymax = pred_upper)) +
  geom_point() + 
  geom_linerange() +
  #facet_wrap(vars(behav), scales = "free") +
  facet_wrap(vars(behav)) +
  geom_point(data = w2, mapping = aes(x = year, y = prob_mean), shape = 1, col = 2,
             inherit.aes = F) 

# The predicted distribution (behav_obspred) matches well the observed (behav_obs)


# Total effects using "obspred" proportions ----------------------------------

# A few inputs
pred_samples <- readRDS("models/behaviour/behaviour_mothers_predictions_attack scenario samples.rds")
undisturbed <- pred_samples$attack[1, 2:(K+1), , ]
undisturbed_95 <- pred_samples$attack[1, 2:(K+1), , ]
for(y in 1:N_years) undisturbed_95[, y, ] <- undisturbed[, 1, ]
disturbed <- pred_samples$no_attack[1, 2:(K+1), , ]
y_agg2 <- aggregate(y_mat_full ~ bdata$year, FUN = mean, na.rm = TRUE)[, -1] |> as.matrix() |> t()

# order prediction samples as y_agg2
obspred_samples <- aperm(pred_dist_samples_year, perm = c(2, 1, 3))

# large array to fill
eff_samples <- abind::abind("Potential" = undisturbed,
                            "Short-term" = undisturbed,
                            "Long-term" = undisturbed,
                            "Total" = undisturbed,
                            along = length(dim(undisturbed)) + 1)

years_unique <- unique(bdata$year)
N_samples <- dim(z_samples)[2]

dimnames(eff_samples) <- list(
  behav = behavs_m_clust,
  year = years_unique, 
  sample = 1:N_samples,
  effect_type = c("Potential", "Short-term", "Long-term", "Total")
)

eff_samples[, , , "Potential"] <- disturbed - undisturbed
eff_samples[, , , "Short-term"] <- obspred_samples - undisturbed
eff_samples[, , , "Long-term"] <- undisturbed - undisturbed_95
eff_samples[, , , "Total"] <- obspred_samples - undisturbed_95

# Get delta index
delta_samples <- apply(eff_samples, 2:4, delta_probs)

# Summarize
delta_summ <- apply(delta_samples, c(1, 3), hdmean, name = "d") |> as.data.frame.table()
delta_wide <- pivot_wider(delta_summ, names_from = "Var1", values_from = "Freq")
delta_wide$year <- as.numeric(as.character(delta_wide$year))

# order effects labels
delta_wide$eff_type_term <- delta_wide$effect_type
delta_wide$eff_type_term <- plyr::revalue(delta_wide$effect_type, replace = c(
  "Potential" = "Short-term"
))
with(delta_wide, table(effect_type, eff_type_term))

delta_wide$eff_type_predobs <- plyr::revalue(delta_wide$effect_type, replace = c(
  "Short-term" = "Observed", "Long-term" = "Observed", "Total" = "Observed"
))

with(delta_wide, table(eff_type_predobs, eff_type_term))

delta_wide$eff_type_term <- factor(delta_wide$eff_type_term, 
                                   levels = c("Short-term", "Long-term", "Total"))
delta_wide$eff_type_predobs <- factor(delta_wide$eff_type_predobs, 
                                   levels = c("Potential", "Observed"))

# save table to plot
saveRDS(delta_wide, "behaviour/files/behaviour_mothers_predictions_total effects.rds")
delta_wide <- readRDS("behaviour/files/behaviour_mothers_predictions_total effects.rds")

# Plot 
ggplot(delta_wide, 
       aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper)) +
  geom_point(size = 2) + 
  geom_linerange(alpha = 0.6) +
  facet_nested_wrap(vars(eff_type_predobs, eff_type_term), nrow = 1) +
  ylim(0, max(delta_wide$d_upper * 1.05)) + 
  theme(panel.grid.minor = element_blank()) + 
  ggtitle("Attacks effects on mothers behaviour") + 
  ylab("Proportion of time in altered behaviour (d)") +
  xlab("Year")
ggsave("behaviour/figures/behaviour_mothers_ts effects total.png",
       height = 10, width = 17, units = "cm")

# WARNING: the total effect is the sum iff the sum is computed before computing
# the delta index. That's because some differences offset others. 


# RAN ALL ABOVE, CONTINUE HERE  -------------------------------------------


# Total effects using observed proportions ----------------------------------

eff_samples2 <- eff_samples

# as y_agg2 has no posterior samples, we need to loop
for(i in 1:N_samples) {
  eff_samples2[, , i, "Short-term"] <- y_agg2 - undisturbed[, , i]
  eff_samples2[, , i, "Total"] <- y_agg2 - undisturbed_95[, , i]
}
eff_samples2[, , , "Long-term"] <- undisturbed - undisturbed_95
eff_samples2[, , , "Sum"] <- eff_samples2[, , , "Short-term"] + eff_samples2[, , , "Long-term"]

# Get delta index
delta_samples2 <- apply(eff_samples2, 2:4, delta_probs)

# Summarize
delta_summ2 <- apply(delta_samples2, c(1, 3), hdmean, name = "d") |> as.data.frame.table()
delta_wide2 <- pivot_wider(delta_summ2, names_from = "Var1", values_from = "Freq")
delta_wide2$year <- as.numeric(as.character(delta_wide2$year))

# plot
ggplot(delta_wide2, aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper)) +
  geom_point(size = 2) + 
  geom_linerange(alpha = 0.6) +
  facet_wrap(vars(effect_type), nrow = 1) +
  ylim(0, max(delta_wide2$d_upper * 1.1)) + 
  theme(panel.grid.minor = element_blank()) + 
  ggtitle("Attacks effects on mothers behaviour") + 
  ylab("Proportion of time in altered behaviour (d)")

par(mfrow = c(2, 2))
plot(delta_wide2$d_mean[delta_wide2$effect_type == "Total"] ~ 
       delta_wide2$d_mean[delta_wide2$effect_type == "Sum"])
abline(0, 1) # OK, it's the sum, it's the same.

plot(delta_wide2$d_lower[delta_wide2$effect_type == "Total"] ~ 
       delta_wide2$d_lower[delta_wide2$effect_type == "Sum"])
abline(0, 1) # OK, it's the sum, it's the same.

plot(delta_wide2$d_upper[delta_wide2$effect_type == "Total"] ~ 
       delta_wide2$d_upper[delta_wide2$effect_type == "Sum"])
abline(0, 1) # OK, it's the sum, it's the same.
par(mfrow = c(1, 1))

# It's also the same.

# Compare short-term effects with data and simulation ---------------------

delta_both <- rbind(
  delta_wide,
  delta_wide2
)

delta_both$calc_type <- rep(c("Simulated", "Observed"), 
                            each = nrow(delta_wide))


ggplot(delta_both, aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper,
                        colour = calc_type, shape = calc_type)) +
  geom_point(size = 2,
             position = position_dodge2(width = 1)) + 
  geom_linerange(alpha = 1,
                 position = position_dodge2(width = 1)) +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0, end = 0.5) +
  facet_wrap(vars(effect_type), nrow = 2) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
  

plot(delta_wide$d_mean[delta_wide$effect_type == "Short-term"] ~ 
       delta_wide2$d_mean[delta_wide2$effect_type == "Short-term"],
     ylab = "Simulated", xlab = "Observed")
abline(0, 1)
# Similar, but with simulated data it's more conservative.



# The sum issue -----------------------------------------------------------

# WARNING: the total effect is the sum iff the sum is computed before computing
# the delta index. That's because some differences offset others. 
# Explore this later.

# Hypothesis: some short and long term effects offset (they have sign)

# Summarize
eff_summ <- apply(eff_samples, c(1, 2, 4), hdmean, name = "eff") |> as.data.frame.table()
eff_wide <- pivot_wider(eff_summ, names_from = "Var1", values_from = "Freq")
eff_wide$year <- as.numeric(as.character(eff_wide$year))

# details for the plot
eff_wide$period <- factor("new", levels = c("old", "new"))
eff_wide$period[eff_wide$year == 1995] <- "old"

effghost1 <- eff_wide[eff_wide$year == 1995, ]
effghost1$year <- 1994.5
effghost2 <- eff_wide[eff_wide$year == 1995, ]
effghost2$year <- 1995.5
effghosts <- rbind(effghost1, effghost2)

eff_wide <- separate_behav(eff_wide)

# plot
selected <- c("Short-term", "Long-term", "Total")

wd <- 0.5
ggplot(eff_wide[eff_wide$effect_type %in% selected, ], 
       aes(x = year, y = eff_mean, ymin = eff_lower, ymax = eff_upper,
           colour = effect_type, shape = effect_type)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 1, size = 0.3) +
  geom_linerange(alpha = 0.6,
                 position = position_dodge2(width = wd)) +
  geom_point(size = 2, alpha = 0.75, 
             position = position_dodge2(width = wd)) + 
  geom_point(data = effghosts, mapping = aes(x = year, y = eff_lower), alpha = 0) +
  facet_nested(rows = vars(behav_raw, position), cols = vars(period), 
               scales = "free", space = "free") +
  theme(panel.grid.minor = element_blank(),
        strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        legend.position = "bottom") + 
  scale_shape_manual(values = c(16, 17, 0, 15)) + 
  scale_color_viridis(discrete = TRUE, begin = 0, end = 0.7, option = "B") +
  ggtitle("Attacks effects on mothers behaviour") + 
  ylab("Proportion difference") +
  scale_x_continuous(breaks = c(1995, 2005, 2008, 2012, 2016)) 
ggsave("behaviour/figures/behaviour_mothers_ts effects short and long.png",
       height = 22, width = 17, units = "cm")

# Indeed, some short and long term effects offset (they have sign).

# Trials: a measure for time of altered behaviour ------------------------

# get max obs effect by year
eff_sums <- aggregate(abs(eff_mean) ~ year + effect_type, data = effs, FUN = sum)
aggregate(abs(eff_mean) ~ year + effect_type, data = effs, FUN = max)


aggregate(eff_mean ~ year + effect_type, data = effs, 
          FUN = function(x) mean(abs(x)) * K)

aggregate(eff_mean ~ year + effect_type, data = effs, 
          FUN = function(x) sum(abs(x)) / 2)

# Test
# discrepancy measure to compare 2 categorical distributions (p vectors).
# being p_i and p_j 2 prob vectors, the discrepancy delta is defined as
# delta = sum(abs(p_i - p_j)) / 2
# It takes 1 if the distributions are completely different (0 overlap),
# and 0 if they are exactly the same 

lpsims <- matrix(rnorm(1e6 * K, sd = 0.00001), 1e6, K)
psims <- apply(lpsims, 1, softmax) |> t()
dsims <- diff(psims) |> abs()
sumsims <- rowSums(dsims) / 2
hist(sumsims, xlim = c(0, 1))


# Trials delta discrepancy for categorical distrib ------------------------

# qué proporción del tiempo pasan haciendo cosas distintas muestras salidas 
# de la misma distribución multinomial?
# Multinom uniforme
delta_sim <- 
      apply(rmultinom(1e6, 1, prob = softmax(rep(0, K))) - 
            rmultinom(1e6, 1, prob = softmax(rep(0, K))),
      2, delta_probs)
hist(delta_sim, xlim = c(0, 1))      
# La prob de que sea una muestra igual a la otra es 1/K si las distrib son unif.

ppp <- softmax(rnorm(K, sd = 20)) # más variabilidad
ppp <- c(1, rep(0, K-1))          # degenerada.
delta_sim <- 
  apply(rmultinom(1e6, 1, prob = ppp) - 
          rmultinom(1e6, 1, prob = ppp),
        2, delta_probs)
hist(delta_sim, xlim = c(0, 1))      
# mientras menos homogénea es la distrib, mayor es la p de que sean iguales, porque 
# la p se concentra en una sola clase.
# En el caso trivial, en una distri degenerada, las muestras son todas iguales y
# delta_prob es cero.

# Se puede pensar como proporción de tiempo la delta_prob?




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
    
    pred_arr[, , y, i] <- apply(X_pred %*% coef, 1, softmax) |> t()
  }
  setTxtProgressBar(pb, i)
}

# summarize and paste rows

tmat_data <- apply(pred_arr, 1:3, hdmean) |> as.data.frame.table()
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


tmat_data$behav_raw_prev <- factor(tmat_data$behav_raw_prev, levels = c(behavs_m_short),
                                   labels = c(behavs_labels))
tmat_data$position_prev <- factor(tmat_data$position_prev, levels = c("S", "UW"))

tmat_data$behav_raw_post <- factor(tmat_data$behav_raw_post, levels = c(behavs_m_short),
                                   labels = c(behavs_labels))
tmat_data$position_post <- factor(tmat_data$position_post, levels = c("S", "UW"))

tmat_data$period <- "2004 to 2018"
tmat_data$period[tmat_data$year_real == "1995"] <- "1995"
tmat_data$period <- factor(tmat_data$period, levels = c("1995", "2004 to 2018"))

saveRDS(tmat_data, "behaviour/files/behaviour_mothers_predictions_tmat data.rds")

# All years
vir <- "D"
ggplot(tmat_data, aes(x = z, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
                      colour = period, fill = year_real, group = year_real,
                      alpha = period, size = period)) +
  geom_line() +
  scale_alpha_manual(values = c(1, 0.6)) +
  scale_size_manual(values = c(1, 0.6)) +
  scale_color_viridis(option = vir, discrete = TRUE, name = "Year", end = 0.7) +
  #scale_color_manual(values = c("red", "gray20"), name = "Year") +
  #geom_ribbon(alpha = 0.2, color = NA) +
  #facet_grid(cols = vars(behav_post), rows = vars(behav_prev)) +
  facet_nested(cols = vars(behav_raw_post, position_post), 
               rows = vars(behav_raw_prev, position_prev)) +
  #ggtitle("Model 30", subtitle = "spline + raneff") +
  #ylim(0, 0.85) +
  theme(panel.grid.minor = element_blank()) +
  guides(size = "none", alpha = "none") +
  ylab("Transition probabilities")
ggsave("behaviour/figures/behaviour_mothers_tmat all years.png", height = 20, width = 25,
       units = "cm")

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

ggsave("behaviour/figures/behaviour_mothers_tmat few years.png", height = 20, width = 25,
       units = "cm")


# Model fit assessment: fitted probabilities and proportions --------------

# I don't use cumulative probability (randomized) residuals because even 
# under an intercept-only model they would have a uniform distribution. 
# I will compare the observed and predicted proportions as a function of z.
# For simplicity, I will only consider the consecutive observed data.

# Extract parameters and define some necessary variables
z_samples <- readRDS("models/behaviour/behaviour_mothers_model_samples_z.rds")
alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
year_cols <- rep(1:N_years, K*K)
N_samples <- nrow(alpha_samples)

# get z values corresponding to observed consecutive behaviours
z_cons <- z_samples[cons_id - 1, ]
y_prev <- y_mat_cons # behaviour in time t-1, like z_cons, it's y_mat[cons_id - 1, ]
y_resp <- y_mat[cons_id, ]
year_num <- bdata$year_num[cons_id]

# I initially tried to make a large array and fill it, but it had 2 problems:
# the RAM usage would go up incredibly when accessing the large array to fill
# it, and that fill step was really time consuming. 
# Using sapply with simplify = "array" appears to solve the problem.
# An interesting feature is that an array takes almost half the size of a list.

pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
probs_arr <- sapply(1:N_samples, simplify = "array", FUN = function(s) {
  # design matrix
  X <- cbind(y_prev, y_prev * z_cons[, s])
  
  # Compute probability matrix by year for sample s
  probs_s <- do.call("rbind", lapply(1:N_years, function(y) {
    # get parameters
    alpha <- matrix(alpha_samples[s, which(year_cols == y)], K, K)
    beta <- matrix(beta_samples[s, which(year_cols == y)], K, K)
    coef <- rbind(alpha, beta)
    year_rows <- sdata$year_start[y] : sdata$year_end[y]
    
    # Compute probabilities
    linpred <- X[year_rows, ] %*% coef
    return(apply(linpred, 1, softmax) |> t())
  }))
  
  setTxtProgressBar(pb, s)
  return(probs_s)
})
gc()
format(object.size(probs_arr), units = "Gb") # 6.9 Gb
# reorder and save
str(probs_arr)
dimnames(probs_arr) <- list(
  row = 1:nrow(y_resp),
  K = behavs_m_clust,
  sample = 1:N_samples
)
str(probs_arr)
saveRDS(probs_arr, "models/behaviour/behaviour_mothers_model_samples_fitted_probabilities.rds")
#when writing to disk it uses almost 8 GB ram extra, and takes a long time.
rm(probs_arr); gc()
probs_arr <- readRDS("models/behaviour/behaviour_mothers_model_samples_fitted_probabilities.rds")
# It consumes the expected RAM when loaded.

# Compare observed vs fitted.
# Classify observations and predictions based on z grouped using only its 
# posterios means. Steps:
# compute the posterior mean by obs for z,
# order z_mean and take that order into the data and probs_arr,
# average all between clases (for probs_arr, it's an aggregation by sample),
# summarize aggregations across posterior samples for probs_arr,
# plot by behaviour category both proportions as a function of z_means.

z_cons_mean <- rowMeans(z_cons)

# order z and all.
z_order <- order(z_cons_mean)
hist(log(z_cons_mean / (1 - z_cons_mean)), breaks = 30) # logit is good

z_ord <- z_cons_mean[z_order]
y_resp_ord <- y_resp[z_order, ]
probs_ord <- probs_arr[z_order, , ]

# group...
N_groups <- 20
group <- cut(z_ord, N_groups)

# aggregate by group
z_agg <- aggregate(z_ord ~ group, FUN = mean)[, -1]
y_resp_agg <- aggregate(y_resp_ord ~ group, FUN = mean)[, -1]

# TOO ram consuming
# tapply_wrapper <- function(x) tapply(x, INDEX = group, FUN = mean)
# check:
# a <- array(rnorm(100 * 4 * 10), dim = c(100, 4, 10))
# apply(a, 2:3, tapply_wrapper)
# probs_agg_samples <- apply(probs_ord, 2:3, tapply_wrapper) # uses a lot of RAM
# probs_agg_samples |> str

probs_agg_samples <- sapply(1:N_samples, function(s) {
  print(s)
  tmp <- aggregate(probs_ord[, , s] ~ group, FUN = mean)[, -1] |> as.matrix()
}, simplify = "array") # it also uses a lot of RAM, but not as much as tapply
probs_agg_samples |> str()
rm(probs_ord); rm(probs_arr); gc() # it was good!
dimnames(probs_agg_samples) <- list(
  z_agg = z_agg,
  behav = behavs_m_clust,
  sample = NULL
)

probs_agg <- apply(probs_agg_samples, 1:2, hdmean, name = "p") |> as.data.frame.table()
probs_agg_wide <- pivot_wider(probs_agg, names_from = "Var1", values_from = "Freq")
probs_agg_wide$prop_type <- "predicted"

colnames(y_resp_agg) <- behavs_m_clust
rownames(y_resp_agg) <- rownames(probs_agg) <- rownames(z_agg) <- NULL
ddd <- cbind(z_agg, y_resp_agg)

ddd$prop_type <- rep("observed", N_groups)
dlong <- pivot_longer(ddd, which(colnames(ddd) %in% behavs_m_clust),
                      values_to = "p_mean", names_to = "behav")

# join with probs
dlong$p_lower <- NA
dlong$p_upper <- NA
dlong <- dlong[, names(probs_agg_wide)]
dlong2 <- rbind(probs_agg_wide, dlong)
dlong2 <- separate_behav(dlong2)
dlong2$prop_type <- factor(dlong2$prop_type, levels = c("predicted", "observed"))
dlong2$z_agg <- dlong2$z_agg |> as.character() |> as.numeric()

ggplot(dlong2, aes(x = z_agg, y = p_mean, colour = prop_type, shape = prop_type,
                  alpha = prop_type, ymin = p_lower, ymax = p_upper)) +
  geom_linerange(size = 2, alpha = 0.5) +
  geom_point(size = 2.5, stroke = 0.5) +
  #scale_color_viridis(discrete = TRUE, end = 0.5, option = "D") +
  scale_color_manual(values = c("black", "red")) +
  scale_shape_manual(values = c(19, 1)) +
  scale_alpha_manual(values = c(0.7, 1)) +
  facet_grid(cols = vars(exposure), rows = vars(behav_raw), 
             scales = "free") +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ylab("Behaviour proportion or probability") +
  xlab("z")

saveRDS(dlong2, "models/behaviour/figures/behaviour_mothers_model fit assessment.rds")
ggsave("behaviour/figures/behaviour_mothers_model fit assessment.png", 
       height = 22, width = 17, units = "cm")

# The prediction is good, but there is a clear overdispersion (poor coverage)


# Model fit assessment: likelihood vs null_likelihood ---------------------

# Compare the log-lik distribution between this model and the model of the mean.
# Here I only use the consecutive observed behaviours, not the UW nor jumps.

# Compute likelihood for fitted model
alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
year_cols <- rep(1:N_years, K*K)
N_samples <- nrow(alpha_samples)
years_unique <- unique(bdata$year)

# matrix to fill with average likelihood by behaviour and iteration
N_cons <- length(sdata$cons_id)
like <- matrix(NA, N_cons, N_samples)

# temporal matrix to fill in every iter
prob_mat <- matrix(NA, N_cons, K)

pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
for(i in 1:N_samples) {
#for(i in 1:100) {
  X <- cbind(
    sdata$y_mat[sdata$cons_id_prev, ],
    sdata$y_mat[sdata$cons_id_prev, ] * z_samples[sdata$cons_id_prev, i]
  )
  
  for(y in 1:N_years) {
    #y = 1
    alpha <- matrix(alpha_samples[i, which(year_cols == y)], K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)], K, K)
    coef <- rbind(alpha, beta)
    
    prob_mat[sdata$year_start[y] : sdata$year_end[y], ] <- 
      apply(X[sdata$year_start[y] : sdata$year_end[y], ] %*% coef, 1, softmax) |> t()
    
  }
  
  for(k in 1:K) {
    like[sdata$y_vec[cons_id] == k, i] <- prob_mat[sdata$y_vec[cons_id] == k, k]
  }

  setTxtProgressBar(pb, i)
}

#View(like)
saveRDS(like, "models/behaviour/behaviour_mothers_model_samples_likelihood.rds")
# memory trials
# lilili <- readRDS("models/behaviour/behaviour_mothers_model_samples_likelihood.rds")
# format(object.size(lilili), units = "Mb")
# dim(lilili)

# Average likelihood (merging obs) by behaviour and sample
like_avg <- matrix(NA, N_samples, K)
for(k in 1:K) {
  #k = 1
  like_avg[, k] <- colMeans(like[sdata$y_vec[cons_id] == k, ])
}
like_avg <- as.data.frame(like_avg)
names(like_avg) <- behavs_m_clust
nice_like <- pivot_longer(like_avg, 1:6, values_to = "likelihood", names_to = "behav")
nice_like$model <- "Nice model"

# Fit null model
counts <- as.data.frame(matrix( table(sdata$y_vec[cons_id]), 1, K))
counts$V
do.call("cbind", counts)
ddd <- data.frame(y = do.call("cbind", counts), 
                  size = length(sdata$y_vec[cons_id]))
str(ddd)
bprior <- brms::prior_string("normal(0, 3)", class = "Intercept")
bm_null <- brms::brm(cbind(y.V1, y.V2, y.V3, y.V4, y.V5, y.V6) | trials(size) ~ 1, 
                     family = brms::multinomial(), 
                     data = ddd,
                     chains = 4, cores = 4,
                     prior = bprior)
brms::stancode(bm_null)
summary(bm_null)

null_samples <- as.matrix(bm_null)[, -K]
colnames(null_samples)
null_params <- apply(cbind(rep(0, 4000), null_samples), 1, softmax) |> t()
null_params <- as.data.frame(null_params)
names(null_params) <- behavs_m_clust
null_like <- pivot_longer(null_params, 1:6, values_to = "likelihood", names_to = "behav")
null_like$model <- "Null model"

likeplot <- rbind(nice_like, null_like)

likeplot <- separate(likeplot, "behav", into = c("behav_raw", "position"), 
                        sep = "_", remove = FALSE)
likeplot$behav_raw <- factor(likeplot$behav_raw, levels = c(behavs_m_short),
                                labels = c(behavs_labels))
likeplot$position <- factor(likeplot$position, levels = c("S", "UW"),
                               labels = c("Surface", "Under water"))

ggplot(likeplot, aes(x = likelihood, color = model, fill = model)) +
  geom_density(alpha = 0.3) +
  facet_grid(rows = vars(position), cols = vars(behav_raw)) +
  xlim(0.05, 0.48) + 
  theme(panel.grid.minor = element_blank()) +
  xlab("Average likelihood across observations") + 
  ylab("Density")

ggsave("behaviour/figures/behaviour_mothers_model likelihood vs null.png", height = 15, width = 23,
       units = "cm")


# Prediction under attacks scenarios (varying attack) ----------------------

# Compute behavior probabilities under a continuous attack scenario starting
# from undisturbed state (steady state). Use inc_b, inc_m, inc_c.
# Do for all years and then average years. 

# subset years 2013-2018
years_unique <- unique(bdata$year)
yidx <- which(years_unique %in% c(2013, 2015:2018))
year_cols <- rep(1:N_years, K*K)

# samples for alpha and beta
alpha_samples <- as.matrix(bm, pars = "alpha")
beta_samples <- as.matrix(bm, pars = "beta")
N_samples <- nrow(alpha_samples)

increase <- abind::abind(
  as.matrix(bm, pars = "increase_b"),
  as.matrix(bm, pars = "increase_m"),
  as.matrix(bm, pars = "increase_c"),
  along = 3
)
str(increase)
dimnames(increase) <- list(
  iter = 1:N_samples,
  year = 1:N_years, 
  scenario = c("Both", "Mother", "Calf")
)

# create prediction array
N_int_pred <- 12 * 1.5 # 1.5 h

pred_attack <- array(NA, dim = c(N_int_pred, K + 1, 3, N_years, N_samples), 
                     dimnames = list(time = 0:(N_int_pred -1),
                                     behav = c("z", behavs_m_clust),
                                     scenario = dimnames(increase)[[3]],
                                     year = years_unique,
                                     sample = 1:N_samples))

str(pred_attack)

# here I will align z with the behavior probabilities in the same interval.
# So, now it's z[t] the one determining the transition probs from t-1 to t.
# The model was written differently to get a differente alignment to
# perform matrix multiplication.
# Here z depends on the previos attacks (which are constant in these  
# scenarios).

pb <- txtProgressBar(min = 0, max = N_samples, style = 3)
for(i in 1:N_samples) {
  #i = 1
  
  # loop over years
  for(y in yidx) {
    #y = 11
    alpha <- matrix(alpha_samples[i, which(year_cols == y)], K, K)
    beta <- matrix(beta_samples[i, which(year_cols == y)], K, K)
    
    # Initialize: z = 0 and behaviour = steady state for undisturbed
    pred_attack[1, "z", , y, i] <- 0
    tmat_0 <- apply(alpha, 1, softmax) |> t()
    mc <- new("markovchain", states = behavs_m_clust, transitionMatrix = tmat_0,
              name = "mc")
    ss <- steadyStates(mc)
    pred_attack[1, 2:(K+1), , y, i] <- ss
    
    # loop over scenarios and time
    for(sc in dimnames(pred_attack)[[3]]) {
      for(t in 2:nrow(pred_attack)) {
        #t = 2
        #sc = "Both"
        pred_attack[t, "z", sc, y, i] <- pred_attack[t-1, "z", sc, y, i] +
          increase[i, y, sc] * (1 - pred_attack[t-1, "z", sc, y, i])
        tmat <- apply(alpha + beta * pred_attack[t, "z", sc, y, i], 1, softmax) |> t()
        
        pred_attack[t, 2:(K+1), sc, y, i] <- pred_attack[t-1, 2:(K+1), sc, y, i] %*% tmat
      }
    } # end loop over scenarios
    
  }# end loop over years
  setTxtProgressBar(pb, i)
}# end loop over posterior samples

# remove empty years and save:
pred_attack <- pred_attack[, , , yidx, ]
str(pred_attack)
anyNA(pred_attack)
saveRDS(pred_attack, "models/behaviour/behaviour_mothers_predictions_attack scenario samples_bmc.rds")

# average years
p1 <- apply(pred_attack, c(1, 2, 3, 5), mean)
str(p1)
# summarize posterior samples
pred_agg_tmp <- apply(p1, 1:3, hdmean, name = "p") |> as.data.frame.table()
pred_agg <- pivot_wider(pred_agg_tmp, names_from = "Var1", values_from = "Freq")
pred_agg <- separate_behav(pred_agg)
pred_agg$min <- as.numeric(as.character(pred_agg$time)) * 5

ggplot(pred_agg[complete.cases(pred_agg), ],
       aes(x = min, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(alpha = 0.1, color = NA) +
  facet_nested_wrap(vars(behav_raw, exposure), nrow = 2) + 
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.7) +
  theme(panel.grid.minor = element_blank())

# Export table
pred_agg$mc <- "Mothers"
saveRDS(pred_agg, "behaviour/files/behaviour_mothers_predictions_attack scenario table_bmc.rds")


