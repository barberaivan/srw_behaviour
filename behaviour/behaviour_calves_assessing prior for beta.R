# Assessing the prior for beta parameters in calves model.

# In the calves model, although the increase and decay parameters were 
# constrained, the beta parameters were sensitive to the prior. A wide prior as
# beta ~ Normal(0, 2) would allow for large beta values and very low increase and
# decrease parameters. In this scenarios the distributions for increase and decrease
# were bimodal, suggesting that a slow z dynamic was mixing with the fast one. 
# We fitted 3 models, varying the prior sd for beta in {1, 1.5, 2}. We chose the 
# widest prior that would make a unimodal posterior for the parameters controlling
# the z dynamics, which was beta ~ Normal(0, 1.5). Here we compare the betas for
# calves using that model with the betas for mothers in the same years, to check 
# whether betas for mother and calves lay in similar ranges. (They do).


# Packages ----------------------------------------------------------------

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
library(magrittr)
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
bdata <- read.csv("data/srw_behaviour_data.csv",
                  stringsAsFactors = FALSE, sep = ",")


# Load models -------------------------------------------------------------

mm <- readRDS("models/behaviour/behaviour_mothers_model_samples.rds")
mc <- readRDS("models/behaviour/behaviour_calves_model_samples15.rds") 

# calves betas
beta_c <- as.matrix(mc, "beta")
beta_c <- beta_c[, colSums(beta_c) != 0]

# mothers betas
beta_m <- as.matrix(mm, pars = "beta")
K <- 8
year_cols <- rep(1:length(unique(bdata$year)), K*K)
beta_m <- beta_m[, year_cols %in% 11:15]
colnames(beta_m) <- year_cols[year_cols %in% 11:15]
beta_m <- beta_m[, colSums(beta_m) != 0]


ncol(beta_m) / 5
ncol(beta_c) # OK


# get means by year in mothers.
beta_m2 <- sapply(as.character(11:15), function(y) {
  beta_m[, colnames(beta_m) == y]
}, simplify = "array")
str(beta_m2)

beta_mm <- apply(beta_m2, 1:2, mean) 
str(beta_mm)
str(beta_c)

# Compare -----------------------------------------------------------------

# Get the sd across betas for mothers:
sd_samples <- apply(beta_mm, 1, sd)
plot(density(sd_samples, from = 0), xlim = c(0.5, 2))
abline(v = 1.5, lty = 2, col = 4) # 1.5 was one of the priors used for calves.

sd_samples_c <- apply(beta_c, 1, sd)
lines(density(sd_samples_c, from = 0), col = 2)

get_map(sd_samples)
posterior::summarise_draws(data.frame(sd_samples = sd_samples))
# variable    mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>      <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# sd_samples  1.25   1.25 0.118 0.118  1.06  1.45  1.00    2397.    4440.

posterior::summarise_draws(data.frame(sd_samples = sd_samples_c))
# variable    mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>      <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# sd_samples  1.32   1.32 0.119 0.119  1.13  1.52  1.00    3811.    5174.

# Very similar distributions.

# The widths of beta posterior distributions, measured with their sds, 
# in the calves model are just a bit wider than in the mothers' model 
# for the same years, so the prior on beta for calves is probably not 
# overshrinking.

# Plot posteriors for beta

dlist_c <- lapply(1:ncol(beta_c), function(i) density(beta_c[, i]))
dlist_m <- lapply(1:ncol(beta_mm), function(i) density(beta_mm[, i]))

par(mfrow = c(2, 1))

plot(dlist_c[[1]], xlim = c(-4, 4), ylim = c(0, 1.5), main = "Calves")
for(i in 2:ncol(beta_c)) lines(dlist_c[[i]], col = i)
curve(dnorm(x, 0, 1.5), lty = 2, lwd = 2, add = TRUE)
curve(dnorm(x, 0, mean(sd_samples)), lty = 1, lwd = 4, add = TRUE, col = 2)

plot(dlist_m[[1]], xlim = c(-4, 4), ylim = c(0, 1.5), main = "Mothers")
for(i in 2:ncol(beta_mm)) lines(dlist_m[[i]], col = i)
curve(dnorm(x, 0, 2), lty = 2, lwd = 2, add = TRUE)

par(mfrow = c(1, 1))


