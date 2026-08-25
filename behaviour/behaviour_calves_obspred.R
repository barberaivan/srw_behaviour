# Calves: behaviour predicted at the estimated z, and the total effects ------
#
# The calves' equivalent of the "Predicted behaviours for estimated z values
# (behav_obspred)" and "Total effects using obspred proportions" sections of
# behaviour/behaviour_calves_analysis.R, which were written but never run (their
# saveRDS() calls are commented out there, and the objects are missing from the
# store). Figure 5 needs them, so they are computed here.
#
# What it does: starting each follow from its first (observed or imputed)
# behaviour, it simulates the behaviour of every later interval from the fitted
# transition matrix of that year evaluated at the z estimated for the previous
# interval. The result is the behaviour distribution the model predicts *under
# the attacks the calves actually received*, free of the model's own prediction
# error — which is what the short-term effect is compared against (see Methods,
# "Quantifying the effect of attacks").
#
# The recursion is the same as in the analysis script; the only difference is
# that the posterior samples are advanced together rather than one at a time,
# which turns a loop of 11,829 x 8,000 iterations into one of 11,829. The draws
# are therefore not identical to a run of the original triple loop, but they come
# from the same distribution.
#
# Writes:
#   models/behaviour/behaviour_calves_model_samples_predict_observed_attack.rds
#   models/behaviour/behaviour_calves_model_samples_predict_observed_attack_years.rds
#   behaviour/files/behaviour_calves_model_pred dist wide.rds
#   behaviour/files/behaviour_calves_predictions_total effects.rds
#
# Run from the repo root (open srw_behaviour.Rproj). Needs ~4 GB of RAM and a
# few minutes.

library(rstan)
library(bayestestR)
library(tidyr)
library(abind)

set.seed(342)


# Functions ---------------------------------------------------------------

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

# delta index: 1 - overlap between two behaviour distributions
delta_probs <- function(x) sum(abs(x)) / 2

behavs_c_clust <- c("R_S", "R_UW",
                    "ST_S", "ST_UW",
                    "FT_S", "FT_UW",
                    "IPA_S", "IPA_UW")
K <- length(behavs_c_clust)


# Data and samples --------------------------------------------------------

bdata <- readRDS("models/behaviour/behaviour_calves_data with imputed behaviours.rds")
z_samples <- readRDS("models/behaviour/behaviour_calves_model_samples_z.rds")

stopifnot(nrow(bdata) == nrow(z_samples))

bc <- readRDS("models/behaviour/behaviour_calves_model_samples.rds")
alpha_samples <- as.matrix(bc, pars = "alpha")
beta_samples <- as.matrix(bc, pars = "beta")
rm(bc); gc()

N_samples <- nrow(alpha_samples)
stopifnot(N_samples == ncol(z_samples))

years_unique <- unique(bdata$year)
N_years <- length(years_unique)

# alpha and beta are stored as alpha[year, from, to] with the year index running
# fastest, so the columns of one year, taken in order, fill a K x K matrix by
# column: entry (from, to) sits at column (to - 1) * K + from
year_cols <- rep(1:N_years, K * K)

# follow limits: bdata is ordered by year, follow and t
runs <- rle(bdata$follow_num)
int_end <- cumsum(runs$lengths)
int_start <- int_end - runs$lengths + 1
N_fol <- length(int_start)
stopifnot(N_fol == length(unique(bdata$follow)))


# Simulate behaviour at the estimated z ------------------------------------

# one S x K^2 matrix of alpha and of beta per year
A <- lapply(1:N_years, function(y) alpha_samples[, which(year_cols == y)])
B <- lapply(1:N_years, function(y) beta_samples[, which(year_cols == y)])
rm(alpha_samples, beta_samples); gc()

pred <- matrix(NA_integer_, nrow(bdata), N_samples)
pred[int_start, ] <- rep(as.integer(bdata$behavc_full[int_start]), N_samples)

s_idx <- seq_len(N_samples)
years_i <- bdata$year_num

pb <- txtProgressBar(min = 0, max = N_fol, style = 3)
for (f in 1:N_fol) {
  for (i in (int_start[f] + 1):int_end[f]) {
    y <- years_i[i]
    prev <- pred[i - 1, ]
    zprev <- z_samples[i - 1, ]

    # linear predictor of every destination behaviour, for every sample
    eta <- matrix(NA_real_, N_samples, K)
    for (k in 1:K) {
      cols <- prev + (k - 1) * K
      eta[, k] <- A[[y]][cbind(s_idx, cols)] + B[[y]][cbind(s_idx, cols)] * zprev
    }

    # softmax by row
    p <- exp(eta - apply(eta, 1, max))
    p <- p / rowSums(p)

    # one categorical draw per sample, by inversion
    cp <- p
    for (k in 2:K) cp[, k] <- cp[, k - 1] + p[, k]
    u <- runif(N_samples)
    pred[i, ] <- as.integer(rowSums(cp < u) + 1L)
  }
  setTxtProgressBar(pb, f)
}
close(pb)

saveRDS(pred,
        "models/behaviour/behaviour_calves_model_samples_predict_observed_attack.rds")


# Summarise by year --------------------------------------------------------

pred_dist_samples_year <- array(NA_real_, dim = c(N_years, K, N_samples),
                                dimnames = list(year = years_unique,
                                                behav = behavs_c_clust,
                                                sample = 1:N_samples))
for (y in 1:N_years) {
  rows <- which(bdata$year_num == y)
  n <- length(rows)
  sub <- pred[rows, ]
  for (k in 1:K) pred_dist_samples_year[y, k, ] <- colSums(sub == k) / n
}

saveRDS(pred_dist_samples_year,
        "models/behaviour/behaviour_calves_model_samples_predict_observed_attack_years.rds")

pred_dist <- apply(pred_dist_samples_year, 1:2, hdmean, name = "pred") |>
  as.data.frame.table()
pred_dist_wide <- pivot_wider(pred_dist, names_from = "Var1", values_from = "Freq")
pred_dist_wide$year <- as.numeric(as.character(pred_dist_wide$year))
saveRDS(pred_dist_wide, "behaviour/files/behaviour_calves_model_pred dist wide.rds")

rm(pred); gc()


# Total effects ------------------------------------------------------------
# (identical to the mothers' section of the same name; the baseline year for the
# long-term and total effects is the first calves year, 2013, not 1995, so those
# two are not comparable with the mothers' and are not used in the paper)

pred_samples <- readRDS("models/behaviour/behaviour_calves_predictions_attack scenario samples.rds")
undisturbed <- pred_samples$attack[1, 2:(K + 1), , ]
disturbed <- pred_samples$no_attack[1, 2:(K + 1), , ]
undisturbed_base <- undisturbed
for (y in 1:N_years) undisturbed_base[, y, ] <- undisturbed[, 1, ]
rm(pred_samples); gc()

obspred_samples <- aperm(pred_dist_samples_year, perm = c(2, 1, 3))

eff_samples <- abind("Potential" = undisturbed,
                     "Short-term" = undisturbed,
                     "Long-term" = undisturbed,
                     "Total" = undisturbed,
                     along = length(dim(undisturbed)) + 1)

dimnames(eff_samples) <- list(
  behav = behavs_c_clust,
  year = years_unique,
  sample = 1:N_samples,
  effect_type = c("Potential", "Short-term", "Long-term", "Total")
)

eff_samples[, , , "Potential"]  <- disturbed - undisturbed
eff_samples[, , , "Short-term"] <- obspred_samples - undisturbed
eff_samples[, , , "Long-term"]  <- undisturbed - undisturbed_base
eff_samples[, , , "Total"]      <- obspred_samples - undisturbed_base

delta_samples <- apply(eff_samples, 2:4, delta_probs)
delta_summ <- apply(delta_samples, c(1, 3), hdmean, name = "d") |>
  as.data.frame.table()
delta_wide <- pivot_wider(delta_summ, names_from = "Var1", values_from = "Freq")
delta_wide$year <- as.numeric(as.character(delta_wide$year))

# same labelling as the mothers' table: the "potential" effect is a short-term
# effect too, and everything that is not potential is realised ("observed")
delta_wide$eff_type_term <- plyr::revalue(delta_wide$effect_type, replace = c(
  "Potential" = "Short-term"
))
delta_wide$eff_type_predobs <- plyr::revalue(delta_wide$effect_type, replace = c(
  "Short-term" = "Observed", "Long-term" = "Observed", "Total" = "Observed"
))
delta_wide$eff_type_term <- factor(delta_wide$eff_type_term,
                                   levels = c("Short-term", "Long-term", "Total"))
delta_wide$eff_type_predobs <- factor(delta_wide$eff_type_predobs,
                                      levels = c("Potential", "Observed"))

saveRDS(delta_wide, "behaviour/files/behaviour_calves_predictions_total effects.rds")
print(as.data.frame(delta_wide))
