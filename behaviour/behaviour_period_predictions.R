# Attack-scenario predictions averaged over periods -------------------------
#
# Post-processing of the attack-scenario prediction samples produced by
# behaviour/behaviour_{mothers,calves}_analysis.R ("Prediction under attacks
# scenarios"). No model is refitted here.
#
# Those objects hold, for every posterior sample and every year, the behaviour
# probabilities of a whale that starts undisturbed and is then attacked without
# pause ("Persistent attacks"), and of a whale that starts fully disturbed and is
# then left alone ("Attacks cessation"), over 31 five-minute intervals.
#
# The manuscript groups years into periods, so here the years of a period are
# averaged. The average is taken **within each posterior sample**, and the
# posterior is summarised (mean and 95 % HDI) only at the very end: averaging the
# posterior means of each year instead would give the same mean but the wrong
# credible interval.
#
# Writes:
#   behaviour/files/behaviour_mothers_predictions_attack scenarios table_periods.rds
#   behaviour/files/behaviour_calves_predictions_attack scenarios table_periods.rds
#   behaviour/files/behaviour_mothers_predictions_attack scenarios table_avg.rds
#   behaviour/files/behaviour_calves_predictions_attack scenarios table_avg.rds
#
# The two *_avg tables are the same computation over the years shared by mothers
# and calves (2013, 2015-2018); they are what plots/plots_script.R reads for its
# Figure 1.
#
# Run from the repo root (open srw_behaviour.Rproj). Needs ~3 GB of RAM.

library(bayestestR)      # highest density intervals
library(tidyr)


# Functions ---------------------------------------------------------------
# (the same summaries used by the analysis scripts)

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

behavs_short  <- c("R", "ST", "FT", "IPA")
behavs_labels <- c("Rest", "Slow travel", "Fast travel", "In place activity")

separate_behav <- function(df, col_name = "behav") {
  df <- separate(df, col_name, into = c("behav_raw", "exposure"),
                 sep = "_", remove = FALSE)
  df$behav_raw <- factor(df$behav_raw, levels = behavs_short,
                         labels = behavs_labels)
  df$exposure <- factor(df$exposure, levels = c("S", "UW"),
                        labels = c("Surface", "Under water"))
  return(df)
}


# Periods -----------------------------------------------------------------
# As asked in the manuscript. A period keeps only the years actually sampled,
# and periods with no sampled year for a given target are dropped (mothers get
# five periods, calves three).

periods <- list(
  "1995"      = 1995,
  "2004-2010" = 2004:2010,
  "2011-2013" = 2011:2013,
  "2014-2019" = 2014:2019,
  "2020-2025" = 2020:2025
)


# Averaging ---------------------------------------------------------------

# Average the years of `years_keep` within each posterior sample, then summarise.
# arr: [time, behav (z first), year, sample]
average_years <- function(arr, years_keep) {
  ys <- as.character(years_keep)
  # mean over years, keeping time, behaviour and sample. Adding the year slices
  # one by one is much faster than apply() over an array of this size.
  avg <- arr[, , ys[1], ]
  if (length(ys) > 1) {
    for (k in ys[-1]) avg <- avg + arr[, , k, ]
    avg <- avg / length(ys)
  }
  out <- apply(avg, 1:2, hdmean, name = "prob") |> as.data.frame.table()
  return(out)
}

# Build the table of one target for an arbitrary grouping of years.
# groups: named list of year vectors; the name becomes the `period` column.
period_table <- function(pred_samples, groups, mc) {
  years_have <- as.numeric(dimnames(pred_samples$attack)[[3]])

  out <- do.call(rbind, lapply(names(groups), function(g) {
    years_keep <- intersect(groups[[g]], years_have)
    if (length(years_keep) == 0) return(NULL)

    at <- average_years(pred_samples$attack, years_keep)
    no <- average_years(pred_samples$no_attack, years_keep)
    at$scenario <- "Attack"
    no$scenario <- "After attack"

    tab <- rbind(at, no)
    tab$period <- g
    tab$n_years <- length(years_keep)
    tab$years <- paste(years_keep, collapse = ", ")
    return(tab)
  }))

  out <- pivot_wider(out, names_from = "Var1", values_from = "Freq")
  out$scenario <- factor(out$scenario, levels = c("Attack", "After attack"))
  out$period <- factor(out$period, levels = names(groups))
  out$time <- as.numeric(as.character(out$time))
  out$min <- out$time * 5
  out$mc <- mc
  out <- separate_behav(out)
  return(as.data.frame(out))
}


# Mothers -----------------------------------------------------------------

pred_samples <- readRDS("models/behaviour/behaviour_mothers_predictions_attack scenario samples.rds")

pm <- period_table(pred_samples, periods, "Mothers")
saveRDS(pm, "behaviour/files/behaviour_mothers_predictions_attack scenarios table_periods.rds")

# the years shared with the calves, for plots/plots_script.R (Figure 1)
pm_avg <- period_table(pred_samples, list("2013-2018" = c(2013, 2015:2018)),
                       "Mothers")
saveRDS(pm_avg, "behaviour/files/behaviour_mothers_predictions_attack scenarios table_avg.rds")

rm(pred_samples); gc()


# Calves ------------------------------------------------------------------

pred_samples <- readRDS("models/behaviour/behaviour_calves_predictions_attack scenario samples.rds")

pc <- period_table(pred_samples, periods, "Calves")
saveRDS(pc, "behaviour/files/behaviour_calves_predictions_attack scenarios table_periods.rds")

pc_avg <- period_table(pred_samples, list("2013-2018" = c(2013, 2015:2018)),
                       "Calves")
saveRDS(pc_avg, "behaviour/files/behaviour_calves_predictions_attack scenarios table_avg.rds")

rm(pred_samples); gc()

# which years ended up in each period
print(unique(pm[, c("period", "n_years", "years")]))
print(unique(pc[, c("period", "n_years", "years")]))
