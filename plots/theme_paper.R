# Shared theme and behaviour coding for the manuscript figures ---------------
#
# Sourced by the figure scripts in plots/. It holds the things every figure
# needs, so they are defined once: the behaviour category coding, a helper to
# split a behaviour label into class and exposure, the palettes, and the ggplot
# theme (theme_mine(), the same one defined in plots/plots_script.R).
#
# Usage, from the repo root:  source("plots/theme_paper.R")

library(ggplot2)
library(tidyr)
library(viridis)


# Behaviour coding --------------------------------------------------------
# (the same coding used by the analysis scripts)

behavs_clust  <- c("R_S", "R_UW",
                   "ST_S", "ST_UW",
                   "FT_S", "FT_UW",
                   "IPA_S", "IPA_UW")
behavs_short  <- c("R", "ST", "FT", "IPA")
behavs_labels <- c("Rest", "Slow travel", "Fast travel", "In place activity")

# separate a behaviour column into class and exposure
separate_behav <- function(df, col_name = "behav") {
  df <- separate(df, col_name, into = c("behav_raw", "exposure"),
                 sep = "_", remove = FALSE)
  df$behav_raw <- factor(df$behav_raw, levels = behavs_short,
                         labels = behavs_labels)
  df$exposure <- factor(df$exposure, levels = c("S", "UW"),
                        labels = c("Surface", "Under water"))
  return(df)
}


# Palettes ----------------------------------------------------------------

# the two attack scenarios (persistent attacks / attacks cessation), and the
# disturbed-undisturbed-observed triplet, as used everywhere in the paper
colors_scenario <- viridis(2, end = 0.5)
colors_condition <- c("Undisturbed" = viridis(2, end = 0.5)[1],
                      "Disturbed"   = viridis(2, end = 0.5)[2],
                      "Observed"    = viridis(3, end = 0.93)[3])

# mothers / calves, when they are keyed by colour instead of by facet
colors_mc <- c("Mothers" = "#2D3184", "Calves" = "#E07B39")

# behaviour classes
colors_behav <- setNames(viridis(4, end = 0.9), behavs_labels)


# Theme -------------------------------------------------------------------
# from https://rpubs.com/mclaire19/ggplot2-custom-themes
# (kept identical to the theme_mine() in plots/plots_script.R)

theme_mine <- function() {
  font <- "Arial"
  marg <- 2 # figure margin in mm

  theme_bw() %+replace%

    theme(
      panel.grid.minor = element_blank(),

      plot.title = element_text(family = font, size = 16, hjust = -0.1,
                                vjust = 1),

      axis.title = element_text(family = font, size = 12),

      # separates the y axis title from the numbers
      axis.title.y = element_text(
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),

      axis.text = element_text(family = font, size = 9),

      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),

      strip.text = element_text(size = 10, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2, 0, 1.2, 0, "mm")),
      strip.text.y = element_text(margin = margin(0, 1.2, 0, 1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),

      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

# the same, with light strips and black strip text
theme_mine2 <- function() {
  font <- "Arial"
  theme_mine() %+replace%

    theme(
      strip.text = element_text(size = 10, family = font, color = "black"),
      strip.text.x = element_text(margin = margin(1.2, 0, 1.2, 0, "mm")),
      strip.text.y = element_text(margin = margin(0, 1.2, 0, 1.2, "mm")),
      strip.background = element_rect(fill = "white", color = "white")
    )
}

theme_set(theme_mine())
