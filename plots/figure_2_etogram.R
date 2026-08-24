# Figure 2: etogram of southern right whales by year -----------------------
#
# Two panels (rows): mothers on top, calves below. Each bar is a year and shows
# the observed proportion of time spent in each behaviour category, marginal to
# (i.e. not conditional on) gull attacks. Proportions are computed on the data
# with the behaviours that were not observed imputed by the behaviour models
# (behavm_full / behavc_full), so every interval contributes.
#
# Behaviour class (rest, slow travel, fast travel, in-place activity) is keyed by
# colour, and exposure (surface / under water) by transparency.
#
# Years run consecutively, so the years without sampling show up as gaps. The
# 1995 gap is too long to draw, so 1995 sits in its own panel, as in the other
# time-series figures (plots/plots_script.R, "steady states by year").
#
# Run from the repo root (open srw_behaviour.Rproj).

library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(viridis)
library(ggh4x)           # facet_nested


# Behaviour coding --------------------------------------------------------
# (same coding used by the analysis scripts)

behavs_clust  <- c("R_S", "R_UW",
                   "ST_S", "ST_UW",
                   "FT_S", "FT_UW",
                   "IPA_S", "IPA_UW")
behavs_short  <- c("R", "ST", "FT", "IPA")
behavs_labels <- c("Rest", "Slow travel", "Fast travel", "In-place activity")

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


# Data --------------------------------------------------------------------

bm <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")
bc <- readRDS("models/behaviour/behaviour_calves_data with imputed behaviours.rds")

# observed proportions by year, using the imputed behaviour column
props_by_year <- function(df, behav_col, mc) {
  b <- factor(as.character(df[[behav_col]]), levels = behavs_clust)
  tab <- table(df$year, b)
  p <- tab / rowSums(tab)
  out <- as.data.frame(p, stringsAsFactors = FALSE)
  names(out) <- c("year", "behav", "prop")
  out$year <- as.numeric(out$year)
  out$mc <- mc
  return(out)
}

d <- rbind(props_by_year(bm, "behavm_full", "Mothers"),
           props_by_year(bc, "behavc_full", "Calves"))

d$behav <- factor(d$behav, levels = behavs_clust)
d <- separate_behav(d)
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))

# 1995 in its own panel: the 9-year gap that follows it would eat the figure
d$period <- factor("new", levels = c("old", "new"))
d$period[d$year == 1995] <- "old"

# ghost rows widening the 1995 panel, so its bar is not squeezed against the
# panel borders (same trick as in plots_script.R)
ghosts <- d[d$period == "old", ][1, ]
ghosts <- rbind(within(ghosts, year <- 1994.5),
                within(ghosts, year <- 1995.5))
ghosts$prop <- 0


# Plot --------------------------------------------------------------------

colors <- viridis(4, end = 0.9)
names(colors) <- behavs_labels

p <- ggplot(d, aes(x = year, y = prop, fill = behav_raw, alpha = exposure)) +
  geom_col(width = 0.85, colour = NA, position = position_stack(reverse = TRUE)) +
  geom_blank(data = ghosts) +
  facet_nested(rows = vars(mc), cols = vars(period),
               scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors, name = "Behaviour") +
  scale_alpha_manual(values = c("Surface" = 1, "Under water" = 0.45),
                     name = "Exposure") +
  scale_x_continuous(breaks = sort(unique(d$year))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(override.aes = list(alpha = 1), order = 1),
         alpha = guide_legend(override.aes = list(fill = "gray30"), order = 2)) +
  ylab("Proportion of time") +
  xlab("Year") +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.y = element_text(angle = 270),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = 0, b = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/figure_2_etogram.png", plot = p,
       height = 14, width = 17, units = "cm", dpi = 400)
