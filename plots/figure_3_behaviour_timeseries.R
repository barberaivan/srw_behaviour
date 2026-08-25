# Figure 3: behaviour by year, observed and estimated -----------------------
#
# Time series, year by year, of the probability of each of the eight behaviour
# categories under three conditions:
#
#   Undisturbed  steady state of the transition matrix at z = 0 (a whale that
#                has been free of attacks long enough to have recovered)
#   Disturbed    steady state at z = 1 (a whale under sustained attack)
#   Observed     the proportion of intervals in which the behaviour was recorded,
#                using the model-imputed behaviour column so every interval counts
#
# Both reference states are counterfactual, estimated quantities: no whale is
# ever observed at exactly z = 0 or z = 1 (see the Methods).
#
# The marker in the manuscript asks for two layouts, so two versions are saved:
#
#   version A (..._faceted_vNN.png)  mothers on the left, calves on the right,
#                                    merged with patchwork; y scales are shared
#                                    row by row between the two halves
#   version B (..._dodged_vNN.png)   mothers and calves in the same panels,
#                                    keyed by point shape and line type, with
#                                    the three conditions dodged and joined by
#                                    lines
#
# 1995 sits in its own panel column: the nine-year gap that follows it would
# otherwise eat the figure (same trick as in plots/plots_script.R).
#
# Inputs (written by behaviour/behaviour_{mothers,calves}_analysis.R):
#   behaviour/files/behaviour_mothers_predictions disturbed and undisturbed by year.rds
#   behaviour/files/behaviour_calves_predictions disturbed and undisturbed by year.rds
#
# Run from the repo root (open srw_behaviour.Rproj).

source("plots/theme_paper.R")
library(ggh4x)           # facet_nested
library(patchwork)


# Data --------------------------------------------------------------------

lm_ <- readRDS("behaviour/files/behaviour_mothers_predictions disturbed and undisturbed by year.rds")
lc_ <- readRDS("behaviour/files/behaviour_calves_predictions disturbed and undisturbed by year.rds")

dm <- lm_$plot_data
dc <- lc_$plot_data
dm$mc <- "Mothers"
dc$mc <- "Calves"

d <- rbind(dm, dc)
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))
d$scenario <- factor(d$scenario, levels = c("Disturbed", "Undisturbed", "Observed"))

# ghost rows widening the 1995 panel so its points are not squeezed against the
# panel borders (they are already in the mothers' object)
ghosts <- lm_$ghosts
ghosts$scenario <- factor(as.character(ghosts$scenario),
                          levels = levels(d$scenario))
ghosts$mc <- factor("Mothers", levels = levels(d$mc))

# shared y range per row (behaviour class x exposure), so that the mothers and
# the calves halves of version A can be read against each other
ylims <- do.call(rbind, lapply(
  split(d, list(d$behav_raw, d$exposure), drop = TRUE),
  function(x) data.frame(
    behav_raw = x$behav_raw[1], exposure = x$exposure[1],
    prob_mean = c(min(c(x$prob_lower, x$prob_mean), na.rm = TRUE),
                  max(c(x$prob_upper, x$prob_mean), na.rm = TRUE))
  )
))
ylims$prob_mean[ylims$prob_mean < 0] <- 0

dodge <- position_dodge2(width = 0.3)

# the invisible rows that impose those limits, one set per period panel so that
# they never widen the x range of a panel
ylim_rows <- function(dat) {
  xs <- tapply(dat$year_num, droplevels(dat$period), min)
  do.call(rbind, lapply(names(xs), function(pp) {
    out <- ylims
    out$period <- factor(pp, levels = levels(dat$period))
    out$year_num <- xs[[pp]]
    out
  }))
}


# Version A: mothers and calves as separate blocks ------------------------

panel_ts <- function(dat, ghost_rows = NULL, show_y = TRUE, show_strip_y = TRUE) {
  p <- ggplot(dat,
              aes(x = year_num, y = prob_mean, ymin = prob_lower,
                  ymax = prob_upper, shape = scenario, fill = scenario)) +
    geom_point(position = dodge, size = 2.4) +
    geom_linerange(position = dodge, show.legend = FALSE) +
    # invisible points that impose the shared row-wise y range
    geom_blank(data = ylim_rows(dat), inherit.aes = FALSE,
               mapping = aes(x = year_num, y = prob_mean)) +
    facet_nested(rows = vars(behav_raw, exposure), cols = vars(period),
                 scales = "free", space = "free_x") +
    scale_fill_manual(values = colors_condition) +
    scale_shape_manual(values = c("Disturbed" = 23, "Undisturbed" = 22,
                                  "Observed" = 21)) +
    theme(strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 270),
          legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
          legend.box.margin = margin(t = 0, unit = "mm")) +
    ylab(NULL) + xlab(NULL)

  if (!is.null(ghost_rows) && nrow(ghost_rows) > 0) {
    p <- p + geom_point(data = ghost_rows,
                        mapping = aes(x = year_num, y = prob_lower), alpha = 0)
  }
  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  if (!show_strip_y) {
    p <- p + theme(strip.text.y = element_blank(),
                   strip.background.y = element_blank())
  }
  return(p)
}

pm <- panel_ts(d[d$mc == "Mothers", ], ghost_rows = ghosts,
               show_y = TRUE, show_strip_y = FALSE) +
  scale_x_continuous(breaks = c(1995, 2005, 2010, 2015, 2020, 2025)) +
  ggtitle("A. Mothers")

pc <- panel_ts(d[d$mc == "Calves", ], show_y = FALSE, show_strip_y = TRUE) +
  scale_x_continuous(breaks = c(2015, 2020, 2025)) +
  ggtitle("B. Calves") +
  guides(fill = "none", shape = "none")

fig3a <- (pm + pc) +
  plot_layout(widths = c(19, 9), guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, hjust = 0, vjust = 1))

fig3a <- wrap_elements(fig3a) +
  labs(tag = "Behaviour probability") +
  theme(plot.tag = element_text(size = 12, angle = 90),
        plot.tag.position = "left")

fig3a <- fig3a / grid::textGrob("Year", gp = grid::gpar(fontsize = 12)) +
  plot_layout(heights = c(50, 1))

ggsave("plots/figure_3_behaviour_timeseries_faceted_v01.png", plot = fig3a,
       height = 23, width = 18, units = "cm", dpi = 400)


# Version B: mothers and calves in the same panels -------------------------

fig3b <- ggplot(d,
                aes(x = year_num, y = prob_mean, ymin = prob_lower,
                    ymax = prob_upper, shape = mc, fill = scenario,
                    colour = scenario, linetype = mc, group = interaction(scenario, mc))) +
  geom_line(position = dodge, linewidth = 0.35, alpha = 0.8) +
  geom_linerange(position = dodge, show.legend = FALSE) +
  geom_point(position = dodge, size = 2.2, colour = "gray20", stroke = 0.3) +
  geom_point(data = ghosts, inherit.aes = FALSE,
             mapping = aes(x = year_num, y = prob_lower), alpha = 0) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(period),
               scales = "free", space = "free_x") +
  scale_fill_manual(values = colors_condition, name = NULL) +
  scale_colour_manual(values = colors_condition, name = NULL) +
  scale_shape_manual(values = c("Mothers" = 21, "Calves" = 24), name = NULL) +
  scale_linetype_manual(values = c("Mothers" = "solid", "Calves" = "22"),
                        name = NULL) +
  scale_x_continuous(breaks = c(1995, 2005, 2010, 2015, 2020, 2025)) +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, linetype = 0)),
         colour = "none",
         shape = guide_legend(order = 2, override.aes = list(fill = "gray50")),
         linetype = guide_legend(order = 2)) +
  ylab("Behaviour probability") +
  xlab("Year") +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270),
        legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
        legend.box.margin = margin(t = 0, unit = "mm"))

ggsave("plots/figure_3_behaviour_timeseries_dodged_v01.png", plot = fig3b,
       height = 23, width = 18, units = "cm", dpi = 400)
