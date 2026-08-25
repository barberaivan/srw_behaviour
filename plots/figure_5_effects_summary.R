# Figure 5: summary of the effects of gull attacks on behaviour -------------
#
# The four effects defined in the Methods, year by year, each measured with the
# discrepancy index delta — half the sum of the absolute differences between the
# probabilities of the eight behaviour categories, read as the proportion of time
# the whale spends doing something other than what it would do in the reference
# state:
#
#   Potential,  short-term  disturbed vs undisturbed, same year
#   Observed,   short-term  behaviour predicted under the attacks actually
#                           received vs undisturbed, same year
#   Observed,   long-term   undisturbed of the year vs undisturbed of 1995
#   Observed,   total       predicted under observed attacks vs undisturbed 1995
#
# Mothers and calves are keyed by colour. The long-term and total effects are
# mothers only: they are measured against 1995 and the calves' series starts in
# 2013, with no pre-increase baseline to compare against.
#
# Inputs:
#   behaviour/files/behaviour_mothers_predictions_total effects.rds
#     (behaviour/behaviour_mothers_analysis.R)
#   behaviour/files/behaviour_calves_predictions_total effects.rds
#     (behaviour/behaviour_calves_obspred.R)
#
# Run from the repo root (open srw_behaviour.Rproj).

source("plots/theme_paper.R")
library(ggh4x)           # facet_nested_wrap


# Data --------------------------------------------------------------------

dm <- readRDS("behaviour/files/behaviour_mothers_predictions_total effects.rds")
dc <- readRDS("behaviour/files/behaviour_calves_predictions_total effects.rds")

dm$mc <- "Mothers"
dc$mc <- "Calves"

# the calves' long-term and total effects are measured against 2013, not 1995,
# so they are not the same quantity as the mothers' and are left out
dc <- dc[dc$eff_type_term == "Short-term", ]

d <- rbind(as.data.frame(dm), as.data.frame(dc)[, names(dm)])
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))

dodge <- position_dodge(width = 1.2)


# Plot --------------------------------------------------------------------

fig5 <- ggplot(d, aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper,
                      colour = mc, fill = mc, shape = mc)) +
  geom_linerange(position = dodge, alpha = 0.6, linewidth = 0.4) +
  geom_point(position = dodge, size = 2.2, alpha = 0.9, stroke = 0.3,
             colour = "gray20") +
  facet_nested_wrap(vars(eff_type_predobs, eff_type_term), nrow = 1) +
  scale_colour_manual(values = colors_mc) +
  scale_fill_manual(values = colors_mc) +
  scale_shape_manual(values = c("Mothers" = 21, "Calves" = 24)) +
  scale_x_continuous(breaks = c(2000, 2010, 2020)) +
  ylim(0, max(d$d_upper) * 1.05) +
  ylab("Proportion of time in\naltered behaviour (δ)") +
  xlab("Year") +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
        legend.box.margin = margin(t = -3, unit = "mm"))

ggsave("plots/figure_5_effects_summary_v01.png", plot = fig5,
       height = 10, width = 18, units = "cm", dpi = 400)
