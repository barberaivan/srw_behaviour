# Figure 4: behaviour under two attack scenarios, by period ----------------
#
# For every period, the behaviour probabilities of a whale that is attacked
# without pause from an undisturbed start ("Persistent attacks") and of a whale
# that is left alone from a fully disturbed start ("Attacks cessation"), over the
# first 70 minutes.
#
# Periods are 1995, 2004-2010, 2011-2013, 2014-2019 and 2020-2025 for mothers,
# and the last three for calves (their series starts in 2013). Single-year
# periods are plain model predictions; multi-year periods are the average of the
# years sampled in the period, averaged **within each posterior sample** — the
# posterior is summarised only at the end. That averaging is done by
# behaviour/behaviour_period_predictions.R, which writes the tables read here.
#
# The marker in the manuscript asks for two layouts, so two versions are saved:
#
#   version A (..._patchwork_vNN.png)  mothers and calves as two separate plots
#                                      merged with patchwork, tagged A and B,
#                                      with heights in proportion to their
#                                      number of periods (5 and 3)
#   version B (..._nested_vNN.png)     one plot, rows nested as
#                                      mother/calf x period
#
# Run from the repo root (open srw_behaviour.Rproj).

source("plots/theme_paper.R")
library(ggh4x)           # facet_nested
library(patchwork)


# Data --------------------------------------------------------------------

pm <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table_periods.rds")
pc <- readRDS("behaviour/files/behaviour_calves_predictions_attack scenarios table_periods.rds")

d <- rbind(pm, pc)
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))
d$scenario <- factor(d$scenario, levels = c("Attack", "After attack"),
                     labels = c("Persistent attacks", "Attacks cessation"))

# the z rows are not behaviour probabilities; the first 70 min is what the
# other figures of the paper show
d <- d[d$behav != "z" & d$min <= 70, ]

ymax <- max(d$prob_upper) * 1.02


# Building block ----------------------------------------------------------

panel_scenarios <- function(dat) {
  ggplot(dat,
         aes(x = min, y = prob_mean, ymin = prob_lower, ymax = prob_upper,
             colour = scenario, fill = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.5) +
    geom_ribbon(colour = NA, alpha = 0.18) +
    scale_colour_manual(values = colors_scenario) +
    scale_fill_manual(values = colors_scenario) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 0.8, by = 0.2),
                       expand = c(0.02, 0)) +
    scale_x_continuous(breaks = seq(0, 60, by = 20)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25),
          # the exposure labels are wider than the panels at size 10
          strip.text.x = element_text(size = 8.5,
                                      margin = margin(1.2, 0, 1.2, 0, "mm")),
          strip.text.y = element_text(angle = 270),
          legend.text = element_text(margin = margin(l = -1, r = 3, unit = "mm")),
          legend.box.margin = margin(t = -2, unit = "mm")) +
    ylab(NULL) + xlab(NULL)
}


# Version A: mothers and calves merged with patchwork ----------------------

pm_plot <- panel_scenarios(d[d$mc == "Mothers", ]) +
  facet_nested(rows = vars(period), cols = vars(behav_raw, exposure)) +
  ggtitle("A. Mothers")

pc_plot <- panel_scenarios(d[d$mc == "Calves", ]) +
  facet_nested(rows = vars(period), cols = vars(behav_raw, exposure)) +
  ggtitle("B. Calves") +
  # the column strips are already named on top of the mothers block
  theme(strip.text.x = element_blank(),
        strip.background.x = element_blank())

fig4a <- (pm_plot / pc_plot) +
  plot_layout(heights = c(5, 3), guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, hjust = 0, vjust = 1))

fig4a <- wrap_elements(fig4a) +
  labs(tag = "Behaviour probability") +
  theme(plot.tag = element_text(size = 12, angle = 90),
        plot.tag.position = "left")

fig4a <- fig4a /
  grid::textGrob("Time since beginning or cessation of attacks (min)",
                 gp = grid::gpar(fontsize = 12)) +
  plot_layout(heights = c(50, 1))

ggsave("plots/figure_4_attack_scenarios_patchwork_v01.png", plot = fig4a,
       height = 22, width = 19, units = "cm", dpi = 400)


# Version B: one plot, rows nested as mother/calf x period -----------------

fig4b <- panel_scenarios(d) +
  facet_nested(rows = vars(mc, period), cols = vars(behav_raw, exposure)) +
  ylab("Behaviour probability") +
  xlab("Time since beginning or cessation of attacks (min)")

ggsave("plots/figure_4_attack_scenarios_nested_v01.png", plot = fig4b,
       height = 22, width = 19, units = "cm", dpi = 400)
