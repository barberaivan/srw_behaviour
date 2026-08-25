# Figure 6: attacks, behaviour and calf mortality through the years ---------
#
# Three stacked panels sharing the year axis:
#
#   (1) gull attacks: the proportion of intervals of the year in which the
#       mother, the calf or both were attacked (left axis), and how hard the
#       attacks were — the mean number of attacks per interval, on the mother and
#       on the calf separately (right axis);
#   (2) behaviour: the proportion of time spent resting and travelling (slow plus
#       fast travel), each split by exposure, which gives four curves. In-place
#       activity is left out, so the four curves do not add up to one;
#   (3) calf mortality: EMPTY PLACEHOLDER, see the note below.
#
# Two versions are saved, differing only in panel 2:
#
#   version A (..._mothers_vNN.png)  panel 2 with mothers only
#   version B (..._both_vNN.png)     panel 2 with mothers and calves, keyed by
#                                    point shape and line type
#
# ---------------------------------------------------------------------------
# NOTE FOR MERI (panel 3): which calf mortality series should go here — the dead
# calves counted by the Southern Right Whale Health Monitoring Program, as in
# Piotto et al. (2024)? And do we have mortality data up to 2025, or does the
# series still stop earlier? The panel is drawn empty until we know; the moment
# the numbers are available, drop them in the `mortality` data frame below.
# ---------------------------------------------------------------------------
#
# Attacks are counted on the mothers' dataset, which covers the whole 1995-2025
# series (the calves' behaviour, not the attacks on calves, is what starts in
# 2013). Behaviour proportions use the model-imputed behaviour column, so every
# interval contributes.
#
# Run from the repo root (open srw_behaviour.Rproj).

source("plots/theme_paper.R")
library(patchwork)


# Data --------------------------------------------------------------------

bm <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")
bc <- readRDS("models/behaviour/behaviour_calves_data with imputed behaviours.rds")

year_breaks <- c(1995, 2000, 2005, 2010, 2015, 2020, 2025)
year_range <- range(bm$year)

# Behaviour was not recorded between 1996 and 2003, nor in 2014, 2020, 2022 and
# 2024. Adding the missing years as NA breaks the lines at the gaps instead of
# drawing a segment across a decade.
add_year_gaps <- function(df, by, from = min(df$year)) {
  grid <- expand.grid(year = from:max(df$year),
                      stringsAsFactors = FALSE)
  keys <- unique(df[, by, drop = FALSE])
  grid <- merge(grid, keys, by = NULL)
  out <- merge(grid, df, by = c("year", by), all.x = TRUE)
  return(out[order(out$year), ])
}


# Panel 1: attacks ---------------------------------------------------------

# a       = attack on the mother, on the calf or on both, in the interval
# nam/nac = number of attacks received by the mother / by the calf
att <- aggregate(cbind(a, nam, nac) ~ year, bm, mean)

occurrence <- data.frame(year = att$year, value = att$a,
                         variable = "Attack occurrence (either)")

# the two intensity series are drawn against a second axis; `k` maps them onto
# the occurrence scale and the secondary axis undoes it
k <- max(c(att$nam, att$nac)) / max(att$a)

intensity <- rbind(
  data.frame(year = att$year, value = att$nam / k, variable = "Attacks on mother"),
  data.frame(year = att$year, value = att$nac / k, variable = "Attacks on calf")
)

att_long <- rbind(occurrence, intensity)
att_long <- add_year_gaps(att_long, "variable")
att_long$variable <- factor(att_long$variable,
                            levels = c("Attack occurrence (either)",
                                       "Attacks on mother", "Attacks on calf"))

colors_att <- c("Attack occurrence (either)" = "gray15",
                "Attacks on mother" = unname(colors_mc["Mothers"]),
                "Attacks on calf" = unname(colors_mc["Calves"]))

p1 <- ggplot(att_long, aes(x = year, y = value, colour = variable,
                           shape = variable, linetype = variable)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  scale_colour_manual(values = colors_att) +
  scale_shape_manual(values = c(21, 16, 17)) +
  scale_linetype_manual(values = c("solid", "22", "22")) +
  scale_x_continuous(limits = year_range, breaks = year_breaks) +
  scale_y_continuous(
    name = "Proportion of intervals\nwith attacks",
    sec.axis = sec_axis(~ . * k, name = "Attacks per interval")
  ) +
  xlab(NULL) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.width = unit(8, "mm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# Panel 2: grouped behaviour classes ---------------------------------------

# rest and travel (slow + fast), each split by exposure; in-place activity is
# not shown
group_behaviour <- function(df, behav_col, mc) {
  b <- as.character(df[[behav_col]])
  grp <- ifelse(grepl("^R_", b), "Rest",
                ifelse(grepl("^ST_|^FT_", b), "Travel", NA))
  exp <- ifelse(grepl("_UW$", b), "Under water", "Surface")
  keep <- !is.na(grp)
  tab <- table(df$year[keep], paste(grp[keep], exp[keep], sep = ", "))
  # denominator: all intervals of the year, in-place activity included
  n_year <- as.numeric(table(df$year))
  out <- as.data.frame(tab / n_year, stringsAsFactors = FALSE)
  names(out) <- c("year", "group", "prop")
  out$year <- as.numeric(out$year)
  out$mc <- mc
  return(out)
}

beh <- rbind(add_year_gaps(group_behaviour(bm, "behavm_full", "Mothers"),
                           c("group", "mc")),
             add_year_gaps(group_behaviour(bc, "behavc_full", "Calves"),
                           c("group", "mc")))
beh$mc <- factor(beh$mc, levels = c("Mothers", "Calves"))
beh$group <- factor(beh$group, levels = c("Rest, Surface", "Rest, Under water",
                                          "Travel, Surface", "Travel, Under water"))

colors_group <- setNames(viridis(4, end = 0.9), levels(beh$group))

panel_behaviour <- function(dat, with_calves) {
  p <- ggplot(dat, aes(x = year, y = prop, colour = group, fill = group)) +
    scale_colour_manual(values = colors_group) +
    scale_fill_manual(values = colors_group) +
    scale_x_continuous(limits = year_range, breaks = year_breaks) +
    ylab("Proportion of time") +
    xlab(NULL) +
    theme(panel.grid.minor = element_blank(),
          legend.key.width = unit(8, "mm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  if (with_calves) {
    p <- p +
      geom_line(aes(linetype = mc), linewidth = 0.4) +
      geom_point(aes(shape = mc), size = 2) +
      scale_linetype_manual(values = c("Mothers" = "solid", "Calves" = "22")) +
      scale_shape_manual(values = c("Mothers" = 16, "Calves" = 17)) +
      guides(colour = guide_legend(order = 1), fill = "none",
             linetype = guide_legend(order = 2),
             shape = guide_legend(order = 2))
  } else {
    p <- p +
      geom_line(linewidth = 0.4) +
      geom_point(size = 2) +
      guides(fill = "none")
  }
  return(p)
}

p2_mothers <- panel_behaviour(beh[beh$mc == "Mothers", ], with_calves = FALSE)
p2_both <- panel_behaviour(beh, with_calves = TRUE)


# Panel 3: calf mortality (placeholder) ------------------------------------

# Empty on purpose: the series to plot is still to be decided (see the note at
# the top). Fill `mortality` with columns year and value and the panel draws
# itself.
mortality <- data.frame(year = numeric(0), value = numeric(0))

p3 <- ggplot(mortality, aes(x = year, y = value)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  annotate("text", x = mean(year_range), y = 0.5,
           label = "series to be decided: see the note for Meri",
           colour = "#008080", size = 3) +
  scale_x_continuous(limits = year_range, breaks = year_breaks) +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Calf mortality") +
  xlab("Year") +
  theme(panel.grid.minor = element_blank())


# Assemble ------------------------------------------------------------------

assemble <- function(p2) {
  (p1 / p2 / p3) +
    plot_layout(heights = c(1, 1, 0.8), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right",
          legend.justification = "left",
          plot.tag = element_text(size = 11, face = "plain"),
          plot.margin = margin(1, 2, 1, 2, unit = "mm"))
}

ggsave("plots/figure_6_attacks_behaviour_mortality_mothers_v02.png",
       plot = assemble(p2_mothers),
       height = 19, width = 18, units = "cm", dpi = 400)

ggsave("plots/figure_6_attacks_behaviour_mortality_both_v02.png",
       plot = assemble(p2_both),
       height = 19, width = 18, units = "cm", dpi = 400)
