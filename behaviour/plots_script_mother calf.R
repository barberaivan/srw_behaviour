library(tidyverse); theme_set(theme_bw())
library(viridis)
library(ggh4x)

# Attack scenario prediction for mother and calves ------------------------

dm <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table object.rds")
dm_avg <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table object_avg.rds")
dc_avg <- readRDS("behaviour/files/behaviour_calves_predictions_attack scenarios table object_avg.rds")

years_single <- c(1995, 2004)
dm_avg$year <- dc_avg$year <- "2013-2018"

# reorder columns:
dm_avg <- dm_avg[, names(dm)]
dc_avg <- dc_avg[, names(dm)]

# mother calf column
dm$mc <- "Mothers"
dm_avg$mc <- "Mothers"
dc_avg$mc <- "Calves"

d <- rbind(dm[dm$year %in% years_single, ],
           dm_avg, dc_avg)
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))

filt <- with(d, 
  min <= 70 & 
  behav != "z"
)

colors <- viridis(2, end = 0.5)

# Plot
ggplot(d[filt, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = scenario)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_nested(cols = vars(behav_raw, exposure), rows = vars(mc, year)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behaviour probability") +
  xlab("Time (min)")
# ggsave("behaviour/figures/behaviour_both_predictions_attack scenarios_a.png", 
#        height = 14, width = 23, units = "cm")

ggplot(d[filt, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = scenario)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc, year)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank()) +
  ylab("Behaviour probability") +
  xlab("Time (min)")
# ggsave("behaviour/figures/behaviour_both_predictions_attack scenarios_b.png", 
#        height = 23, width = 17, units = "cm")


# Total effects for mothers and calves ------------------------------------

em <- readRDS("behaviour/files/behaviour_mothers_predictions_total effects.rds")
ec <- readRDS("behaviour/files/behaviour_calves_predictions_total effects.rds")

em$mc <- "Mothers"
ec$mc <- "Calves"

ec <- ec[, names(em)]

eff <- rbind(em, ec)
eff$mc <- factor(eff$mc, levels = c("Mothers", "Calves"))

# data for attacks:
bdata <- readRDS("models/behaviour/behaviour_mothers_data with imputed behaviours.rds")
# (use last_not)

# Plot 
ggplot(eff, 
       aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper,
           colour = mc)) +
  geom_point(size = 2,
             position = position_dodge2(width = 0.9)) + 
  geom_linerange(alpha = 0.6,
                 position = position_dodge2(width = 0.9)) +
  facet_wrap(vars(effect_type), nrow = 2) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") + 
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.7) +
  ylab("Proportion of time in altered behaviour (d)") +
  xlab("Year") 
# ggsave("behaviour/figures/behaviour_both_effects year.png",
#        height = 15, width = 17, units = "cm")


# Observed effects as a function of attacks

names(bdata)
atd <- do.call("data.frame", 
               aggregate(cbind(a, am_marg, ac_marg) ~ year,
                         data = bdata[bdata$last_not, ],
                         FUN = mean))

names(atd) <- c("year", "On mothers or calves", "On mothers", "On calves")
# merge
effat <- left_join(eff, atd, by = "year")
efflong <- pivot_longer(effat, which(names(effat) %in% 
                                     names(atd)[2:4]),
                        names_to = "attack_var", values_to = "attack")
#efflong <- efflong[efflong$effect_type %in% c("Short-term", "Potential"), ]
efflong <- efflong[efflong$effect_type %in% c("Short-term"), ]
efflong$attack_var <- factor(efflong$attack_var, levels = names(atd)[2:4])


ggplot(efflong, 
       aes(x = attack, y = d_mean, ymin = d_lower, ymax = d_upper)) +
  geom_smooth(method = "lm", se = F, lty = 2, size = 0.5, colour = "gray") +
  geom_point(size = 2) + 
  geom_linerange(alpha = 0.6) +
  facet_grid(cols = vars(mc), rows = vars(attack_var)) +
  #facet_wrap(vars(attack_var), nrow = 1, scales = "free_x") +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ylim(0, 0.3) +
  ylab("Proportion of time in altered behaviour (d)") +
  xlab("Attack frequency") 
# ggsave("behaviour/figures/behaviour_both_effects short-term attack.png",
#        height = 20, width = 17, units = "cm")


# z dynamics for mothers and calves ---------------------------------------

zdyn <- rbind(
  readRDS("behaviour/files/behaviour_mothers_predictions_attack scenario table_bmc.rds"),
  readRDS("behaviour/files/behaviour_calves_predictions_attack scenario table_bmc.rds")
)
zdyn$mc <- factor(zdyn$mc, levels = c("Mothers", "Calves"))

ggplot(zdyn[complete.cases(zdyn) & zdyn$min <= 60, ],
       aes(x = min, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(alpha = 0.1, color = NA) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc)) +
  #facet_nested(cols = vars(behav_raw, exposure), rows = vars(mc)) +  
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.7) +
  theme(panel.grid.minor = element_blank()) +
  labs(colour = "Attack scenario", fill = "Attack scenario")

# ggsave("behaviour/figures/behaviour_both_predictions_attack scenarios_bmc.png", 
#        height = 23, width = 15, units = "cm")



