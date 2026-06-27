library(tidyverse); theme_set(theme_bw())
library(viridis)

oc <- readRDS("attack/calves/calves_occurrence_prediction table.rds")
om <- readRDS("attack/mothers/mothers_occurrence_prediction table.rds")
cc <- readRDS("attack/calves/calves_count_prediction table.rds")
cm <- readRDS("attack/mothers/mothers_count_prediction table.rds")

nam <- names(oc)
om <- om[, nam]
cc <- cc[, nam]

names(cm)[names(cm) == "behavm_clust"] <- "behav"
cm <- cm[, nam]
names(cm)

d <- rbind(oc, om, cc, cm)
d$var <- factor(d$var, levels = c("Attack probability", "Attack number"))
d$mc <- factor(d$mc, levels = c("Mothers", "Calves"))
d$behav_raw2 <- factor(d$behav_raw, labels = c("R", "ST", "FT", "IPA"))
summary(d)

ggplot(d, 
       aes(x = behav_raw2, y = mu_mean, ymin = mu_lower, ymax = mu_upper,
           colour = exposure)) +
  geom_point(position = position_dodge2(width = 0.4), 
             size = 3, shape = 19, alpha = 0.7) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  geom_point(mapping = aes(x = behav_raw2, y = mu_obs,
                           colour = exposure),
             position = position_dodge2(width = 0.4), 
             size = 3, shape = 1, alpha = 1) +
  
  scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, name = NULL) +
  facet_grid(rows = vars(var), cols = vars(mc), scales = "free_y") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 15),
        legend.position = "bottom") +
  ylab(element_blank()) +  
  xlab("Whales' behaviour")

ggsave("attack/attack_figure_all.png", 
       width = 17, height = 13, units = "cm")
