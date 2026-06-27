# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(ggh4x)           # facet_nested
library(gridExtra)       # grid.arrange (order plots)
library(grid)            # textGrob
library(egg)             # ggarrange

# Custom ggplot theme -----------------------------------------------------
# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() { 
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),               
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),               
      
      # para separar el eje y de los nros
      axis.title.y = element_text(             
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),
      
      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size
      
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),
      
      strip.text = element_text(size = 10, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      
      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())


theme_mine2 <- function() { 
  font <- "Arial"
  theme_mine() %+replace%    #replace elements we want to change
    
    theme(
      strip.text = element_text(size = 10, family = font, color = "black"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "white", color = "white"),
    )
}

# Function to extract legend and xlab  ---------------------------------------
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g_xlab <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "xlab")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Figure 1: attack scenario for mothers and calves ------------------------


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
             behav != "z")

d$scenario <- factor(d$scenario, levels = levels(dm$scenario),
                     labels = c("Persistent attacks", "Attacks cessation"))

colors <- viridis(2, end = 0.5)


# Set some theming options, we can use `element_blank()`
a <- element_rect(fill = "gray10")
b <- element_rect(fill = "darkgray")
backgrounds <- list(a, b, a, a, a, b)

a_t <- element_text(colour = "white")
b_t <- element_text(colour = "black")
texts <- list(a_t, b_t, a_t, a_t, a_t, b_t)
# Or we could use `NULL` to use the global theme
# texts <- list(element_text(colour = "red"), NULL, element_text(face = "bold"))

data_rect <- data.frame(mc = factor("Calves", levels = c("Mothers", "Calves")),
                        year = factor("2013-2018", levels = c("1995", "2004", "2013-2018")))

# plot without separation between mothers and calves
ggplot(d[filt, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = scenario)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.2) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc, year),
               nest_line = element_line(colour = "white", size = 0.4),
               strip = strip_nested(
                 background_x = backgrounds,
                 text_x = texts,
                 clip = "off"
               )) +
  
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        # panel.background = element_rect(fill = "black"),
        strip.text.y = element_text(angle = 270),
        legend.text = element_text(margin = margin(l = -1, r = 3, unit = "mm")),
        legend.box.margin = margin(t = -2, unit = "mm")) +
  ylab("Behaviour probability") +
  xlab("Time since beginning or cessation of attacks (min)")


# Lo mismo pero uniendo con ggarrange mothers y calves, para separarlos más

filt_m <- with(d, min <= 70 & behav != "z" & mc == "Mothers")
filt_c <- with(d, min <= 70 & behav != "z" & mc == "Calves")

pm <- 
ggplot(d[filt_m, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = scenario)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc, year),
               nest_line = element_line(colour = "white", size = 0.4),
               strip = strip_nested(
                 text_y = list(element_blank(), element_text(colour = "white")),
                 background_y = list(element_blank(), element_rect(fill = "white", color = "white")), 
                 by_layer_y = TRUE
               )) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2),
                     limits = c(0, 0.8),
                     expand = c(0.05, 0)) +
  
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        strip.text.y = element_text(angle = 270, size = 0.001),
        legend.text = element_text(margin = margin(l = -1, r = 3, unit = "mm"),
                                   size = 10),
        legend.box.margin = margin(t = -2, unit = "mm")) +
  ylab("Behaviour probability") +
  xlab("Time since beginning or cessation of attacks (min)")
# pm

pc <- 
  ggplot(d[filt_c, ], 
         aes(x = min, y = prob_mean, ymin = prob_lower, 
             ymax = prob_upper, 
             colour = scenario, fill = scenario,
             linetype = scenario)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc, year),
               nest_line = element_line(colour = "white", size = 0.4)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2),
                     limits = c(0, 0.8),
                     expand = c(0.05, 0)) +
  
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        strip.text.y = element_text(angle = 270),
        legend.text = element_text(margin = margin(l = -1, r = 3, unit = "mm")),
        legend.box.margin = margin(t = -2, unit = "mm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab(NULL)
pc

# Unir con egg::ggarrange sin legend ni eje x
joined <- egg::ggarrange(
  pm + theme(axis.title.x = element_blank(), 
             legend.position = "none",
             plot.margin = margin(2, 0, 0, 2, unit = "mm")), 
  pc + theme(axis.title.x = element_blank(), 
             legend.position = "none",
             plot.margin = margin(2, 2, 0, 0, unit = "mm")), 
  nrow = 1, widths = c(3, 1)
  
  # La label queda fuera del plot al guardar, no sé acomodarlo.
  # labels = c("        A", "B"),
  # label.args = list(gp = grid::gpar(face = "plain"),
  #                   hjust = -0.3)
)


# agregar texto y legenda usando grid.arrange()
m_legend <- g_legend(pm)
x_label <- textGrob("Time since beginning or cessation of attacks (min)", 
                    gp = gpar(fontsize = 12),
                    hjust = 0.48)#, vjust = -1.8

joined2 <- grid.arrange(
  joined, x_label, 
  arrangeGrob(ggplot() + theme_void(), m_legend,
              nrow = 1, widths = c(0.5, 10)),
  nrow = 3, heights = c(20, 0.7, 1)
)

# ggsave("plots/behaviour_both_predictions_attack scenarios.png", 
#        plot = joined2,
#        height = 23, width = 17, units = "cm")



# Figure S1 mothers attack scenario for all years --------------------------

predictions <- readRDS("behaviour/files/behaviour_mothers_predictions_attack scenarios table object.rds")

# Plots
predictions$scenario <- factor(predictions$scenario, 
                               levels = levels(predictions$scenario),
                               labels = c("Persistent attacks", 
                                          "Attacks cessation"))

colors <- viridis(2, end = 0.5)

filt1 <- with(predictions, 
              min <= 70 & 
              behav != "z")

ggplot(predictions[filt1, ], 
       aes(x = min, y = prob_mean, ymin = prob_lower, 
           ymax = prob_upper, 
           colour = scenario, fill = scenario,
           linetype = exposure)) + 
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.15) +
  #facet_grid(rows = vars(year), cols = vars(behav_raw)) +
  facet_nested(cols = vars(behav_raw, exposure), rows = vars(year)) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270),
        legend.text = element_text(margin = margin(r = 3, l = -1, unit = "mm"))) +
  ylab("Behaviour probability") +
  xlab("Time since beginning or cessation of attacks (min)") #+

# ggsave("plots/behaviour_mothers_predictions_attack scenarios all years.png",
#        height = 35, width = 20, units = "cm")


# Figure 2: steady states by year mothers ---------------------------------

l1 <- readRDS("behaviour/files/behaviour_mothers_predictions disturbed and undisturbed by year.rds")
ppp <- l1$plot_data
ghosts <- l1$ghosts

colors3 <- c(viridis(2, end = 0.5), 
             viridis(3, end = 0.93)[3])

ppp$scenario <- factor(ppp$scenario, 
                       levels = c("Disturbed", "Undisturbed", "Observed"))

ggplot(ppp,#[filt4, ], 
       aes(x = year_num, 
           y = prob_mean, ymin = prob_lower, ymax = prob_upper, 
           shape = scenario, fill = scenario)) + 
  geom_point(position = position_dodge2(width = 0.3), size = 3) +
  geom_linerange(position = position_dodge2(width = 0.3), show.legend = FALSE) +
  geom_point(data = ghosts, mapping = aes(x = year_num, y = prob_lower), alpha = 0) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(period), 
               #scale = "free", space = "free") +
               scales = "free", space = "free_x",
               nest_line = element_line(colour = "white", size = 0.4)
  )+
  scale_fill_manual(values = colors3) +
  scale_shape_manual(values = c(23, 22, 21)) +
  # guides(shape = "none") +
  ylab("Behaviour probability") +
  xlab("Year") +
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270),
        legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
        legend.box.margin = margin(t = 0, unit = "mm")) +
        # esto mejora la seaparación entre texto y key de la legend
  scale_x_continuous(breaks = c(1995, 2005, 2008, 2012, 2016)) 

# ggsave("plots/behaviour_mothers_ts_free_y.png", height = 21, width = 17, unit = "cm")


# Figure 3: effects on mothers' and calves' behaviour ---------------------


delta_wide <- readRDS("behaviour/files/behaviour_mothers_predictions_total effects.rds")

# Plot 
# ggplot(delta_wide, 
#        aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper)) +
#   geom_point(size = 2) + 
#   geom_linerange(alpha = 0.6) +
#   facet_nested_wrap(vars(eff_type_predobs, eff_type_term), nrow = 1) +
#   ylim(0, max(delta_wide$d_upper * 1.05)) + 
#   theme(panel.grid.minor = element_blank()) + 
#   ggtitle("Attacks effects on mothers behaviour") + 
#   ylab("Proportion of time in altered behaviour (d)") +
#   xlab("Year")


ggplot(delta_wide, 
       aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper,
           colour = eff_type_predobs, shape = eff_type_predobs,
           fill = eff_type_predobs)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_linerange(alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.45, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.45, direction = -1) +
  scale_shape_manual(values = c(22, 21)) +
  facet_wrap(vars(eff_type_term), nrow = 1) +
  ylim(0, max(delta_wide$d_upper * 1.05)) + 
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
        legend.box.margin = margin(-10,0,0,0)) + 
  ylab("Proportion of time in\naltered behaviour (d)") +
  xlab("Year")

# effects on calves' behaviour

delta_wide_c <- readRDS("behaviour/files/behaviour_calves_predictions_total effects.rds")
names(delta_wide_c)
names(delta_wide)
delta_wide_c$eff_type_term <- "Short-term"
delta_wide_c$eff_type_predobs <- as.character(delta_wide_c$effect_type)
delta_wide_c$eff_type_predobs[delta_wide_c$eff_type_predobs == "Short-term"] <- "Observed"

# mothers and calves effects
delta_wide_c$mc <- "Calves"
delta_wide$mc <- "Mothers"

dbind <- rbind(delta_wide, delta_wide_c)
dbind$mc <- factor(dbind$mc, levels = c("Mothers", "Calves"))

ggplot(dbind, 
       aes(x = year, y = d_mean, ymin = d_lower, ymax = d_upper,
           colour = eff_type_predobs, shape = eff_type_predobs,
           fill = eff_type_predobs)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_linerange(alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.45, direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.45, direction = -1) +
  scale_shape_manual(values = c(22, 21)) +
  facet_grid(cols = vars(eff_type_term), rows = vars(mc)) +
  ylim(0, max(delta_wide$d_upper * 1.05)) + 
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(margin = margin(r = 3, l = -2, unit = "mm")),
        legend.box.margin = margin(-10,0,0,0),
        strip.text.y = element_text(angle = 270)) + 
  ylab("Proportion of time in altered behaviour (d)") +
  xlab("Year")


# ggsave("plots/behaviour_mothers and calves_ts effects total.png",
#        height = 11, width = 17, units = "cm")



# Figure S2: joint behaviour mothers and calves ---------------------------

zdyn <- rbind(
  readRDS("behaviour/files/behaviour_mothers_predictions_attack scenario table_bmc.rds"),
  readRDS("behaviour/files/behaviour_calves_predictions_attack scenario table_bmc.rds")
)
zdyn$mc <- factor(zdyn$mc, levels = c("Mothers", "Calves"))
# zdyn$scenario <- factor(zdyn$scenario, levels = levels(zdyn$scenario),
#                         labels = c("Both", "Mothers", "Calves"))

ggplot(zdyn[complete.cases(zdyn) & zdyn$min <= 60, ],
       aes(x = min, y = p_mean, ymin = p_lower, ymax = p_upper,
           colour = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(alpha = 0.1, color = NA) +
  facet_nested(rows = vars(behav_raw, exposure), cols = vars(mc),
               nest_line = element_line(colour = "white", size = 0.4)) +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.7,
                      name = "Attacks on") +
  scale_fill_viridis(discrete = TRUE, option = "A", end = 0.7,
                     name = "Attacks on") +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270),
        legend.position = "right",
        legend.title = element_text()) +
  labs(colour = "Attack scenario", fill = "Attack scenario") + 
  ylab("Behaviour probability") + 
  xlab("Time since beginning of attacks (min)")

ggsave("plots/behaviour_both_predictions_attack scenarios_bmc.png", 
       height = 23, width = 11, units = "cm")




# Lo que sigue: -----------------------------------------------------------

# - Behav -

# Casi mismos plots para calves de behav,
# Ajuste de los modelos,
# Comportamiento conjunto de madres y crías (casi todo hecho en plots_script_mother calf.R)

# - Attacks -
# Resumir todo en un solo plot:
#   (attack prob, attack number) ~ behaviour,
#   con variables en filas y mc en columnas. 
#     (usar 2 facet wraps y luego rbindear con ggarrange)
#     Poner media estimada e histograma en attack number,
#     y proporción obs y estimada en attack prob.
#   
# Suplementario: cambio en attack p en función del behav entre años 
# (hacer plot similar a ts de behav para madres)
# Quizás haga falta correr ese modelo de ataque again
 
# Quizás se pueda simplificar todo e ignorar la recursividad en attack p.



# Attacks plots ------------------------------------------------------------

# mothers data (Using spline model for mothers -occurrence-)

mo_list <- readRDS("attack/mothers_spline years/model_gaviots_occurrence_02 predictions and data for plots.rds")
names(mo_list)

mc_list <- readRDS("attack/mothers_spline years/model_gaviots_count_02 predictions for plot.rds")
names(mc_list)

mo_pred <- mo_list$pred_steady_avg
mo_data <- mo_list$data_avg

mc_pred <- mc_list$count_means
mc_data <- mc_list$count_data
# aggregate count data 
mc_data_agg1 <- aggregate(nam ~ follow + year + behav_raw + exposure, 
                          mc_data, mean)
mc_data_agg2 <- aggregate(nam ~ year + behav_raw + exposure, 
                          mc_data_agg1, mean)
mc_data_agg3 <- aggregate(nam ~ behav_raw + exposure, 
                          mc_data_agg2, mean)


# calves
cc_data <- readRDS("models/attack/calves/calves_count_data_subset.rds")
cc_pred <- readRDS("attack/calves/calves_count_prediction table.rds")

co_data <- readRDS("models/attack/calves/calves_occurrence_data_subset.rds")
co_pred <- readRDS("attack/calves/calves_occurrence_prediction table.rds")

# aggregate occurrence data 
co_data_agg1 <- aggregate(ac_marg ~ follow_id + year_id + behav_raw + exposure, 
                          co_data, mean)
co_data_agg2 <- aggregate(ac_marg ~ year_id + behav_raw + exposure, 
                          co_data_agg1, mean)
co_data_agg3 <- aggregate(ac_marg ~ behav_raw + exposure, 
                          co_data_agg2, mean)

# aggregate count data 
cc_data_agg1 <- aggregate(nac ~ follow + year + behav_raw + exposure, 
                          cc_data, mean)
cc_data_agg2 <- aggregate(nac ~ year + behav_raw + exposure, 
                          cc_data_agg1, mean)
cc_data_agg3 <- aggregate(nac ~ behav_raw + exposure, 
                          cc_data_agg2, mean)


# merge data means
names(mo_data)[which(names(mo_data) == "am_marg")] <- "mean"
names(co_data_agg3)[which(names(co_data_agg3) == "ac_marg")] <- "mean"
mo_data$mc <- "Mothers"
co_data_agg3$mc <- "Calves"
o_data <- rbind(mo_data, co_data_agg3)
o_data$variable <- "Attack probability"

names(mc_data_agg3)[which(names(mc_data_agg3) == "nam")] <- "mean"
names(cc_data_agg3)[which(names(cc_data_agg3) == "nac")] <- "mean"
mc_data_agg3$mc <- "Mothers"
cc_data_agg3$mc <- "Calves"
c_data <- rbind(mc_data_agg3, cc_data_agg3)
c_data$variable <- "Attacks number"

at_data <- rbind(o_data, c_data)


# merge predictions
names(mo_pred)[grep("p_", names(mo_pred))] <- c("mu_lower", "mu_mean", "mu_upper")
mo_pred$mc <- "Mothers"
o_pred <- rbind(mo_pred, co_pred[, names(mo_pred)])
o_pred$variable <- "Attack probability"

names(mc_pred)[grep("p_", names(mo_pred))] <- c("mu_lower", "mu_mean", "mu_upper")
mc_pred <- mc_pred[, names(mo_pred)[1:5]]
cc_pred <- cc_pred[, names(mc_pred)]
mc_pred$mc <- "Mothers"
cc_pred$mc <- "Calves"
c_pred <- rbind(mc_pred, cc_pred)
c_pred$variable <- "Attacks number"

at_pred <- rbind(o_pred, c_pred)

# relevel factors
at_pred$mc <- factor(at_pred$mc, levels = c("Mothers", "Calves"))
at_pred$variable <-  factor(at_pred$variable, 
                            levels = c("Attack probability", 
                                       "Attacks number"))
at_data$mc <- factor(at_data$mc, levels = c("Mothers", "Calves"))
at_data$variable <-  factor(at_data$variable, 
                            levels = c("Attack probability", 
                                       "Attacks number"))

# merge all
at_data2 <- at_data
names(at_data2)[which(names(at_data) == "mean")] <- "mu_mean"
at_data2$mu_lower <- NA
at_data2$mu_upper <- NA
at_data2$data_type <- "Observed"

at_pred$data_type <- "Predicted"
at_data2 <- at_data2[, names(at_pred)]

at_all <- rbind(at_data2, at_pred)
at_all$data_type <- factor(at_all$data_type, levels = c("Observed", "Predicted"))

# plot
ggplot(at_all, mapping = aes(x = behav_raw, 
                   y = mu_mean, ymin = mu_lower, ymax = mu_upper,
                   colour = exposure,
                   shape = data_type)) +
  # prediction
  geom_linerange(data = at_all[at_all$data_type == "Predicted", ],
                 position = position_dodge2(width = 0.4),
                 size = 3, alpha = 0.5, 
                 show.legend = F) +
  geom_point(data = at_all[at_all$data_type == "Predicted", ],
             position = position_dodge2(width = 0.4), size = 3,
             stroke = 0.5) +

  # data
  geom_point(data = at_all[at_all$data_type == "Observed", ],
             position = position_dodge2(width = 0.4), size = 3) +
  
  # scale_shape_manual(values = c(17, 22)) + 
  scale_shape_manual(values = c(20, 23)) + 
  scale_color_viridis(option = "A", end = 0.5, discrete = TRUE, name = "Position") +
  facet_grid(rows = vars(variable), cols = vars(mc), scales = "free_y",
             switch = "y") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8),
        axis.title.x = element_text(size = 12, margin = margin(t = 3, unit = "mm")),
        
        text = element_text(size = 15),
        strip.text.y = element_text(angle = 270, colour = "black", size = 12),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        legend.box.margin = margin(t = -3, unit = "mm"),
        legend.text = element_text(margin = margin(l = -3, r = 5, unit = "mm"), 
                                   size = 10),
        legend.title = element_text(size = 11)) +
  ylab(NULL) +  
  xlab("Whales' behaviour") +
  guides(colour = guide_legend(order = 1, 
                               override.aes = list(shape = 15, size = 4),
                               title = "Exposure",
                               title.position = "top",
                               title.vjust = -1,
                               title.hjust = 0.4), 
         shape = guide_legend(order = 2, 
                              title = "Data type",
                              title.position = "top",
                              title.vjust = -1))


ggsave("plots/attack figure.png", width = 17, height = 12, units = "cm")
