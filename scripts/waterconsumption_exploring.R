# This script aims to explore changes water consumption in orexin-cre mice

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()
library(ggpubr)
library(car)       # for Levene’s Test
library(tidyverse)
library(emmeans)

#format plot
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

#orexin-cre mice####
##water intake####
water <- read_csv("../data/WATER_CONSUMPTION.csv") %>% 
  mutate(DATE = lubridate::mdy(DATE)) %>% 
  filter(AIM=="UNCERTAINTY") %>% 
  group_by(ID) %>% 
  mutate(WATER_GR = (WATER_START_G - WATER_END_G)) %>% 
  mutate(delta_measurement = (DATE - lag(DATE)),
          corrected_waterintake_gr = (WATER_GR / as.numeric(delta_measurement))) %>% #DAILY WATER INTAKE
    mutate(GROUP = case_when(
    ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_UNCERT",
    ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_CERT",
    ID %in% c(307,316, 319) ~ "OREXIN_CRE_CONTROL_CERT"
  ))
n_distinct(water$ID) #here we know there is 14 animals

# Calculate mean water intake per date and group
mean_intake <- water %>%
  group_by(DATE, GROUP) %>%
  summarise(mean_water = mean(corrected_waterintake_gr, na.rm = TRUE), .groups = "drop") %>% 
  drop_na()

# Plot
plot <- ggplot(water, aes(x = DATE, y = corrected_waterintake_gr, group = ID)) +
  geom_line(aes(color = ID), alpha = 0.4) +  # Individual lines
  geom_point(aes(color = ID), size = 1.5, alpha = 0.6) +  # Individual points
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.5) +  # ID labels
  geom_line(data = mean_intake, aes(x = DATE, y = mean_water), 
            color = "black", size = 1.2, inherit.aes = FALSE) +  # Mean line
  labs(
    x = "Date",
    y = "Daily Water Intake (g)",
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text = element_text(family = "Helvetica", size = 13),
    axis.title = element_text(family = "Helvetica", size = 14),
    legend.position = "none"
  ) +
  facet_wrap(~GROUP)

plot

 
 ##one way anova ####
 # Check normality by group (Shapiro-Wilk test)
 water %>%
   group_by(GROUP) %>%
   summarise(p_value = shapiro.test(WATER_GR)$p.value) #all data normally distributed 
 # Check equal variances (Levene’s Test)
 leveneTest(WATER_GR ~ GROUP, data = water) #p = 0.6 
 #one way anoca
 anova_result <- aov(WATER_GR ~ GROUP, data = water)
 summary(anova_result)
 
 emmeans_result <- emmeans(anova_result, pairwise ~ GROUP)
 summary(emmeans_result)
 
 