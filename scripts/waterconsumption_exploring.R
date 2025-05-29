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
##body weight####
water <- read_csv("../data/WATER_CONSUMPTION.csv") %>% 
  filter(AIM=="UNCERTAINTY") %>% 
  mutate(WATER_GR = (WATER_START_G - WATER_END_G)) %>% 
  #mutate(delta_measurement = DATE - lag(DATE)),
         # corrected_waterintake_gr = WATER_GR / as.numeric(delta_measurement)) #DAILY WATER INTAKE
    mutate(GROUP = case_when(
    ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_UNCERT",
    ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_CERT",
    ID %in% c(307,316, 319) ~ "OREXIN_CRE_CONTROL_CERT"
  ))
n_distinct(water$ID) #here we know there is 14 animals

water %>%
  group_by(GROUP) %>% 
  summarise(mean = mean(WATER_GR, na.rm = TRUE)) 

#one way anova 
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


plot <- water  %>%
  ggplot(aes(GROUP, WATER_GR, group = ID)) +
  geom_point() +
  geom_line() +
   geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "group",
    y = "Water intake (grams)/2 days")+
  format.plot+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
 plot
 