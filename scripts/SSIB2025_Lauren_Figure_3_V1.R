# We aim to evaluate cumulative food intake in NZO mice after 12 weeks of LFD

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models
library(ggpubr)  

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

#FI data import####

FI_data <-read_csv("~/Documents/GitHub/data/data/FI.csv") #data import

FI_data <- FI_data %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
 filter(DATE <= "2025-02-24") %>% #DAY 1 OF RESTRICTION 
  filter(DIET == "LFD") %>%
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(ID != 3715)%>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  drop_na(INTAKE_GR) %>% 
  filter(corrected_intake_kcal <= 50)%>% # we have two animals that ate more than 75 kcal in one day so we believe is a systematic error 
  mutate(DATE_rel = DATE - first(DATE),
         cumulativeFI = cumsum(corrected_intake_kcal),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  mutate(relative_weeks = as.integer(DATE_rel / 7) + 1) %>% 
  filter(DATE_rel <= 100)


plot_1 <- FI_data %>%
  ggplot(aes(x = DATE_rel, y =cumulativeFI)) +
  geom_point(aes(group = ID), alpha = 1) +
  geom_line(aes(group = ID), alpha = 1) +
  labs(
    x = "Weeks",
    y = " Food intake (kcal)"
  ) #+
 # geom_text(aes(label = ID))

plot_1

