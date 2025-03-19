# We aim to test the rate of change of FI over 12 weeks of LFD in NZO mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models

fi_data <- read_csv("~/Documents/Github/data/scripts/SSIB_2025_Lauren/fi_data.csv") %>% 
filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>% 
  group_by(ID) %>%
  mutate(start_date = ifelse(COMMENTS == "DAY_7_SABLE_DAY_1_LFD", DATE, NA)) %>%
  fill(start_date, .direction = "down") %>%  # Propagate the date down within each ID group
  filter(DATE >= start_date) %>%
  select(-start_date) %>%  # Remove the helper column
  drop_na(corrected_intake_kcal) %>% 
  select(ID,DATE,corrected_intake_kcal) %>% 
  filter(is.finite(corrected_intake_kcal)) %>%  #to eliminate inf values 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTION PERIOD FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal)) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols="cumsum_kcal") %>% 
  select(
    ID, time, name, value
  ) 
# Fit a linear mixed model: cumulative intake ~ time with random effect for ID
model <- lmer(value ~ time + (1 | ID), data = fi_data)
# Summary of the model
summary(model)

#The model estimated the rate of food intake as 156.99 kcal/day across all animals. 
#Given that we analyzed 23 animals, this corresponds to an average increase of 6.82 kcal/day per mouse (SE = 0.76, t = 205.9).

