#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)

#bw data import####
BW_data <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = BW - first(BW),
    date_rel = DATE - first(DATE)
  ) %>% 
  filter(COHORT==0) %>% #We just want clarity guys
  filter(!ID %in% c(7875, 7874)) #we want to eliminate the ID were used for antibody titration

ggplot(BW_data, aes(DATE, BW)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA) +
  geom_smooth() +
  theme_minimal()

#fI data import####
FI_data <- read_csv("~/Documents/GitHub/data/data/FI.csv") %>% 
  filter(COHORT==0) %>% #We just want clarity guys
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  group_by(ID) %>%
 mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0)) %>% #we jump the days in which there is no data such as the day of the echoMRI
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
    date_rel = DATE - first(DATE),
    FI_cum =cumsum(corrected_intake_kcal)) %>% 
  filter(!ID %in% c(7875, 7874))  #we want to eliminate the ID were used for antibody titration

ggplot(FI_data, aes(date_rel, FI_cum)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA) +
  geom_smooth() +
  theme_minimal()
