#This script aim to explore changes in BW after DREADD injection in orexin-cre mice
#we aim to split the orexin cre males into three groups:
#n=2 Males orexin - cre with inhibitory DREADD - uncertainty
#n=2 Males orexin - cre with control DREADD - uncertainty
#n=2 Males orexin - cre with control DREADD - no uncertainty

#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

#Data import####
BW_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import

BW_data <- BW_data %>% 
  filter(COHORT %in% c(11, 17)) %>% 
  arrange(DATE) %>% 
  mutate(GROUP = case_when(
    ID %in% c(305, 306, 307, 314, 316, 319, 322, 519,524) ~ "CONTROL_NO_UNCERT",
    ID %in% c(297, 313, 318, 320,520,523) ~ "INH_DREADD_UNCERT",
    ID %in% c(321, 323, 325, 522,525) ~ "CONTROL_UNCERT",
    )) %>% 
  group_by(ID) %>% 
  mutate(
    surgery_date = min(
      DATE[COMMENTS %in% c(
        "BILATERAL_CONTROL_DREADD_INJECTION_SURGERY",
        "BILATERAL_INH_DREADD_INJECTION_SURGERY"
      )],
      na.rm = TRUE
    )
  ) %>% 
  filter(DATE >= surgery_date) %>% 
  drop_na(GROUP) %>% 
  filter(!(COHORT == 11 & DATE >= as.Date("2025-05-09"))) %>%
  mutate(day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))) 

#plot####
plot <- ggplot(BW_data, aes(x = day_rel, y = BW)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = ID), alpha = 0.7) +
  theme_minimal() +
  ylab("Body weight (g)") +
  xlab("Date") +
  facet_wrap(~ID)

plot

