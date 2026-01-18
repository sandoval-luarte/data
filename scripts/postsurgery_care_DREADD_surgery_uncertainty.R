#This script aim to explore changes in BW after DREADD injection in orexin-cre
#we aim split the orexin cre males into three groups:
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
  filter(COHORT ==17) %>% 
  mutate(GROUP = case_when(
    ID %in% c(519,524) ~ "CONTROL_NO_UNCERT",
    ID %in% c(520,523) ~ "INH_DREADD_UNCERT",
    ID %in% c(522,525) ~ "CONTROL_UNCERT",
    )) 

#plot####
plot <- ggplot(BW_data, aes(x = DATE, y = BW)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = ID), alpha = 0.7) +
  theme_minimal() +
  ylab("Body weight (g)") +
  xlab("Date") +
  facet_wrap(~ID)

plot

