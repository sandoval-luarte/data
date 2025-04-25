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
  ) 
ggplot(BW_data, aes(date_rel, bw_rel)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_grid(~STRAIN)

