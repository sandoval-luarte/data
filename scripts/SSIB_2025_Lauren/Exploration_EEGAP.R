#BW anf fi exploration all cohort

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()


BW_data <- read_csv("../data/BW.csv") %>% 
  group_by(ID) %>% 
  mutate(bw_rel = BW - first(BW), date_rel= DATE - first(DATE)) %>% 
  ggplot(aes(date_rel,bw_rel)) +
  geom_point()+
  geom_text(aes(label = ID)) +
  facet_wrap(COHORT~STRAIN)

BW_data

