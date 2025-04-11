#This script aim to explore BW changes after virus injection in orexin-cre and wt mice

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr) #to use drop_na()
library(lme4)


#BW####
  
  BW_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import
  
  BW_data <- BW_data %>% 
    filter(COHORT == 11) %>% 
    filter(ID %in% c(307, 316, 318, 320)) %>% 
    group_by(ID) %>% 
     mutate(day_rel = BW - first(BW))
  
  n_distinct(BW_data$ID) #here we know there is 27 animals in total n=18 orexin cre and 9 c57
  
  plot <- BW_data %>% 
    ggplot(aes(DATE, day_rel, group = ID, color = SEX)) +
    geom_point() +
    geom_line() +
    facet_wrap(STRAIN~SEX)+
    geom_text(data = BW_data,
              aes(label = ID), 
              hjust = -0.2, vjust = 0.5, 
              size = 3, show.legend = FALSE)  #check BW 316

 