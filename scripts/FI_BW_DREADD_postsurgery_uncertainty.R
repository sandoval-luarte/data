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
    filter(SEX == "M") %>% #We just want to follow the animals with surgery 
    filter(!ID %in% c(298, 315)) %>% #died in surgery
    filter(!ID %in% c(308, 317, 324)) %>% #C57BL6J not used in surgery
    group_by(ID) %>% 
    arrange(DATE) %>% 
     mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
            body_lag = lag(BW) - BW) %>% 
    group_by(ID) %>% 
    mutate(day_rel = DATE - first(DATE))
  
  n_distinct(BW_data$ID) #here we know there is 14 animals
  
# Subset rows where COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE"
highlight_data <- BW_data %>%
  filter(COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE")

plot <- BW_data %>% 
  ggplot(aes(DATE, body_lag, group = ID, color = SEX)) +
  geom_point() +
  geom_line() +
  facet_wrap(~STRAIN) +
  geom_text(data = BW_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE) +
  # Highlight using vertical lines
  geom_vline(data = highlight_data, 
             aes(xintercept = as.numeric(DATE)), 
             linetype = "dashed", color = "red", alpha = 0.7)
plot

#echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT == 11) %>% 
  filter(SEX == "M") %>% 
  filter(!ID %in% c(298, 315)) %>% #died in surgery
  filter(!ID %in% c(308, 317, 324)) %>% #C57BL6J not used in surgery
  select(ID,Fat,Lean,Weight,n_measurement,adiposity_index)

plot_echo <- echoMRI_data %>% 
  ggplot(aes(n_measurement, Lean, group = ID)) +
  geom_line() +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo


 