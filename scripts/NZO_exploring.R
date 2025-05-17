# This script aims to explore BW and FI in NZO female mice

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()

#bw####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 )) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = lag(BW) - BW) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals

# Subset rows where COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE"
highlight_data <- BW_data %>%
  filter(COMMENTS == "RESTRICTED_DAY_1")

plot <- BW_data %>% 
  ggplot(aes(DATE, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ID) +
  geom_smooth()+
  # Highlight using vertical lines
  geom_vline(data = highlight_data, 
    aes(xintercept = as.numeric(DATE)), 
   linetype = "dashed", color = "red", alpha = 0.7)
plot
