#This script aims to evaluate changes in BW, FI and body comp of animals after 5 days of IP injections of RTIOXA 43 or RTIOXA 47
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()

#body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(9, 12)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(8079, 8076, 8074, 8100, 8103, 8095) ~ "Control",
      ID %in% c(8078, 8077, 8075) ~ "RTI_47_10mgkg",
      ID %in% c(8099, 8096, 8102) ~ "RTI_43_10mgkg_Y",
      ID %in% c(8098, 8101, 8097) ~ "RTI_43_10mgkg_MED"
    )
  ) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals


plot <- BW_data %>% 
  ggplot(aes(day_rel, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP, scales = "free_x") +  # free x-axis scales
  geom_smooth() +
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + 
  labs(
    x = "Date",
    y = "BW (grams)"
  )

plot

