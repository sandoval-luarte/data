# This script aims to explore changes in body weight, food intake and body comp
# in middle age NZO and C57 after different stages of feeding:
#1: Basal, 2: peak obesity, 3: Acute body, weight loss 4: BW maintenance
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()
library(ggpubr)


##body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
         ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726,7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
                   7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
         ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 46 animals,
                        # n=22 NZO, n=24 C57


plot <- BW_data %>%
  ggplot(aes(day_rel, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~STRAIN*SEX*GROUP) +
 #  geom_smooth()+
#  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "day_rel",
    y = "BW (grams)")

plot

#body comp analysis####
##echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT > 1 & COHORT < 6) %>% #Just NZO females
  arrange(Date) %>% 
  filter(!ID %in% c(3712, 3715)) %>% #died during study 
 mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726,7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
              7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted"))%>% 
  group_by(ID) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index,GROUP) %>% 
  mutate(day_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean)) %>% 
  left_join(BW_data %>% select(ID, DIET_FORMULA,SEX,STRAIN), by = "ID") 

NZO_injections <-echoMRI_data %>%
  filter(STRAIN =="NZO/HlLtJ") %>% 
  slice_tail(n = 1) %>% 
  mutate(INJECTION = case_when(
    ID %in% c(3714,
              3728,
              3725,
              3724,
              3727,
              3720,
              3707,
              3709,
              3711,
              3713,
              3706) ~ "Vehicle",
    ID %in% c(3710,
              3729,
              3708,
              3723,
              3721,
              3722,
              3726,
              3719,
              3718,
              3716,
              3717) ~ "rtioxa_47"))
  

###adiposity index####
plot_echo <- NZO_injections %>%
  ggplot(aes(INJECTION, adiposity_index, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
   geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label
  
plot_echo

###injection assignation NZO####
plot_inject <- echoMRI_data %>%
  ggplot(aes(day_rel, adiposity_index, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP*SEX) +
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label

plot_inject

