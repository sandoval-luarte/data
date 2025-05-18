# This script aims to explore BW, FI and locomotion in NZO female mice after different stages of feeding:
#1. Basal, 2: peak obesity, 3: Acute ody weight loss 4: BW maintenance

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
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
         ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
         ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals

# Subset rows where COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE"
highlight_data <- BW_data %>%
  filter(COMMENTS == "RESTRICTED_DAY_1")

plot <- BW_data %>% 

  ggplot(aes(day_rel, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
   geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) # Add ID labels+  
#   # Highlight using vertical lines
#   geom_vline(data = highlight_data, 
#     aes(xintercept = as.numeric(DATE)), 
#    linetype = "dashed", color = "red", alpha = 0.7)
 plot

#fi####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 )) %>% #died during study
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>%
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         date_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) 

n_distinct(BW_data$ID) #here we know there is 22 animals

# Subset rows where COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE"
highlight_data <- FI_cum %>%
  filter(COMMENTS == "RESTRICTED_DAY_1")

plot <- FI_data %>% 
 filter(corrected_intake_kcal <= 40) %>% 
  ggplot(aes(DATE,corrected_intake_kcal , group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_smooth()+
  # Highlight using vertical lines
  geom_vline(data = highlight_data, 
             aes(xintercept = as.numeric(DATE)), 
             linetype = "dashed", color = "red", alpha = 0.7)
plot

#sable####
sable_data <- readRDS("../data/sable_downsampled_data.rds") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 ))   #died during study

#locomotion####

sable_data_locomotion <- sable_data %>% 
  filter(parameter =="AllMeters")

