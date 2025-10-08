#TOS 2025 LL CS TASK 1
#AIM> Building empirical estimates for BW regain
#y predicts (BW regain period)
#y1 , y1= (BW regain - BW peak obesity)/ BW peak obesity
#y2 , y2= (adiposity index (AI) BW regain - AI peak obesity)/ AI peak obesity
#y3 , y3= ROC BW = (((BW regain - BW maintenance)/BW maintenance))/days 
#y4 , y3= ROC AI = (((AI regain - AI BW maintenance)/AI BW maintenance))/days 

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models
library(ggpubr)  
library(lmerTest)
library(emmeans)

#BW DATA####

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO F AND C57 (M and F)
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726)) %>% # ad lib NZO
  filter(!ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                    7882, 7883)) %>% # ad lib C57
  #filter(!ID %in% c(3723, 3724, 3725)) %>% # CAGE 5 ISSUES NZO AND 3724 CAGE 6 ISSUES 
  #filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
  #filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    day_rel = DATE - first(DATE)) %>% 
  mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
      STRAIN == "C57BL6/J" & DATE == as.Date("2025-02-26") ~ "peak obesity",
      STRAIN == "C57BL6/J" & day_rel == 209 ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & day_rel == 161 ~ "BW loss",
      STRAIN == "C57BL6/J" & day_rel == 300 ~ "BW maintenance", 
      STRAIN == "NZO/HlLtJ" & day_rel == 196 ~ "BW maintenance",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    ) ) %>% 
  drop_na(STATUS)

BW_data %>%
  group_by(STRAIN, SEX,STATUS) %>%
  summarise(
    n_ID = n_distinct(ID)
  )
 #20 animals, 12 NZO, 8 C57. 16 F in total and 4 males in total

#ADIPOSITY INDEX DATA####
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO F AND C57 (M & F)
  group_by(ID) %>%
  arrange(Date) %>%
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726)) %>% # ad lib NZO
  filter(!ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                    7882, 7883)) %>% # ad lib C57
  #filter(!ID %in% c(3723, 3724, 3725)) %>% # CAGE 5 ISSUES NZO AND 3724 CAGE 6 ISSUES
  #filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
 # filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index,STRAIN,DIET_FORMULA, SEX) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-03-07") ~ "peak obesity",
      STRAIN == "NZO/HlLtJ" & Date == as.Date("2025-02-20") ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX=="F" & Date == as.Date("2025-03-31") ~ "BW loss",
      STRAIN == "C57BL6/J" & ID==7872 & Date == as.Date("2025-04-21") ~ "BW loss",
      STRAIN == "C57BL6/J" & SEX=="M" & Date == as.Date("2025-04-21") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-06-05") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & Date == as.Date("2025-05-27") ~ "BW maintenance",
      STRAIN == "C57BL6/J" & Date %in% as.Date(c("2025-08-28", #7878 last day of echo
                                                 "2025-09-01", #7872 last day of echo
                                                 "2025-09-02", #7874 last day of echo
                                                 "2025-09-05",#7877 last day of echo
                                                 "2025-08-27", #7863 last day of echo
                                                 "2025-08-27" #7861, 7863 last day of echo
                                                 )
      ) ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & Date %in% as.Date(c("2025-07-08", #3708 last day of echo
                                                  "2025-07-09",#3710 last day of echo
                                                  "2025-07-14",#3714 last day of echo
                                                  "2025-07-17", #3720, 3721, 3722 last day of echo
                                                  "2025-07-22", #3727, 3728, 3729 last day of echo
                                                  "2025-07-21" #3723, 3724, 3725 last day of echo
                                                  
      ) 
      ) ~ "BW regain",
      
      TRUE ~ NA_character_
    )
  ) %>% 
  drop_na(STATUS)

echoMRI_data  %>%
  group_by(STRAIN, SEX,STATUS) %>%
  summarise(
    n_ID = n_distinct(ID)
  )

##y1####
y1_data <- BW_data %>%
  filter(STATUS %in% c("peak obesity", "BW regain")) %>%  # keep only needed stages
  select(ID, STRAIN, SEX, STATUS, BW) %>%
  pivot_wider(names_from = STATUS, values_from = BW) %>%
  mutate(
    y1 = ( `BW regain` - `peak obesity` ) / `peak obesity`
  )
y1_data

y1_data  %>%
  group_by(STRAIN, SEX) %>%
  summarise(
    n_ID = n_distinct(ID)
  )

write_csv(y1_data, "../data/y1_data.csv") # Save as CSV

## y2 ####
y2_data <- echoMRI_data %>%
  filter(STATUS %in% c("peak obesity", "BW regain")) %>%
  select(ID, STRAIN, SEX, STATUS, adiposity_index) %>%
  pivot_wider(names_from = STATUS, values_from = adiposity_index) %>%
  mutate(
    y2 = (`BW regain` - `peak obesity`) / `peak obesity`
  )
y2_data

y2_data  %>%
  group_by(STRAIN, SEX) %>%
  summarise(
    n_ID = n_distinct(ID)
  )


write_csv(y2_data, "../data/y2_data.csv") # Save as CSV


##y3####

y3_data <- BW_data %>%
  filter(STATUS %in% c("BW maintenance", "BW regain")) %>% 
  select(ID, STRAIN, SEX, STATUS, BW, DATE) %>%
  pivot_wider(names_from = STATUS, values_from = c(BW, DATE)) %>%
  mutate(
    days = as.numeric(`DATE_BW regain` - `DATE_BW maintenance`),
    y3 = ((`BW_BW regain` - `BW_BW maintenance`) / `BW_BW maintenance`) / days
  )
y3_data

y3_data  %>%
  group_by(STRAIN, SEX) %>%
  summarise(
    n_ID = n_distinct(ID)
  )

write_csv(y3_data, "../data/y3_data.csv") # Save as CSV

## y4 ####
y4_data <- echoMRI_data %>%
  filter(STATUS %in% c("BW maintenance", "BW regain")) %>% 
  select(ID, STRAIN, SEX, STATUS, adiposity_index, Date) %>%
  pivot_wider(names_from = STATUS, values_from = c(adiposity_index, Date)) %>%
  mutate(
    days = as.numeric(`Date_BW regain` - `Date_BW maintenance`),
    y4 = ((`adiposity_index_BW regain` - `adiposity_index_BW maintenance`) / 
            `adiposity_index_BW maintenance`) / days
  )
y4_data

y4_data  %>%
  group_by(STRAIN, SEX) %>%
  summarise(
    n_ID = n_distinct(ID)
  )

write_csv(y4_data, "../data/y4_data.csv") # Save as CSV




