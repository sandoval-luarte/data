#we aim to evaluate changes in 24h TEE in food restricted Female C57 and NZO after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#TEE####
# build the summarized dataset
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  filter(SEX == "F") %>% # Just females from both strains
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726)) %>% # ad lib NZO
  filter(!ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                    7882, 7883)) %>% # ad lib C57
  filter(!ID %in% c(3723, 3724, 3725)) %>% # CAGE 5 ISSUES NZO AND 3724 CAGE 6 ISSUES 
  filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
  filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "baseline",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3","SABLE_DAY_4") ~ "peak obesity",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7","SABLE_DAY_8") ~ "BW loss",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12") ~ "BW maintenance", 
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  # keep both complete days
  filter(is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop")

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )


ggplot(sable_TEE_data, aes(x = SABLE, y = tee,  group = ID)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 2, alpha = 0.5) +
  geom_text(aes(label = ID), size = 2.5, show.legend = FALSE) +
  facet_wrap(~ID) +
  labs(y = "TEE (kcal/day)") +
  theme_minimal()

