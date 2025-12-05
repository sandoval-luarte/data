#We aim to evaluate changes in locomotion after 5 days of injections of RTIOXA 43 in females and males C57

#Libraries
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels
library(lme4)
library(stringr)
library(lme4)

zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#locomotion after RTIOXA 43 injections ----
##mean of day 1 and 2 light collapsed ----

sable_loc_data <- sable_dwn %>% 
  filter(COHORT == 6) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
         SEX = case_when(
           startsWith(as.character(ID), "1") ~ "M",
           startsWith(as.character(ID), "2") ~ "F",
           TRUE ~ NA_character_
         )) %>% 
  mutate(DRUG= case_when(
    ID==1001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==1001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==1004 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==1004 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1005 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1006 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1006 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==1006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==1007 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==1007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1007 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==2001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==2002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==2002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==2003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2003 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==2004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2004 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2004 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==2005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2005 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==2005 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2006 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==2006 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2007 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2007 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
  )) %>% 
  filter(grepl("AllMeters_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,SEX) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX,DRUG) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter(!ID %in% c(1003,2001), is_complete_day == 1, complete_days %in% c(1,2)) #ID 1003 DIED DURING EXPERIMENT AND ID 2005 WAS IN CAGE 5

sable_loc_data_minutes <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,SEX) %>%    # group per animal, drug, and day
  mutate(
    locomotion = value - lag(value),       # change in meters
    locomotion = if_else(locomotion < 0, 0, locomotion),  # remove negative jumps
    moving_min = if_else(locomotion > 0, 1, 0)             # 1 min per movement
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    # total minutes moved per day
    total_distance = sum(locomotion),      # total meters per day
    .groups = "drop"
  ) %>%
  group_by(ID, DRUG,SEX) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  )

sable_loc_data_minutes%>%
  group_by(DRUG,SEX) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

sable_loc_data_minutes<- sable_loc_data_minutes %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_medchem", "RTIOXA_43_donated")),
    SEX = factor(SEX)
  )


sable_loc_data_minutes_movers <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       
  group_by(ID, DRUG, complete_days, SEX) %>%    
  mutate(
    locomotion = value - lag(value),       
    locomotion = if_else(locomotion < 0, 0, locomotion),  
    moving_min = if_else(locomotion > 0, 1, 0)             
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    
    total_distance = sum(locomotion),      
    .groups = "drop"
  ) %>%
  group_by(ID, DRUG, SEX) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  
    avg_distance = mean(total_distance),
    .groups = "drop"
  ) %>%
  # add the mover column based on vehicle group
  group_by(ID, SEX) %>%
  mutate(
    mover = case_when(
      avg_moving_hr[DRUG == "vehicle"] > 10 ~ "high",
      avg_moving_hr[DRUG == "vehicle"] <= 10 ~ "low"
    )
  ) %>%
  ungroup()


sable_loc_data_minutes_movers<- sable_loc_data_minutes_movers %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_medchem", "RTIOXA_43_donated")),
    SEX = factor(SEX)
  )


## Plot locomotion in meters ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_distance , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ SEX) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_medchem" = "darkgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "total locomotion (meters in 24h)", x = "")

# Plot time spent moving in min ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ SEX) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_medchem" = "darkgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")

# Plot time spent moving in min MOVERS ----
ggplot(sable_loc_data_minutes_movers, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ mover) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_medchem" = "darkgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")


#I will leave just the data of RTIOXA-43 obtained from medchem

#locomotion after RTIOXA 43 injections ----
##mean of day 1 and 2 light collapsed ----

sable_loc_data <- sable_dwn %>% 
  filter(COHORT == 6) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
         SEX = case_when(
           startsWith(as.character(ID), "1") ~ "M",
           startsWith(as.character(ID), "2") ~ "F",
           TRUE ~ NA_character_
         )) %>% 
  mutate(DRUG= case_when(
    ID==1001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43",
    ID==1001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
    ID==1003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
    ID==1004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43",
    ID==1004 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
    ID==1006 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43",
    ID==1007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1007 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43",
    ID==2001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
    ID==2001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43",
    ID==2002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43",
    ID==2003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2004 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
    ID==2005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2005 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43",
    ID==2006 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43",
    ID==2007 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43",
  )) %>% 
  filter(grepl("AllMeters_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,SEX) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX,DRUG) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  filter(!ID %in% c(1003,2001), is_complete_day == 1, complete_days %in% c(1,2)) #ID 1003 DIED DURING EXPERIMENT AND ID 2005 WAS IN CAGE 5

sable_loc_data_minutes <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,SEX) %>%    # group per animal, drug, and day
  mutate(
    locomotion = value - lag(value),       # change in meters
    locomotion = if_else(locomotion < 0, 0, locomotion),  # remove negative jumps
    moving_min = if_else(locomotion > 0, 1, 0)             # 1 min per movement
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    # total minutes moved per day
    total_distance = sum(locomotion),      # total meters per day
    .groups = "drop"
  ) %>%
  group_by(ID, DRUG,SEX) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  )

sable_loc_data_minutes%>%
  group_by(DRUG,SEX) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

sable_loc_data_minutes<- sable_loc_data_minutes %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43")),
    SEX = factor(SEX)
  )


sable_loc_data_minutes_movers <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       
  group_by(ID, DRUG, complete_days, SEX) %>%    
  mutate(
    locomotion = value - lag(value),       
    locomotion = if_else(locomotion < 0, 0, locomotion),  
    moving_min = if_else(locomotion > 0, 1, 0)             
  ) %>%
  drop_na() %>% 
  summarise(
    total_moving_min = sum(moving_min),    
    total_distance = sum(locomotion),      
    .groups = "drop"
  ) %>%
  group_by(ID, DRUG, SEX) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  
    avg_distance = mean(total_distance),
    .groups = "drop"
  ) %>%
  # add the mover column based on vehicle group
  group_by(ID, SEX) %>%
  mutate(
    mover = case_when(
      avg_moving_hr[DRUG == "vehicle"] > 10 ~ "high",
      avg_moving_hr[DRUG == "vehicle"] <= 10 ~ "low"
    )
  ) %>%
  ungroup()


sable_loc_data_minutes_movers<- sable_loc_data_minutes_movers %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43")),
    SEX = factor(SEX)
  )


## Plot locomotion in meters ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_distance , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ SEX) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43" = "#66C2A5"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "total locomotion (meters in 24h)", x = "")

# Plot time spent moving in min ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ SEX) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43" = "#66C2A5"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")

# Plot time spent moving in min MOVERS ----
ggplot(sable_loc_data_minutes_movers, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~ mover) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43" = "#66C2A5"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")

#data analysis

lowmovers_data <- sable_loc_data_minutes_movers %>%
  filter(mover == "low")

model_low<- lmer(
  avg_moving_hr ~ DRUG + (1 | ID) ,
  data = lowmovers_data)

summary(model_low)
anova(model_low)
emmeans(model_low, pairwise ~ DRUG , adjust = "tukey")

# split by lights analysis ----
