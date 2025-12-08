#We aim to evaluate changes in locomotion after injections of RTIOXA 47 in females and males C57 and NZO

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
library(car)  # car for Anova(), vif()
library(effsize) # for Cohen's d; install if needed

zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#locomotion after RTIOXA 47 injections ----
##mean of day 1 and 2 light collapsed ----

sable_loc_data <- sable_dwn %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
    SABLE = case_when(
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
      STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!ID %in% c(3715,3712)) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728, 7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729, 7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"),
      HEALTH = case_when(
        ID %in% c(3708, 3720, 3721, 3722, 3725, 3728, 3729,7861, 7863, 7872, 7877, 7878) ~ "non responders to food restriction",
        ID %in% c(3714, 3710, 3723, 3724, 3727, 7865, 7866, 7874) ~ "responders to food restriction",
        TRUE ~ "ad lib (not restricted)"
      )
  ) %>%
  group_by(ID, DRUG,SEX,STRAIN,GROUP,DIET_FORMULA,SABLE,HEALTH) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX,DRUG,STRAIN,DIET_FORMULA,SABLE,HEALTH) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
 filter( is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("baseline", "peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))) %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki")) %>% 
   filter(!ID %in% c(3718, 3719)) #potential diabetic mice
  
 # filter(!ID %in% c(3708, 3720, 3721, 3722, 3725, 3728, 3729,7861, 7863, 7872, 7877, 7878)) #NON RESPONDERS TO FOOD RESTRICTION
  

sable_loc_data_minutes <- sable_loc_data %>%
  filter(SABLE =="BW regain") %>% 
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,SEX,STRAIN,GROUP,DIET_FORMULA,SABLE,HEALTH) %>%    # group per animal, drug, and day
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
  group_by(ID, DRUG,SEX,STRAIN,GROUP,DIET_FORMULA,SABLE,HEALTH) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  ) 

sable_loc_data_minutes%>%
  group_by(STRAIN,SABLE,HEALTH) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

sable_loc_data_minutes<- sable_loc_data_minutes %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    SEX = factor(SEX),
    GROUP = factor(GROUP),
    DIET_FORMULA = factor(DIET_FORMULA),
    STRAIN = factor(STRAIN),
    SABLE = factor(SABLE),
    HEALTH = factor(HEALTH)
  )


## Plot locomotion in meters ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_distance , fill = DRUG)) +
 # geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~HEALTH) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "locomotion (meters in 24h)", x = "")

# data analysis----
#to evaluate if vehicle is different to RTIOXA 47 group in terms of meters/24h

# Quick sample sizes per cell
sable_loc_data_minutes %>%
  group_by(DRUG,HEALTH) %>%
  summarise(n = n(), .groups = "drop") %>%
  print(n = Inf)

# Base model: DRUG effect controlling for GROUP, SEX, STRAIN
model_base <- lm(
  avg_distance ~ DRUG + HEALTH,
  data = sable_loc_data_minutes
)

# Model with interaction between DRUG and GROUP (test whether drug effect depends on feeding group)
model_inter <- lm(
  avg_distance ~ DRUG + HEALTH,
  data = sable_loc_data_minutes
)

emm_group <- emmeans(model_inter, ~ DRUG | HEALTH)
pairs(emm_group)

model_subgroup <- lm(
  avg_distance ~ DRUG + HEALTH,
  data = sable_loc_data_minutes
)
emm_sub <- emmeans(model_subgroup, ~ DRUG | HEALTH)
pairs(emm_sub) #NO EFFECTS IN TOTAL METERS


# Plot time spent moving in min ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
#  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~HEALTH) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)

#data analysis to evaluate if ID 3708 is outlier

subset_data <- sable_loc_data_minutes %>%
  filter(GROUP == "restricted")


ggplot(subset_data, aes(y = avg_moving_hr, x = "")) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
  geom_point(aes(label = ID), position = position_jitter(width = 0.1)) +
  geom_text_repel(aes(label = ID)) +
  labs(y = "Avg Moving Hours", x = "") +
  theme_minimal() #IT SEEMS LIKE ID 3708 IS AN OUTLIER

Q1 <- quantile(subset_data$avg_moving_hr, 0.25)
Q3 <- quantile(subset_data$avg_moving_hr, 0.75)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

subset_data %>%
  mutate(outlier_IQR = avg_moving_hr < lower_bound | avg_moving_hr > upper_bound) 

#conclusion 3708 IS an outlier so lets run the analysis without those ID

sable_loc_data_minutes<- sable_loc_data_minutes %>% 
  filter(!ID == 3708) 
  
# Plot time spent moving in min without ID 3708 ----
ggplot(sable_loc_data_minutes, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
  #  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(~GROUP) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "", title = 'time spent moving in food restricted animals')+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)

#data analysis

# Base model: DRUG effect controlling for GROUP, SEX, STRAIN
model_base <- lm(
  avg_moving_hr  ~ DRUG + GROUP + SEX + STRAIN,
  data = sable_loc_data_minutes
)

# Model with interaction between DRUG and GROUP (test whether drug effect depends on feeding group)
model_inter <- lm(
  avg_moving_hr  ~ DRUG * GROUP + SEX + STRAIN,
  data = sable_loc_data_minutes
)

emm_group <- emmeans(model_inter, ~ DRUG | GROUP)
pairs(emm_group) 

#so RTIOXA-47 decreased time spent moving in ad lib animals 
#we should eliminate 3708 considering that ID as an outlier
#these effects are only significant when we collapsed SEX and STRAIN

emm_sex <- emmeans(model_inter, ~ DRUG | SEX)
pairs(emm_sex)

emm_strain <- emmeans(model_inter, ~ DRUG | STRAIN)
pairs(emm_strain)

restricted_data <- sable_loc_data_minutes %>%
  filter(GROUP == "restricted")

model_restricted <- lm(
  avg_moving_hr ~ DRUG * SEX * STRAIN,
  data = restricted_data
)

emm_restricted <- emmeans(model_restricted, ~ DRUG | SEX * STRAIN)
pairs(emm_restricted) 
#conclusion: there is no effect of RTIOXA 47 in time spent moving 

# mean of day 1 and 2 separated by lights----

sable_loc_data_lights <- sable_dwn %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
         SABLE = case_when(
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
           STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain",
           TRUE ~ NA_character_)) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728, 7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729, 7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )) %>%
  filter(!ID %in% c(3715,3712,3720,3721)) %>% #just testing what happens if we eliminate 20 and 21
  filter(grepl("AllMeters_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,SEX,lights,GROUP,STRAIN,DIET_FORMULA) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX,DRUG,lights,GROUP,STRAIN,DIET_FORMULA) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  filter(!complete_days %in% c(0, 3)) %>% 
  mutate(
    SABLE = factor(SABLE,
                   levels = c("baseline", "peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))) %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki"))

sable_loc_data_minutes_lights <- sable_loc_data %>%
  arrange(ID, complete_days, hr) %>%       # make sure data is ordered
  group_by(ID, DRUG, complete_days,SEX,lights,GROUP,STRAIN,DIET_FORMULA) %>%    # group per animal, drug, and day
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
  group_by(ID, DRUG,SEX,lights,GROUP,STRAIN,DIET_FORMULA) %>%
  summarise(
    avg_moving_hr = mean(total_moving_min/60),  # average across days
    avg_distance = mean(total_distance),
    .groups = "drop"
  )

sable_loc_data_minutes_lights%>%
  group_by(SEX,lights,STRAIN) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

sable_loc_data_minutes_lights<- sable_loc_data_minutes_lights %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    SEX = factor(SEX),
    GROUP = factor(GROUP),
    DIET_FORMULA = factor(DIET_FORMULA),
    STRAIN = factor(STRAIN)
  )

# Plot total locomotion in meters LIGHTS ON AND OFF ----
ggplot(sable_loc_data_minutes_lights, aes(x = DRUG, y =  avg_distance , fill = DRUG)) +
#  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(lights~ SEX*GROUP*STRAIN) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "total locomotion (meters in 24h)", x = "")

## Plot time spent moving in min LIGHTS ON AND OFF ----
ggplot(sable_loc_data_minutes_lights, aes(x = DRUG, y =  avg_moving_hr , fill = DRUG)) +
#  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) +   # connect the same ID across drugs
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +   # bars with mean ± SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,   # individual points
             position = position_jitter(width = 0.1)) +
  facet_grid(lights~ SEX*GROUP*STRAIN) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_47" = "orange"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "time spent moving (hr in 24h)", x = "")+
  geom_text_repel(aes(label = ID),
                  size = 3, alpha = 0.7)





#SPA after RTIOXA 47 injections ----

