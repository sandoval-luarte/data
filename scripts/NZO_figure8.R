#we aim to evaluate changes in 24h locomotion in middle age NZO after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#locomotion####
sable_locomotion_data <- sable_dwn %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1",
                     "SABLE_DAY_2",
                     "SABLE_DAY_3",
                     "SABLE_DAY_4",
                     "SABLE_DAY_5",
                     "SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "baseline",
    sable_idx %in% c("SABLE_DAY_8",
                     "SABLE_DAY_9",
                     "SABLE_DAY_10",
                     "SABLE_DAY_11") ~ "peak obesity",
    sable_idx %in% c("SABLE_DAY_12",
                     "SABLE_DAY_13",
                     "SABLE_DAY_14",
                     "SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16",
                     "SABLE_DAY_17",
                     "SABLE_DAY_18",
                     "SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20",
                     "SABLE_DAY_21",
                     "SABLE_DAY_22",
                     "SABLE_DAY_23") ~ "BW regain",
  )) %>% 
  filter(grepl("AllMeters_*", parameter)) %>% # just to see locomotion
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
  group_by(ID,complete_days,is_complete_day,SABLE,lights) %>% 
  filter(!ID %in% c(3715,3712), is_complete_day ==1, complete_days==2) %>% #3715 died and 3723 has issues with the sable ,3723,3718,3709,3725
  ungroup() %>% 
  group_by(ID,SABLE,lights) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )) %>% 
  mutate(value = max(value))



sable_min_plot_locomotion <- sable_locomotion_data %>% 
  mutate(SABLE = factor(SABLE, 
                        levels = c("baseline", "peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) %>% 
  filter(!(ID %in% c(3723,3725,3718,3709)))


#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        #   strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

plot <- sable_min_plot_locomotion %>%
  ggplot(aes(x = SABLE, y = tee, color = DRUG, group = ID)) +
  
  # individual trajectories
  geom_line(alpha = 0.3) +   
  geom_point(size = 2, alpha = 0.3) +  
  
  # mean ± SD ribbon (need to set 'fill' separately from 'color')
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1), 
    geom = "ribbon", aes(group = DRUG, fill = DRUG), 
    alpha = 0.2, color = NA
  ) +
  
  # mean solid line
  stat_summary(
    fun = mean, geom = "line", aes(group = DRUG, color = DRUG), 
    size = 1.2
  ) +
  
  # mean dashed line (optional, if you want to keep it too)
  # stat_summary(fun = mean, geom = "line", aes(group = DRUG, color = DRUG), 
  #              size = 1.2, linetype = "dashed") +
  
  theme_minimal() +
  labs(y = "locomotion ", color = "Drug", fill = "Drug") +
  facet_wrap(~GROUP) +
  # format.plot+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot


sable_locomotion_max <- sable_dwn %>%
  filter(COHORT %in% c(3, 4, 5)) %>% 
  mutate(
    lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
    SABLE = case_when(
      sable_idx %in% paste0("SABLE_DAY_", 1:7) ~ "baseline",
      sable_idx %in% paste0("SABLE_DAY_", 8:11) ~ "peak obesity",
      sable_idx %in% paste0("SABLE_DAY_", 12:15) ~ "BW loss",
      sable_idx %in% paste0("SABLE_DAY_", 16:19) ~ "BW maintenance",
      sable_idx %in% paste0("SABLE_DAY_", 20:23) ~ "BW regain"
    )
  ) %>% 
  filter(grepl("AllMeters_*", parameter)) %>%
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>% 
  ungroup() %>% 
  group_by(ID, complete_days, is_complete_day, SABLE, lights) %>% 
  filter(!ID %in% c(3715, 3712),
         is_complete_day == 1,
         complete_days == 2) %>% 
  ungroup() %>% 
  group_by(ID, SABLE, lights, hr) %>%   # 🔹 include hr here
  summarise(
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )



