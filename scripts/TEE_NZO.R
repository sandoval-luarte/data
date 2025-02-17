#we aim to evaluate the basal SPA of NZO female prior to LFD and after 10 diets with LFD

#Libraries
library(dplyr) #to open a RDS and use pipe

#Data analysis of for the food restriction

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter((ID %in% c(3710,3714,3720,3721,3722,3723,3724,3725,3727,3728,3729))) %>% #just for restricted animals
 # filter(corrected_intake_gr < 15) %>%  #TYPING MISTAKES
  filter(DATE > "2025-02-01") %>% 
  drop_na()

 fi_plot_gr <- FI_data %>% 
  ggplot(aes(x = DATE, y = corrected_intake_gr, group =as.factor(ID))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(color = as.factor(ID))) +  # Add color for better distinction
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  labs(
    title = "Corrected Intake Over Days",
    x = "Date",
    y = "Corrected Intake (gr.)"
  ) +
  theme_minimal() +
  facet_wrap(~ID)
fi_plot_gr    

#TEE analysis before and after LFD in NZO mice

sable_min_data <- readRDS(file = "../data/sable_downsampled_data.rds") %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% # just to see TEE in kcal first
group_by(ID,sable_idx,lights,date) %>% 
  summarise(sum_tee = mean(abs(value))) %>% 
  mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
  ungroup() %>% 
  group_by(ID,SABLE,lights) %>% 
  summarise(mean_TEE=mean(sum_tee))

# Create the plot
sable_min_plot <- sable_min_data %>%
  ggplot(aes(x = SABLE, y = mean_TEE, color = as.factor(ID))) +  # or just value to see the hourly value
  geom_line(aes(group = ID)) +              # Connect lines across days for each ID
  geom_point(size = 3, alpha = 0.8) +                 # Add individual points
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +
  facet_wrap(~lights) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.5,
    shape = 21,
    color = "black",
    fill = "red",
    position = position_dodge(width = 0.2),  # Mean and SEM as big points,
    aes(group = SABLE)
  ) 
sable_min_plot 

#We conclude cage 5 could have calibration issues


#All meters analysis before and after LFD in NZO mice

sable_min_data <- readRDS(file = "../data/sable_downsampled_data.rds") %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  filter(grepl("AllMeters_*", parameter)) %>% # just to see TEE in kcal first
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(
    value = cumsum(if_else(replace_na(value - lag(value),0)<0, 0, replace_na(value - lag(value),0))),
    zt_time = as.numeric(DateTime - min(DateTime)),
    complete_days = cumsum(if_else(zt_time %% 86400 == 0, 1, 0))
  ) %>% 
  ungroup() %>% 
  group_by(ID,sable_idx,complete_days,lights) %>% 
  summarise(max_allmeters = max(value) - min(value)) %>% 
  mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
  ungroup() %>% 
  group_by(ID,SABLE,lights) %>% 
  summarise(mean_allmeters= mean(max_allmeters)) %>% 
  filter(ID != 3715)

sable_min_data %>% 
  filter(ID == 3706) %>% 
  ggplot(aes(
    DateTime, value
  )) + 
  geom_point()

# Create the plot
sable_min_plot <- sable_min_data %>%
  ggplot(aes(x = SABLE, y = mean_allmeters, color = as.factor(ID))) +  # or just value to see the hourly value
  geom_line(aes(group = ID)) +              # Connect lines across days for each ID
  geom_point(size = 3, alpha = 0.8) +                 # Add individual points
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.5,
    shape = 21,
    color = "black",
    fill = "red",
    position = position_dodge(width = 0.2),  # Mean and SEM as big points,
    aes(group = SABLE)
  ) +
  facet_wrap(~lights)
sable_min_plot 

#We conclude cage 5 could have calibration issues


# mdl_ee <- lme4::lmer(
#   data = sable_hr_data,
#   (mean_hr_tee) ~ lights + (1|ID)
# )
# summary(mdl_ee)
# 
# emmeans::emmeans(
#   mdl_ee,
#   pairwise ~ lights,
#   type = "response"
# )
# 
# sable_hr_plot

