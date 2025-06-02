#we aim to evaluate the metabolic changes that occurs in NZO female before and after 12 weeks with LFD

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
# TEE####
sable_tee_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11"), "After", "Before")) %>% 
    filter(grepl("kcal_hr_*", parameter)) %>% # just to see TEE in kcal first
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
    group_by(ID,complete_days,SABLE,is_complete_day) %>% 
   mutate(tee = sum(value)*(1/60)) %>% 
  filter(!ID %in% c(3715,3723), is_complete_day ==1) %>% #3715 died and something weird happens with 3723 
#  filter(!ID %in% c(3707,3713,3706,3709,3716,3712,3717), is_complete_day ==1) %>% # those animals has weird data
  ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1)

mdl_tee <- lmer(
    data = sable_tee_data,
    tee ~ SABLE + (1|ID)
)
summary(mdl_tee)

emmeans(
    mdl_tee,
    pairwise ~ SABLE,
    type = "response" #here again how can I split by lights>?
)

# Create the plot tee 
# Add ID labels to the lowest TEE values
sable_min_plot_tee <- sable_tee_data %>%
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, tee)) +
  # Bar for mean value
  stat_summary(fun = mean, fill = "gray", alpha = 0.5) +
  # Individual points with transparency
  geom_point(aes(group = ID), alpha = 0.5) +
  # Lines connecting paired observations
  geom_line(aes(group = ID), alpha = 0.5) +
  # Mean with SEM as point and error bars
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.7,
    shape = 21,
    color = "black",
    fill = "red"
  ) +
  # Add text labels for IDs
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5, size = 3) + 
  # Axis labels
  labs(x = NULL, y = "24h TEE") +
  # White background
  theme_classic() #+
  #facet_grid(~lights) #pareciera q son los mismos datos, weird

sable_min_plot_tee
 #something weid happen with ID 3723
