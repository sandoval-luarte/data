#we aim to evaluate changes in 24h TEE in middle age C57 after different stages of feeding:
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
library(ggpubr)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#TEE####
sable_TEE_data <- sable_dwn %>% # Load the data
  filter(COHORT ==2) %>%   #we only want C57 mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1",
                     "SABLE_DAY_2",
                     "SABLE_DAY_3",
                     "SABLE_DAY_4") ~ "peak obesity",
    sable_idx %in% c("SABLE_DAY_5",
                     "SABLE_DAY_6",
                     "SABLE_DAY_7",
                     "SABLE_DAY_8") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_9",
                     "SABLE_DAY_10",
                     "SABLE_DAY_11",
                     "SABLE_DAY_12") ~ "BW maintenance", 
    sable_idx %in% c("SABLE_DAY_13",
                     "SABLE_DAY_14",
                     "SABLE_DAY_15",
                     "SABLE_DAY_16") ~ "BW regain"
  )) %>% 
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
  group_by(ID,complete_days,is_complete_day,SABLE) %>% 
  mutate(tee = sum(value)*(1/60)) %>% 
  filter( is_complete_day ==1, complete_days==2) %>% 
  ungroup() %>% 
  group_by(ID,SABLE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ))
  


sable_min_plot_tee <- sable_TEE_data %>% 
  mutate(SABLE = factor(SABLE, 
                        levels = c("baseline", "peak obesity", "BW loss", 
                                   "BW maintenance", "BW regain"))) %>% 
  filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) #cage 5 in at least one SABLE stage


#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        #   strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

plot <- sable_min_plot_tee %>%
  ggplot(aes(x = SABLE, y = tee, color = DRUG, group = ID)) +
  
  # individual trajectories
  geom_line(alpha = 0.3) +   
  geom_point(size = 2, alpha = 0.3) +  
  
  # add ID labels
  geom_text(aes(label = ID), 
            size = 3,          # font size of labels
            vjust = -0.5,      # vertical position
            alpha = 0.7,       # transparency so they don’t clutter too much
            show.legend = FALSE) +  
  
  # mean ± SD ribbon
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
  
  theme_minimal() +
  labs(y = "24 h TEE", color = "Drug", fill = "Drug") +
  facet_wrap(~GROUP) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot



