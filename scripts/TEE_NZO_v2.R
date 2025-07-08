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

#TEE####
sable_TEE_data <- sable_dwn %>% # Load the data
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
                     "SABLE_DAY_15") ~ "acute weight loss", 
    sable_idx %in% c("SABLE_DAY_16",
                     "SABLE_DAY_17",
                     "SABLE_DAY_18",
                     "SABLE_DAY_19") ~ "weight maintence")) %>% 
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
  filter(!ID %in% c(3715), is_complete_day ==1, complete_days==2) %>% #3715 died 
  ungroup() %>% 
    group_by(ID,SABLE) %>% 
  filter(tee>0, SABLE %in% c("baseline", "peak obesity"))

mdl_tee <- lmer(
    data = sable_TEE_data,
    tee ~ SABLE + (1|ID)
)
summary(mdl_tee)

emmeans(
    mdl_tee,
    pairwise ~ SABLE,
    type = "response" #here again how can I split by lights>?
)

# Create the plot tee 
sable_min_plot_tee <- sable_TEE_data %>% 
  mutate(SABLE = factor(SABLE, levels = c("baseline", "peak obesity","acute weight loss", "weight maintence"))) %>%
  ggplot(aes(SABLE, tee)) +
  stat_summary(
    fun = mean, 
    # geom = "col", 
    fill = "gray", 
    alpha = 0.5) +
  geom_point(aes(group = ID), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.5) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.7,
    shape = 21,
    color = "black",
    fill = "red"
  ) +
  labs(x = NULL, y = "24h TEE ") +
  theme_classic() 
sable_min_plot_tee


  plot <- sable_TEE_data %>% 
  ggplot(aes(complete_days,tee, group=ID))+
  geom_line()+
  facet_wrap(~SABLE)
plot

#delta calculation TEE

delta <- sable_TEE_data %>% 
  ungroup() %>% 
  group_by(ID,SABLE) %>% 
  summarise(tee = mean(tee)) %>% 
  summarise(delta_tee = tee[SABLE =="peak obesity"] - tee[SABLE =="baseline"]) 

write_csv(x = delta,file = "delta_tee.csv")
