#we aim to evaluate the metabolic changes that occurs in NZO female before and after 10 weeks with LFD

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)

#functions####
zt_time <- function(hr){
    return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#SPA####
sable_meters_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
    filter(grepl("AllMeters_*", parameter)) %>% # just to see TEE in kcal first
    ungroup() %>% 
    group_by(ID, SABLE) %>% 
    mutate(
        value = cumsum(if_else(replace_na(value - lag(value),0)<0, 0, replace_na(value - lag(value),0))),
        zt_time = zt_time(hr),
        is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
        complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
    ) %>% 
    ungroup() %>% 
    group_by(ID, complete_days) %>% 
    mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
    ungroup() %>% 
    group_by(ID,complete_days,SABLE,is_complete_day) %>% 
    summarise(meters = abs(max(value) - min(value))) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1) 

sable_meters_data %>% 
    ggplot(aes(
        SABLE, meters
    )) +
    geom_point() +
    geom_line(aes(group=ID)) 

mdl_meters <- lmer(
    data = sable_meters_data,
    meters~SABLE+(1|ID)
)
summary(mdl_meters)

#FOOD INTAKE####

sable_min_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
    filter(grepl("FoodA_*", parameter)) %>% # just to see TEE in kcal first
    ungroup() %>% 
    group_by(ID, SABLE) %>% 
    mutate(
        value = cumsum(if_else(replace_na(value - lag(value),0)>0, 0, replace_na(value - lag(value),0))),
        zt_time = zt_time(hr),
        is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
        complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
    ) %>% 
    ungroup() %>% 
    group_by(ID, complete_days) %>% 
    mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
    ungroup() %>% 
    group_by(ID,complete_days,SABLE,is_complete_day) %>% 
    summarise(intake_gr = abs(max(value) - min(value))) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1) %>% 
    mutate(
        kcal = if_else(SABLE=="Before", intake_gr*3.1, intake_gr*3.82)
    )

mdl_intake <- lmer(
    data = sable_min_data%>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    kcal ~ SABLE + (1|ID)
)
summary(mdl_intake)

emmeans(
    mdl_intake,
    pairwise ~ SABLE,
    type = "response"
)

# Create the plot FI kcal
sable_min_plot <- sable_min_data %>%
    filter(!ID %in% c(3727, 3718, 3711, 3723)) %>% #issues with cage 5
    ggplot(aes(x = SABLE, y = kcal, color = as.factor(ID))) +  # or just value to see the hourly value
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
    )
sable_min_plot 

# TEE####
sable_tee_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
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
    summarise(tee = sum(value)*(1/60)) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1)

mdl_intake <- lm(
    data = sable_tee_data %>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    tee ~ SABLE
)
summary(mdl_intake)

emmeans(
    mdl_intake,
    pairwise ~ SABLE,
    type = "response"
)

# Create the plot tee 
sable_min_plot_tee <- sable_tee_data %>%
    filter(!ID %in% c(3727, 3718, 3711, 3723)) %>% #issues with cage 5
    ggplot(aes(x = SABLE, y = tee, color = as.factor(ID))) +  # or just value to see the hourly value
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
    )
sable_min_plot_tee 

#Change in bw before and after LFD in NZO mice ####

bw <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study
  group_by(ID) %>% 
  mutate(daily_delta=replace_na( (BW - lag(BW))/lag(BW), 0 ),
         cd = cumsum(daily_delta))

bw %>% 
ggplot(aes(x = DATE, y = daily_delta)) + 
  geom_smooth()+
  geom_point() +
  geom_line(aes(group = ID))+
  labs(
    x = "Date",
    y = "Body Weight (BW)"
  ) +
  theme_minimal()



# RQ####
sable_RQ_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
    filter(grepl("RQ_*", parameter)) %>% # just to see TEE in kcal first
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
    summarise(RQ = mean(value)) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1)

mdl_RQ <- lm(
    data = sable_RQ_data %>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    RQ ~ SABLE
)
summary(mdl_intake)

emmeans(
    mdl_RQ,
    pairwise ~ SABLE,
    type = "response"
)

# Create the plot RQ
sable_min_plot_RQ <- sable_RQ_data %>%
    filter(!ID %in% c(3727, 3718, 3711, 3723)) %>% #issues with cage 5
    ggplot(aes(x = SABLE, y = RQ, color = as.factor(ID))) +  # or just value to see the hourly value
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
    )
sable_min_plot_RQ 

# speedmeters####
sable_PedSpeed_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
    filter(grepl("PedSpeed_*", parameter)) %>% # just to see TEE in kcal first
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
    summarise(PedSpeed = mean(value)) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1)

mdl_speed <- lm(
    data = sable_PedSpeed_data %>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    PedSpeed ~ SABLE
)
summary(mdl_speed)

emmeans(
    mdl_speed,
    pairwise ~ SABLE,
    type = "response"
)

# Create the plot PedSpeed
sable_min_plot_PedSpeed <- sable_PedSpeed_data %>%
    filter(!ID %in% c(3727, 3718, 3711, 3723)) %>% #issues with cage 5
    ggplot(aes(x = SABLE, y = PedSpeed, color = as.factor(ID))) +  # or just value to see the hourly value
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
    )
sable_min_plot_PedSpeed


#FI intake over days####
FI_data_ <- read_csv("../data/FI.csv") %>% 
    filter(COHORT %in% c(3, 4, 5))%>%
    drop_na() %>% 
    group_by(ID) 
FI_data_

fi_plot <- FI_data_ %>%
    filter(corrected_intake_gr<15, DIET=="LFD") %>% 
    ggplot(aes(x = as.factor(DATE), y = corrected_intake_gr)) +  
    geom_point() +
    geom_smooth(aes(group=1))
fi_plot



#echoMRI analysis####
echomri_data <- read_csv("../data/echomri.csv") %>% 
    filter(COHORT %in% c(3, 4, 5)) %>%
    filter(ID != 3715) %>%  # Exclude animals that died during the study
  #  filter(Date %in% c("2024-11-12", "2024-11-19", "2024-11-26", "2025-02-20")) %>%   # Keep only selected dates
  select(Date,ID, adiposity_index, Fat, Lean) %>% 
   mutate(STATUS  = if_else(Date %in% c("2024-11-12", "2024-11-19", "2024-11-26"), "before", "after")) 
 
echomri_data_<-echomri_data %>%
    ggplot(aes(x = Date, y = Fat, color = as.factor(ID))) + #Fat
    geom_point(size = 3, alpha = 0.8) 
echomri_data_

mdl_echo <- lm(
    data = echomri_data %>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    Fat ~ STATUS
)
summary(mdl_echo)

emmeans(
    mdl_echo,
    pairwise ~ STATUS,
    type = "response"
)

#NZO mice, better known as butter mice

echomri_data_<-echomri_data %>%
    ggplot(aes(x = Date, y = adiposity_index, color = as.factor(ID))) + #Fat
    geom_point(size = 3, alpha = 0.8) 
echomri_data_

mdl_echo <- lm(
    data = echomri_data %>% filter(!ID %in% c(3727, 3718, 3711, 3723)),
    Lean ~ STATUS
)
summary(mdl_echo)

emmeans(
    mdl_echo,
    pairwise ~ STATUS,
    type = "response"
)

#Data analysis of for the food restriction####

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


