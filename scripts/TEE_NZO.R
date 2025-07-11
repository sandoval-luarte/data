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

#SPA####
sable_meters_data <- sable_dwn %>% # Load the data
    filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11"), "After", "Before")) %>% #CHECK THIS LINE
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
    mutate(meters = abs(max(value) - min(value))) %>% 
  filter(ID != 3715, is_complete_day ==1) %>% 
  ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1) 

sable_meters_data %>% 
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, meters)) +
  # Bar for mean value
  stat_summary(
    fun = mean, 
   # geom = "col", 
    fill = "gray", 
    alpha = 0.5  # Transparency for the bar
  ) +
  
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
  # Axis labels
  labs(x = NULL, y = "24 h SPA (meters)") +
  
  # White background
  theme_classic() +
  facet_grid(~lights)


mdl_meters <- lmer(
    data = sable_meters_data,
    meters~SABLE+(1|ID) 
)
summary(mdl_meters) #how lights play here?

emmeans(
  mdl_meters,
  pairwise ~ SABLE,
  type = "response"
)

#FOOD INTAKE####

sable_min_data <- sable_dwn %>% 
    filter(COHORT %in% c(3, 4, 5)) %>%  
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11"), "After", "Before")) %>% 
    filter(grepl("FoodA_*", parameter)) %>% 
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
    mutate(intake_gr = abs(max(value) - min(value))) %>% 
  filter(ID != 3715, is_complete_day ==1) %>% 
  ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1) %>% 
    mutate(
        kcal = if_else(SABLE=="Before", intake_gr*3.1, intake_gr*3.82)
    )

mdl_intake <- lmer(
    data = sable_min_data,
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
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, kcal)) +
  
  # Bar for mean value
  stat_summary(
    fun = mean, 
    # geom = "col", 
    fill = "gray", 
    alpha = 0.5  # Transparency for the bar
  ) +
  
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
  # Axis labels
  labs(x = NULL, y = "24h kcal food intake ") +
  
  # White background
  theme_classic() +
  facet_grid(~lights) #pareciera q son los mismos datos, weird

sable_min_plot

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
  filter(!ID %in% c(3715), is_complete_day ==1) %>% #3715 died and the rest of the animals were measured in cage 5 
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
sable_min_plot_tee <- sable_tee_data %>% 
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, tee)) +
  
  # Bar for mean value
  stat_summary(
    fun = mean, 
    # geom = "col", 
    fill = "gray", 
    alpha = 0.5  # Transparency for the bar
  ) +
  
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
  # Axis labels
  labs(x = NULL, y = "24h TEE ") +
  
  # White background
  theme_classic() +
  facet_grid(~lights) #pareciera q son los mismos datos, weird

sable_min_plot_tee

# bw####
sable_bw_data <- sable_dwn %>% # Load the data
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10"), "After", "Before")) %>% 
  filter(grepl("BodyMass_*", parameter)) %>% 
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
  summarise(bw = mean(value)) %>% 
  filter(!ID %in% c(3715,3727, 3718, 3711, 3723), is_complete_day ==1) %>% #3715 died and the rest of the animals were measured in cage 5 
 # filter(!ID %in% c(3706,3707,3708, 3714, 3721, 3728)) %>% #issues with body mass sensor 
  ungroup() %>% 
  group_by(SABLE, ID) %>% 
  slice_max(order_by = complete_days, n=1)

mdl_bw <- lmer(
  data = sable_bw_data,
  bw ~ SABLE + (1|ID)
)
summary(mdl_bw)

emmeans(
  mdl_bw,
  pairwise ~ SABLE,
  type = "response"
)

# Create the plot BW
sable_min_plot_bw  <- sable_bw_data %>% 
mutate(SABLE = factor(SABLE, levels = c("Before", "After"))) %>%
  ggplot(aes(SABLE, bw)) +
  
  # Bar for mean value
  stat_summary(
    fun = mean, 
    # geom = "col", 
    fill = "gray", 
    alpha = 0.5  # Transparency for the bar
  ) +
  
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
  # Axis labels
  labs(x = NULL, y = "BW (g) ") +
  
  # White background
   theme_classic()  
sable_min_plot_bw 
  # # Individual points with transparency
  # geom_jitter(aes(group = ID), width = 0.1, alpha = 0.5) +
  # # Add ID labels near the points
  # geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) 


# speedmeters####
sable_PedSpeed_data <- sable_dwn %>% 
    filter(COHORT %in% c(3, 4, 5)) %>%   
    mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
    mutate(SABLE  = if_else(sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10", "SABLE_DAY_11"), "After", "Before")) %>% 
    filter(grepl("PedSpeed_*", parameter)) %>% 
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
   mutate(PedSpeed = mean(value)) %>% 
    filter(ID != 3715, is_complete_day ==1) %>% 
    ungroup() %>% 
    group_by(SABLE, ID) %>% 
    slice_max(order_by = complete_days, n=1)

mdl_speed <- lmer(
    data = sable_PedSpeed_data %>% 
    PedSpeed ~ SABLE + (1|ID)
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

#echoMRI analysis####
#Our aim is to take the slope between the first 5 measurements and the 4 last measurements

echomri_data <- read_csv("../data/echomri.csv") %>% 
    filter(COHORT %in% c(3, 4, 5)) %>%
    filter(ID != 3715) %>% 
  select(Date,ID, adiposity_index, Fat, Lean) %>% 
   mutate(STATUS  = if_else(Date < "2025-01-01", "period_1", "period_2")) 
 
#####fat####
FAT<-echomri_data %>%
    ggplot(aes(x = Date, y = Fat)) + #Fat
  geom_point(aes(group = ID), alpha = 0.5) +
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  geom_line(aes(group = ID), alpha = 0.5) +
  theme_classic() 
FAT #3711 3718 AND 3727 DID NOT RECOVER FROM FASTING

#Alternative model
model_data <- echomri_data %>% 
    group_by(STATUS, ID) %>% 
    mutate(
        time = as.numeric(Date - min(Date))
    )
model <-  lmer(Fat ~ STATUS * time + (1 | ID), data = model_data)
summary(model)

trend <- emtrends(
    model,
    pairwise ~ STATUS,
    var = "time"
)
trend

#NZO mice, better known as butter mice
###Lean####
LEAN<-echomri_data %>%
  ggplot(aes(x = Date, y = Lean)) + #Lean
  geom_point(aes(group = ID), alpha = 0.5) +
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  geom_line(aes(group = ID), alpha = 0.5) +
  theme_classic() 
LEAN

lean_data <- echomri_data %>% 
    group_by(ID, STATUS) %>% 
    mutate(
        time = as.numeric(Date - min(Date))
    )

mdl_echo_lean <- lmer(
  data = lean_data,
  Lean ~ STATUS * time + (1|ID)
)
summary(mdl_echo_lean)

emtrends(
  mdl_echo_lean,
  pairwise ~ STATUS,
  var = "time"
)

#####Adiposity index####
AI<-echomri_data %>%
    ggplot(aes(x = Date, y = adiposity_index)) + #Lean
    geom_point(aes(group = ID), alpha = 0.5) +
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
    geom_line(aes(group = ID), alpha = 0.5) +
    theme_classic() 
AI

AI_data <- echomri_data %>% 
    group_by(ID, STATUS) %>% 
    mutate(
        time = as.numeric(Date - min(Date))
    )

mdl_echo_AI <- lmer(
    data = AI_data,
    adiposity_index ~ STATUS * time + (1|ID)
)
summary(mdl_echo_AI)

emtrends(
    mdl_echo_AI,
    pairwise ~ STATUS,
    var = "time"
)

####Food intake period 1 and 2####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%
  filter(ID != 3715) %>% #DIED
    select(DATE,ID, corrected_intake_kcal, DIET) %>% 
  mutate(STATUS  = if_else(DATE < "2025-01-01", "period_1", "period_2")) %>% 
  drop_na()

#####FI over the period 1 and 2####
FI_data_ <-FI_data %>%
  filter(corrected_intake_kcal<50) %>% 
  filter(DIET=="LFD") %>% 
  ggplot(aes(x = DATE, y = corrected_intake_kcal)) + 
  geom_point(aes(group = ID), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.5) +
  theme_classic() 
FI_data_

#Alternative model for FI
model_data <- FI_data %>% 
  group_by(STATUS, ID) %>% 
  mutate(
    time = as.numeric(DATE - min(DATE))
  )

model <-  lmer(corrected_intake_kcal ~ STATUS * time + (1 | ID), data = model_data)
summary(model)

trend <- emtrends(
  model,
  pairwise ~ STATUS,
  var = "time"
)
trend


#Data analysis of for the food restriction####

FI_data <- read_csv("../data/FI.csv") %>% 
    filter(COHORT %in% c(3, 4, 5)) %>% 
    filter((ID %in% c(3710,3714,3720,3721,3722,3723,3724,3725,3727,3728,3729))) %>% #just for restricted animals
    # filter(corrected_intake_gr < 15) %>%  #TYPING MISTAKES
   # filter(DATE > "2025-02-01") %>% 
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


