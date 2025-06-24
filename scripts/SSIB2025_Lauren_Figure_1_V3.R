# We aim to graph changes in body weight, body comp, FI and TEE in NZO mice after 12 weeks of LFD
#Figure 1 poster SSIB Lauren

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models
library(ggpubr)  
library(lmerTest)
library(emmeans)


format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

##Body weight ####
BW_data <- read_csv("../data/BW_compilate.csv")

# t test
  t.test(BW_data$starting_bw, BW_data$ending_bw, paired = TRUE)

  # Reshape to long format for plotting
  BW_long <- BW_data %>%
    pivot_longer(cols = c(starting_bw, ending_bw),
                 names_to = "timeframe",
                 values_to = "BW") %>% 
    mutate(Time = factor(timeframe, levels = c("starting_bw", "ending_bw"))) 

  ggplot(BW_long, aes(x = Time, y = BW, group = ID)) +
    geom_point(size = 3, color = "black", alpha = 0.1) +      # decrease alpha here
    geom_line(alpha = 0.5, color = "gray40") +
    stat_summary(aes(group = 1), fun = mean, geom = "line",       # add mean line
                 color = "black", size = 1.2, linetype = "solid") +
    stat_summary(aes(group = 1), fun = mean, geom = "point",      # add mean points
                 color = "black", size = 4) +
    labs(x = "Weeks", y = "Body weight (g)") +
    scale_x_discrete(labels = c("starting_bw" = "0", "ending_bw" = "12")) +
    format.plot

#echoMRI data import####

echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
 filter(Date <= "2025-02-24") %>% #DAY 1 OF RESTRICTION 
  filter(ID != 3715) %>% 
  arrange(Date) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index) %>% 
  group_by(ID) %>% 
  mutate(day_rel = n_measurement - first(n_measurement),
         date_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) 

###adiposity index####

# Calculate mean and SEM
adiposity_summary <- echoMRI_data %>%
  group_by(day_rel) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    sem_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n())
  )

plot_2 <- echoMRI_data %>%
  ggplot(aes(x = day_rel*3, y = adiposity_index)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = adiposity_summary,
              mapping = aes(x = day_rel*3, ymin = mean_adiposity - sem_adiposity, ymax = mean_adiposity + sem_adiposity),
              inherit.aes = FALSE, fill = "gray", alpha = 0.6) +
  geom_line(data = adiposity_summary,
           mapping = aes(x = day_rel*3, y = mean_adiposity),
            color = "black", size = 1.2) +
  labs(
    x = "weeks",
    y = "Adiposity Index (fat/lean)"
  )+
  format.plot
plot_2

###Lean mass####

# Calculate mean and SEM
lean_summary <- echoMRI_data %>%
  group_by(day_rel) %>%
  summarise(
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean = sd(Lean, na.rm = TRUE) / sqrt(n())
  )

plot_3 <- echoMRI_data %>%
  ggplot(aes(x = day_rel*3, y = Lean)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = lean_summary, 
              mapping = aes(x = day_rel*3, ymin = mean_lean - sem_lean, ymax = mean_lean + sem_lean), 
              inherit.aes = FALSE, fill = "gray", alpha = 0.6) +
  geom_line(data = lean_summary, 
            mapping = aes(x = day_rel*3, y = mean_lean), 
            color = "black", size = 1.2) +
  labs(
    x = "weeks",
    y = "Lean mass (g)"
  )+
  format.plot

plot_3


###Fat mass####

# Calculate mean and SEM
fat_summary <- echoMRI_data %>%
  group_by(day_rel) %>%
  summarise(
    mean_fat = mean(Fat, na.rm = TRUE),
    sem_fat = sd(Fat, na.rm = TRUE) / sqrt(n())
  )

plot_4 <- echoMRI_data %>%
  ggplot(aes(x = day_rel*3, y = Fat)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = fat_summary, 
              mapping = aes(x = day_rel*3, ymin = mean_fat - sem_fat, ymax = mean_fat + sem_fat), 
              inherit.aes = FALSE, fill = "gray", alpha = 0.6) +
  geom_line(data = fat_summary, 
            mapping = aes(x = day_rel*3, y = mean_fat), 
            color = "black", size = 1.2) +
  labs(
    x = "weeks",
    y = "Fat mass (g)"
  )+
  format.plot

plot_4


#TEE####
#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 
sable_tee_data <- sable_dwn %>% # Load the data
  mutate(numeric_index= as.numeric(str_extract(sable_idx,"[0-9]+"))) %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  filter(numeric_index < 11) %>% 
  mutate(SABLE  = if_else(numeric_index >= 8 , "After", "Before")) %>% 
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
  filter(!ID == 3715, is_complete_day ==1) %>% #3715 died and something weird happend with 3723 
   filter(!ID %in% c(3715,3723), is_complete_day ==1) %>% #3715 died and something weird happend with 3723 
 #   filter(!ID %in% c(3706,3707,3709,3712,3713, 3716,3717), is_complete_day ==1) %>% # those animals has weird data
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

sable_tee_data %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(
    numeric_index = as.numeric(str_extract(sable_idx, "[0-9]+"))
  ) %>%
  summarise(
    ID = ID,
    AFTER = max(numeric_index),
    BEFORE = min(numeric_index)
  ) %>%
  unique() %>%
  view()


# Create the plot tee 
# Add ID labels to the lowest TEE values
plot_5 <- sable_tee_data %>%
  mutate(SABLE = factor(SABLE, levels = c("Before", "After"), labels = c("0", "12"))) %>%
  ggplot(aes(SABLE, tee)) +
  geom_point(aes(group = ID), color = "gray", alpha = 0.3)+
  geom_line(aes(group = ID), alpha = 0.3) + # Lines connecting paired observations
  stat_summary(fun.data = mean_se, 
               geom = "ribbon", fill = "gray", color = NA, 
               aes(group = 1), alpha = 1) +
  stat_summary(fun = mean,
               geom = "line", color = "black", size = 1.2,
               aes(group = 1))+
  # Axis labels
  labs(x = "Weeks", y = "24h TEE") + 
  format.plot

plot_5

#FI data import####

FI_data <-read_csv("~/Documents/GitHub/data/data/FI.csv") #data import

FI_data <- FI_data %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(DATE <= "2025-02-24") %>% #DAY 1 OF RESTRICTION 
  filter(DIET == "LFD") %>%
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(ID != 3715)%>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  drop_na(INTAKE_GR) %>% 
  filter(corrected_intake_kcal <= 50)%>% # we have two animals that ate more than 75 kcal in one day so we believe is a systematic error 
  mutate(DATE_rel = DATE - first(DATE),
         cumulativeFI = cumsum(corrected_intake_kcal),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  mutate(relative_weeks = as.integer(DATE_rel / 7) + 1) %>% 
  filter(relative_weeks < 14) 
  

# Filter data to even-numbered weeks only
FI_data_even <- FI_data %>%
  filter(relative_weeks %% 2 == 0)

# Plot
plot_6 <- ggplot(FI_data_even, aes(x = as.factor(relative_weeks), y = cumulativeFI)) +
  geom_point(aes(group = ID), color = "gray", alpha = 0.3) +
  geom_line(aes(group = ID), alpha = 0.3, color = "gray") +
  stat_summary(fun.data = mean_se,
               geom = "ribbon", fill = "gray", color = NA,
               aes(group = 1), alpha = 1) +
  stat_summary(fun = mean,
               geom = "line", color = "black", size = 1.2,
               aes(group = 1)) +
  labs(x = "weeks", y = "food intake (kcal)") +
  format.plot
plot_6

# Create an arranged plot
ggarrange(plot_1, plot_5,plot_6, plot_4,plot_3, plot_2,
          nrow = 2, ncol=3,
          align = "hv", 
          labels = c("A", "B","C", "D", "E", "F"))


