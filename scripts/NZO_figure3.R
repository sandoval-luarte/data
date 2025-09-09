# This script aims to explore changes fat mass in middle age NZO after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      Date == as.Date("2025-02-20") ~ "peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))
#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        #   strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))


plot <- echoMRI_data  %>%
  ggplot(aes(x = STATUS, y = Fat, color = DRUG, group = ID)) +
  
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
  labs(y = "Fat mass in grams", color = "Drug", fill = "Drug") +
  facet_wrap(~GROUP) +
  format.plot+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot




