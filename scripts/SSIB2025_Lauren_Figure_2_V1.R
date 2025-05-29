# We aim to test the changes in body comp in NZO mice after 12 weeks of LFD

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models
library(ggpubr)  

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))
#echoMRI data import####

echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
 filter(Date <= "2025-02-2024") %>% #DAY 1 OF RESTRICTION 
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

plot_1 <- echoMRI_data %>%
  ggplot(aes(x = day_rel*3, y = adiposity_index)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = adiposity_summary, 
              mapping = aes(x = day_rel*3, ymin = mean_adiposity - sem_adiposity, ymax = mean_adiposity + sem_adiposity), 
              inherit.aes = FALSE, fill = "gray70", alpha = 0.6) +
  geom_line(data = adiposity_summary, 
           mapping = aes(x = day_rel*3, y = mean_adiposity), 
            color = "black", size = 1.2) +
  labs(
    x = "Weeks",
    y = "Adiposity Index"
  )+
  format.plot

plot_1

first_last_measurements <- echoMRI_data %>%
  arrange(ID, day_rel) %>%
  group_by(ID) %>%
  filter(row_number() == 1 | row_number() == n()) %>%
  ungroup()

first_last_measurements <- first_last_measurements %>%
  mutate(day_rel = ifelse(day_rel == 4, 3, day_rel))

# Calculate mean and SEM
adiposity_summary <- first_last_measurements%>%
  group_by(day_rel) %>%
  summarise(
    n = n(),  # number of observations
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    sem_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n())
  )

plot_2 <- first_last_measurements %>%
  ggplot(aes(x = day_rel, y = adiposity_index)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = adiposity_summary, 
              aes(x = day_rel, ymin = mean_adiposity - sem_adiposity, ymax = mean_adiposity + sem_adiposity), 
              inherit.aes = FALSE, fill = "gray70", alpha = 0.6) +
  geom_line(data = adiposity_summary, 
            aes(x = day_rel, y = mean_adiposity), 
            color = "black", size = 1.2) +
  scale_x_continuous(
    breaks = c(0, 3), 
    labels = c("0", "12")
  ) +
  labs(
    x = "Weeks",
    y = "Adiposity Index")+
      format.plot
   
plot_2


###Lean amass####

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
              inherit.aes = FALSE, fill = "gray70", alpha = 0.6) +
  geom_line(data = lean_summary, 
            mapping = aes(x = day_rel*3, y = mean_lean), 
            color = "black", size = 1.2) +
  labs(
    x = "Weeks",
    y = "Lean mass (g)"
  )+
  format.plot

plot_3

# Calculate mean and SEM
lean_summary <- first_last_measurements%>%
  group_by(day_rel) %>%
  summarise(
    n = n(),  # number of observations
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean = sd(Lean, na.rm = TRUE) / sqrt(n())
  )

plot_4 <- first_last_measurements %>%
  ggplot(aes(x = day_rel, y = Lean)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = lean_summary, 
              aes(x = day_rel, ymin = mean_lean - sem_lean, ymax = mean_lean + sem_lean), 
              inherit.aes = FALSE, fill = "gray70", alpha = 0.6) +
  geom_line(data = lean_summary, 
            aes(x = day_rel, y = mean_lean), 
            color = "black", size = 1.2) +
  scale_x_continuous(
    breaks = c(0, 3), 
    labels = c("0", "12")
  ) +
  labs(
    x = "Weeks",
    y = "Lean mass (g)")+
  format.plot

plot_4

# Create an arranged plot
ggarrange(plot_1, plot_2, plot_3, plot_4, 
          nrow = 1, 
          ncol = 2,
          align = "hv", 
          widths = c(0.7, 0.3),  # Increase width of panel A, decrease B
          labels = c("A", "B","C", "D"))


