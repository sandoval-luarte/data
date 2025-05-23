# We aim to test the rate of change of FI over 12 weeks of LFD in NZO mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models

#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))


fi_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>% 
  group_by(ID) %>%
  drop_na(corrected_intake_kcal) %>% 
  select(ID,DATE,corrected_intake_kcal) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal),
         day_rel = DATE - first(DATE)) %>% 
  mutate(relative_weeks = as.integer(day_rel / 7) + 1) %>% 
  filter(relative_weeks <=13)

# Custom function for mean ± SEM
mean_se <- function(x) {
  se <- sd(x) / sqrt(length(x))
  return(c(y = mean(x), ymin = mean(x) - se, ymax = mean(x) + se))
}

plot_1 <- fi_data %>%
  ggplot(aes(x = relative_weeks, y = cumsum_kcal)) +
  geom_point(aes(group = ID), alpha = 0.1) +             # individual points
  geom_line(aes(group = ID), alpha = 0.1) +              # individual lines
  stat_summary(fun.data = mean_se, 
               geom = "ribbon", fill = "grey70", color = NA, 
               aes(group = 1), alpha = 0.7) +
  stat_summary(fun = mean, 
               geom = "line", color = "black", 
               aes(group = 1)) +
  labs(
    x = "Weeks",
    y = "cumulative intake in kcal"
  ) +
  format.plot +
 scale_y_continuous(limits = c(0, 800)) +
 scale_x_continuous(breaks = seq(3, 12, by = 3)) 
plot_1 


# #Backup####
# fi_data <- read_csv("../data/FI.csv") %>% 
# filter(COHORT %in% c(3, 4, 5)) %>% 
#   filter(ID != 3715) %>% 
#   group_by(ID) %>%
#   mutate(start_date = ifelse(COMMENTS == "DAY_7_SABLE_DAY_1_LFD", DATE, NA)) %>%
#   fill(start_date, .direction = "down") %>%  # Propagate the date down within each ID group
#   filter(DATE >= start_date) %>%
#   select(-start_date) %>%  # Remove the helper column
#   drop_na(corrected_intake_kcal) %>% 
#   select(ID,DATE,corrected_intake_kcal) %>% 
#   filter(is.finite(corrected_intake_kcal)) %>%  #to eliminate inf values 
#   mutate(
#     ID = as.factor(ID),
#     time = as.numeric(DATE - min(DATE))
#   ) %>% 
#   filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTION PERIOD FROM THE ANALYSIS
#   filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
#   mutate(cumsum_kcal= cumsum(corrected_intake_kcal)) %>% 
#   mutate(
#     time = scale(time) # this is for model convergence
#   ) %>% 
#   ungroup() %>% 
#   pivot_longer(cols="cumsum_kcal") %>% 
#   select(
#     ID, time, name, value
#   ) 
# # Fit a linear mixed model: cumulative intake ~ time with random effect for ID
# model <- lmer(value ~ time + (1 | ID), data = fi_data)
# # Summary of the model
# summary(model)
# 
# #The model estimated the rate of food intake as 156.99 kcal/day across all animals. 
# #Given that we analyzed 23 animals, this corresponds to an average increase of 6.82 kcal/day per mouse (SE = 0.76, t = 205.9).
# 

#echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>% 
  arrange(Date) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index) %>% 
  mutate(day_rel = n_measurement - first(n_measurement),
         date_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean)) %>% 
 # mutate(relative_weeks = as.integer(date_rel / 7) + 1) %>% 
  filter(date_rel <=100)

###adiposity index####

# Calculate mean and SEM
adiposity_summary <- echoMRI_data %>%
  group_by(day_rel) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    sem_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n())
  )

# Plot with mean and SEM
plot_1 <- echoMRI_data %>%
  ggplot(aes(x = day_rel, y = adiposity_index)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  geom_ribbon(data = adiposity_summary, 
              mapping = aes(x = day_rel, ymin = mean_adiposity - sem_adiposity, ymax = mean_adiposity + sem_adiposity), 
              inherit.aes = FALSE, fill = "gray70", alpha = 0.6) +
  geom_line(data = adiposity_summary, 
            mapping = aes(x = day_rel, y = mean_adiposity), 
            color = "black", size = 1.2) +
  labs(
    x = "measurement",
    y = "Adiposity Index"
  ) +
  format.plot +
#  scale_x_continuous(breaks = seq(0, 12, by = 2.4)) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE) 
  

plot_1

