#This script aim to explore adiposity index (fat/lean mass) in cohort 12 animals prior RTIOXA47 injections

#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggrepel)
#Data import####
echomri_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echomri_data <- echomri_data %>% 
  filter(COHORT ==12) %>% 
  select(ID, adiposity_index) %>% #good so here we checked there is 27 animals in total (cohort 11)
  mutate(GROUP = case_when(
    ID %in% c(8075, 8077, 8078) ~ "RTI_47",
    ID %in% c(8074, 8076, 8079) ~ "VEHICLE"
    )) 
# Summarize the data####
summary_data <- echomri_data %>%
  group_by(GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

#plot####
plot <- ggplot() +
  geom_point(data = echomri_data, aes(x = GROUP, y = adiposity_index), alpha = 0.7) +  # Optional: raw data points
  geom_point(data = summary_data, aes(x = GROUP, y = mean_adiposity), 
             color = "red", size = 3) +
  geom_errorbar(data = summary_data, aes(x = GROUP, 
                                         ymin = mean_adiposity - se_adiposity, 
                                         ymax = mean_adiposity + se_adiposity),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("Adiposity Index") +
  xlab("Group")
plot
# Check if adiposity index differs between groups using ANOVA####
model <- lm(adiposity_index ~ GROUP, data = echomri_data)
anova_result <- anova(model)
print(anova_result) #no diff among groups 