# Data analysis of the adiposity index (fat mass/lean mass) of NZO female mice assigned to two diet regimens: 
# "ad lib" or "restricted" after 9-12 weeks on LFD (D12450Ki, Research Diets).

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr)  



#----FOOD INTAKE----
# My idea is to take the mean of the last 3 days of daily food intake to restrict the animals to ~60% of the obese energy intake
# A complementary idea is to take only the IDs that are assigned to the "restricted" group.

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 6) #Just fisetin mice 
  

FI_data_ <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5))
FI_data

fi_plot <- FI_data %>% 
  # filter(ID %in% c(3710, 3725,3708,3729,3724,3723,3728,3714,3727,3721,3722,3720)) %>% 
  ggplot(aes(x = DATE, y = corrected_intake_kcal, group =as.factor(ID))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(color = ID)) +  # Add color for better distinction
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    size = 0.5,
    shape = 21,
    color = "black",
    fill = "red",
    position = position_dodge(width = 0.2)  # Mean and SEM as big points
  ) +
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  labs(
    title = "Corrected Intake Over Days",
    x = "Date",
    y = "Corrected Intake (kcal)"
  ) +
  theme_minimal() 

fi_plot_gr <- FI_data_ %>% 
  # drop_na() %>% 
  #  filter(ID %in% c(3710, 3725,3708,3729,3724,3723,3728,3714,3727,3721,3722,3720)) %>% 
  filter(DATE > "2025-01-20") %>% 
  ggplot(aes(x = DATE, y = corrected_intake_gr, group =as.factor(ID))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(color = as.factor(ID))) +  # Add color for better distinction
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  labs(
    title = "Corrected Intake Over Days",
    x = "Date",
    y = "Corrected Intake (gr.)"
  ) +
  theme_minimal() 
fi_plot_gr    

bw <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(3, 4, 5))
bw

echo_full <- echomri_data %>% 
  filter(COHORT %in% c(3, 4, 5))

echo_full %>% 
  ggplot(aes(Date, Lean/Weight, color = as.factor(ID))) +
  geom_point() +
  geom_line()

bw %>% 
  ggplot(aes(
    DATE, BW, group = ID, color = as.factor(ID)
  )) +
  geom_point() +
  geom_line()
