# Data analysis of the adiposity index (fat mass/lean mass) of NZO female mice assigned to two diet regimens: 
# "ad lib" or "restricted" after 9-12 weeks on LFD (D12450Ki, Research Diets).

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr)  

#----ADIPOSITY INDEX----

#open echo_mri_data
echomri_data <- read_csv("../data/echomri.csv")

echomri_data_2<- echomri_data %>%
  filter(COHORT > 2 & COHORT < 6) %>%
  filter(Date == "2025-01-15") %>%
  select(ID, adiposity_index) %>%
  mutate(diet_group = if_else(ID %in% c(
    3711, 3718, 3707, 3719, 3709,
    3726, 3712, 3716, 3713, 3717, 3706
   ), "AD_LIB", "RESTRICTED"))  #%>%
#   group_by(diet_group) %>%
#   summarize(animal_count = n(), .groups = "drop")
# 
 echomri_data_3<-echomri_data_2 %>%
  ggplot(aes(x = diet_group, y = adiposity_index, color = ID)) +
  geom_point(size = 3, alpha = 0.8) + # Individual points
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +    # Labels for each point
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    size = 1.5,
    shape = 21,
    color = "black",
    fill = "red",
    position = position_dodge(width = 0.2)     # Mean and SEM as big points
  ) +
  theme_minimal() +     # Aesthetic improvements
  labs(
    title = "Adiposity Index by Diet Group",
    x = "Diet Group",
    y = "Adiposity Index",
    color = "Animal ID"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
echomri_data_3

#check if the adiposity index between AD LIB and RESTRICTED groups are different using t-test
t_test_result <- lm(adiposity_index ~ diet_group, data = echomri_data_2)
# View the result
print(t_test_result)
# Summary
summary(t_test_result)

#----FOOD INTAKE----
# My idea is to take the mean of the last 3 days of daily food intake to restrict the animals to ~60% of the obese energy intake
# A complementary idea is to take only the IDs that are assigned to the "restricted" group.

FI_data <- read_csv("../data/FI.csv") %>% 
    filter(COHORT > 2 & COHORT < 6) %>%
    filter(DATE > "2025-01-1") 

fi_plot <- FI_data %>% 
    filter(ID %in% c(3710, 3725,3708,3729,3724,3723,3728,3714,3727,3721,3722,3720)) %>% 
    ggplot(aes(x = DATE, y = corrected_intake_kcal, group =as.factor(ID))) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(aes(color =ID)) +  # Add color for better distinction
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
    