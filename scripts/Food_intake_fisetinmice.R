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
 # filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID == 7866) %>% 
  ggplot(aes(
    DATE, BW, group = ID, color = as.factor(ID)
  )) +
  geom_point() +
  geom_line() 
bw

echo_full <- echomri_data %>% 
  filter(COHORT ==2)

A <-echo_full %>% 
  filter(DIET_CODE == "D12451i") %>% 
  filter(Date=="2025-02-11") %>% 
  ggplot(aes(Date,adiposity_index, color = as.factor(ID))) +
  geom_point() +
  geom_line()+
  geom_text(aes(label = ID), vjust = -0.5, size = 3, check_overlap = TRUE) +  # Add labels
  facet_wrap(DIET_CODE~SEX)
  A
  
  B <- echo_full %>% 
    filter(DIET_CODE == "D12451i") %>% 
    filter(Date=="2025-02-11") %>% 
    mutate(REGIMEN = if_else(ID %in% c(7873,7875,7876,7879,7860,7867,7864,7862), "ADLIB", "RESTRICTED")) %>% 
    select(c(ID, adiposity_index,DIET_CODE,REGIMEN,SEX))
  
  #check if the adiposity index between RESTRICTED HFD and ADLIB HFD are different using t-test (F)
  t_test_result <- lm(adiposity_index ~ REGIMEN, data = B %>% filter(SEX=="F"))
  # View the result
  print(t_test_result)
  # Summary
  summary(t_test_result)
  
  #check if the adiposity index between RESTRICTED HFD and ADLIB HFD are different using t-test (M)
  t_test_result_M <- lm(adiposity_index ~ REGIMEN, data = B %>% filter(SEX=="M"))
  # View the result
  print(t_test_result_M)
  # Summary
  summary(t_test_result_M)
    
###BW NZO MICE####

bw %>% 
 # filter(SEX =="M") %>% 
 # filter(DATE > "2025-02-01") %>% 
 filter((ID %in% c(3722,3723,3724,3725,3727,3728,3729))) %>% 
  ggplot(aes(
    DATE, BW, group = ID, color = as.factor(ID)
  )) +
  geom_point() +
  geom_line() #+
  #facet_wrap(~SEX)

####FI NZO MICE####
fi_plot_gr    <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter((ID %in% c(3722,3723,3724,3725,3727,3728,3729))) %>% 
  filter(corrected_intake_gr < 15) %>%  #TYPING MISTAKES
  filter(DATE > "2025-02-01") %>% 
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



