#This script aim to explore FI and BW over the weeks in orexin-cre and wt mice

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr) #to use drop_na()
library(lme4)

#FOOD INTAKE####

FI_data <-read_csv("~/Documents/GitHub/data/data/FI.csv") #data import

FI_data <- FI_data %>% 
  filter(COHORT ==11) %>% 
  select(ID,INTAKE_GR,corrected_intake_gr,corrected_intake_kcal,STRAIN,COHORT,SEX,DATE) %>% 
  drop_na()
  

  n_distinct(FI_data$ID) #here we know there is 27 animals in total n=18 orexin cre and 9 c57

  plot <- FI_data %>% 
    ggplot(aes(DATE, corrected_intake_kcal, group = ID, color = SEX)) +
    geom_point() +
    geom_line() +
    facet_wrap(STRAIN~SEX)+
    geom_text(data = FI_data,
              aes(label = ID), 
              hjust = -0.2, vjust = 0.5, 
              size = 3, show.legend = FALSE) #so check the FI for 315, 318, 314
#BW####
  
  BW_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import
  
  BW_data <- BW_data %>% 
    filter(COHORT ==11) 
  
  n_distinct(BW_data$ID) #here we know there is 27 animals in total n=18 orexin cre and 9 c57
  
  plot <- BW_data %>% 
    ggplot(aes(DATE, BW, group = ID, color = SEX)) +
    geom_point() +
    geom_line() +
    facet_wrap(STRAIN~SEX)+
    geom_text(data = BW_data,
              aes(label = ID), 
              hjust = -0.2, vjust = 0.5, 
              size = 3, show.legend = FALSE) #so check the FI for 315, 318, 314


  mixed_model_result <- lmer(BW ~ STRAIN * SEX + (1 | ID), data = BW_data)
  summary(mixed_model_result)
  anova(mixed_model_result) 
  
  #BODY COMPOSITION####
  
  echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import
  
  echo_data <-   echoMRI_data  %>% 
    filter(COHORT ==11) 
  
  n_distinct(echo_data$ID) #here we know there is 27 animals in total n=18 orexin cre and 9 c57
  
  # First, reshape your data to long format
  data_long <- echo_data  %>%
    pivot_longer(cols = c(adiposity_index, Lean, Fat),
                 names_to = "measurement",
                 values_to = "value")
  
  # Calculate the mean and standard error for each group and measurement
  summary_data <- data_long %>%
    group_by(STRAIN,SEX, measurement) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      n_value = n(),
      se_value = sd_value / sqrt(n_value)
    ) %>%
    ungroup()
  
  n_distinct(data_long$ID) #here we know there is 27 animals in total n=18 orexin cre and 9 c57
 
  # Create the bar plot with error bars
  ggplot(summary_data, aes(x = STRAIN, y = mean_value, fill = measurement)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      position = position_dodge(width = 0.9),
      width = 0.2
    ) +
    theme_minimal() +
    facet_wrap(SEX~ measurement, scales = "free_y")
  #It seems like orexin cre males has slightly more fat tissue
  
  # ANOVA for Adiposity Index
  anova_adiposity <- aov(value ~ STRAIN * SEX, data = filter(data_long, measurement == "adiposity_index"))
  summary(anova_adiposity) #only sex effects
  
  # ANOVA for Fat
  anova_fat <- aov(value ~ STRAIN * SEX, data = filter(data_long, measurement == "Fat"))
  summary(anova_fat) #only sex effects
  
  # ANOVA for Lean
  anova_lean <- aov(value ~ STRAIN * SEX, data = filter(data_long, measurement == "Lean"))
  summary(anova_lean) #only sex effects
 
  