#This script aim to explore adiposity index (fat/lean mass) in cohort 14 animals prior RTIOXA43 and RTIOXA47 injections for 5 days
#all drugs will be injected in a concentration of 10 mg-kg

#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggrepel)
#Data import####
echomri_data <-read_csv("~/Documents/GitHub/data/data/echoMRI.csv") #data import

echomri_data <- echomri_data %>% 
  filter(COHORT ==14) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
   mutate(DRUG = case_when(
     ID %in% c(41014,41016,41020, 41010, 41012, 41005, 41003) ~ "vehicle", #n=7 in each group
     ID %in% c(41006, 41015, 41008, 41021, 41001, 41018, 41007) ~ "RTI43",
     ID %in% c(41017, 41019, 41013, 41011, 41004, 41009, 41002) ~ "RTI47"    
     )) %>% 
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, DRUG) %>%
  mutate(
    day_rel = Date - first(Date))

 summary_data <- echomri_data %>%
   group_by(DRUG) %>%
   summarise(
     mean_AI = mean(adiposity_index, na.rm = TRUE), 
     se_AI = sd(adiposity_index, na.rm = TRUE) / sqrt(n()), 
     .groups = 'drop'
   )

#plot----
 plot <- ggplot() +
   geom_point(data = echomri_data,
              aes(x = DRUG, y = adiposity_index),
              alpha = 0.7) +
   geom_point(data = summary_data,
              aes(x = DRUG, y = mean_AI),
              color = "red", size = 3) +
   geom_errorbar(data = summary_data,
                 aes(x = DRUG,
                     ymin = mean_AI - se_AI,
                     ymax = mean_AI + se_AI),
                 width = 0.2, color = "red") +
   geom_text(data = echomri_data,
             aes(x = DRUG, y = adiposity_index, label = ID),
             vjust = -0.5,
             size = 3) +
   theme_minimal() +
   ylab("Adiposity index (fat/lean mass at baseline)") +
   xlab("Drug")
 
 plot
  
# Check if adiposity index differs between groups using ANOVA#### 

#check assumptions first
 
 model <- lm(adiposity_index ~ DRUG, data = echomri_data)
 
 par(mfrow = c(1, 2))
 plot(model, which = 1)  # Residuals vs fitted → variance, Residuals are centered around 0
 plot(model, which = 2)  # Q–Q plot → normality, Residuals are approximately normally distributed
 
 #Model assumptions were assessed by visual inspection of residual and Q–Q plots
 #which indicated approximate normality and homogeneity of variance. 
 #Group differences were analyzed using one-way ANOVA.
 
  
 model <- lm(adiposity_index ~ DRUG, data = echomri_data)
 anova_result <- anova(model)
 print(anova_result) #no diff among groups 

