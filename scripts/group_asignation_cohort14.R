#This script aim to explore adiposity index (fat/lean mass) in cohort 14 animals prior RTIOXA43 and RTIOXA47 injections for 5 days
#all drugs will be injected in a concentration of 10 mg-kg

#FOR NOW I WILL EVALUATE WHATS GOING ON WITH THE BW


#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggrepel)
#Data import####
echomri_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import

echomri_data <- echomri_data %>% 
  filter(COHORT ==14) %>% 
 # select(ID, adiposity_index) %>% #good so here we checked there is 21 animals in total (cohort 14)
 select(ID, DATE, BW) %>% #good so here we checked there is 27 animals in total (cohort 11)
#filter(DATE=="2025-12-15") %>% 
   mutate(GROUP = case_when(
     ID %in% c(41004, 41002, 41020, 41016, 41018, 41021, 41008) ~ "vehicle",
     ID %in% c(41005, 41012, 41014, 41006, 41017, 41010, 41001) ~ "RTI43",
     ID %in% c(41007, 41013, 41003, 41019, 41015, 41011, 41009) ~ "RTI47"    
     )) 

 summary_data <- echomri_data %>%
   group_by(GROUP) %>%
   summarise(
     mean_BW = mean(BW, na.rm = TRUE), #HERE CHANGE BW FOR AI
     se_BW = sd(BW, na.rm = TRUE) / sqrt(n()), #HERE CHANGE BW FOR AI
     .groups = 'drop'
   )

#plot----
 plot <- ggplot() +
   geom_point(data = echomri_data, aes(x = GROUP, y = BW), alpha = 0.7) + #HERE CHANGE BW FOR AI
   geom_point(data = summary_data, aes(x = GROUP, y = mean_BW), #HERE CHANGE BW FOR AI
              color = "red", size = 3) +
   geom_errorbar(data = summary_data, aes(x = GROUP, 
                                          ymin = mean_BW - se_BW, #HERE CHANGE BW FOR AI
                                          ymax = mean_BW + se_BW),#HERE CHANGE BW FOR AI
                 width = 0.2, color = "red") +
   theme_minimal() +
   ylab(" BW in g") + #HERE CHANGE BW FOR AI
   xlab("Group")
  plot
  
  
# Check if BW differs between groups using ANOVA#### #HERE CHANGE BW FOR AI
  
 model <- lm(BW ~ GROUP, data = echomri_data)
 anova_result <- anova(model)
 print(anova_result) #no diff among groups 


ggplot(echomri_data, aes(x = DATE, y = BW)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ ID, scales = "free_y") +
  labs(x = "Date", y = "Body Weight (BW)") +
  theme_minimal()

