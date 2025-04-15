#This script aim to explore adiposity index (fat/lean mass) in wt mice
#Our aim is to test no differences in adiposity index among 3 experimental group
#group 1 n=3: vehice (DMSO 2% and 98% corn oil) IP x 5 days and 2 days off
#group 2 n=3: RTIOXA-43 30 mg/kg IP (Yanan)
#group 3 n=3: RTIOXA-43 30 mg/kg IP (medchemexpress)

library(ggplot2)
library(readr)
library(dplyr)

echomri_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import 

echomri_data <- echomri_data %>% 
  filter(COHORT ==9) %>% 
  select(ID, adiposity_index, n_measurement) %>% 
  mutate(GROUP = case_when(
    ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
    ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
    ID %in% c(8095, 8100, 8103) ~ "CONTROL")) %>% 
  filter(n_measurement ==2)
  
# Summarize the data
summary_data <- echomri_data %>%
  group_by(GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Create the plot
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
# Check if adiposity index differs between groups using ANOVA
model <- lm(adiposity_index ~ GROUP, data = echomri_data)
anova_result <- anova(model)
print(anova_result)
