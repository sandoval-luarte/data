#This script aim to explore adiposity index (fat/lean mass) in wt mice
#Our aim is to test no differences in adiposity index among 3 experimental group
#group 1 n=3: vehice (DMSO 2% and 98% corn oil) IP x 5 days 
#group 2 n=3: RTIOXA-43 30 mg/kg IP (Yanan)
#group 3 n=3: RTIOXA-43 30 mg/kg IP (medchemexpress)

library(ggplot2)
library(readr)
library(dplyr)

echomri_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import 

#la cagada queda cuando voy a procesar el echomri_open_files de data proces.R, no importa los animales 1001-1, 1002-1, 1003-1

echomri_data <- echomri_data %>% 
  filter(COHORT ==9) %>% #why here I can not identify the 1001-1
  select(ID, adiposity_index) %>% 
  mutate(GROUP = if_else(ID %in% c(8101,8099,8102), "rtioxa_43", "control"))
  
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
  xlab("Virus Group")
plot
#check if the adiposity index between strains are different using t-test
t_test_result <- lm(adiposity_index ~ GROUP, data = echomri_data)
# View the result
print(t_test_result)
# Summary
summary(t_test_result)
