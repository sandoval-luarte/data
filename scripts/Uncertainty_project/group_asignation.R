#This script aim to explore adiposity index (fat/lean mass) in orexin-cre and wt mice
#we aim split the orexin cre males into two groups: the first will receive the inhibitory DREADD and the second one will receive contro (mcherry)
#due technical issues not all the orexin cre males could be measured in the echoMRI because the piercing issue

library(ggplot2)
library(readr)
library(dplyr)

echomri_data <-read_csv("../data/echomri.csv") #data import

echomri_data <- echomri_data %>% 
  filter(COHORT ==11) %>% 
  select(ID, adiposity_index,SEX,STRAIN) %>% 
  filter(SEX=="M") %>% 
  filter(STRAIN=="OREXIN-CRE") %>% 
  mutate(VIRUS_GROUP = if_else(ID %in% c(315,320,318), "INHIBITORY_DREADD", "CONTROL"))
  
#In total we have 12 males orexin cre, however here appear just 7 so 5 animals are pending in the echoMRI

# Summarize the data
summary_data <- echomri_data %>%
  group_by(VIRUS_GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Create the plot
plot <- ggplot() +
  geom_point(data = echomri_data, aes(x = VIRUS_GROUP, y = adiposity_index), alpha = 0.7) +  # Optional: raw data points
  geom_point(data = summary_data, aes(x = VIRUS_GROUP, y = mean_adiposity), 
             color = "red", size = 3) +
  geom_errorbar(data = summary_data, aes(x = VIRUS_GROUP, 
                                         ymin = mean_adiposity - se_adiposity, 
                                         ymax = mean_adiposity + se_adiposity),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("Adiposity Index") +
  xlab("Virus Group")
plot
#check if the adiposity index between strains are different using t-test
t_test_result <- lm(adiposity_index ~ VIRUS_GROUP, data = echomri_data)
# View the result
print(t_test_result)
# Summary
summary(t_test_result)
