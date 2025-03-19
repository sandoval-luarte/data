# We aim to compare the change in adiposity index after 12 weeks with LFD in NZO female mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)

#echoMRI analysis####
echo_data <- read_csv("~/Documents/Github/data/scripts/SSIB_2025_Lauren/echoMRI.csv") %>%
  filter(COHORT %in% c(3, 4, 5)) %>%
  filter(ID != 3715) %>%
  select(Date, ID, adiposity_index, Fat, Lean) %>%
  filter(Date %in% as.Date(c("2024-11-12", "2024-11-19", "2024-11-26", "2024-12-03", 
                             "2025-02-07", "2025-02-20"))) %>%
  filter(!(ID %in% c(3722, 3723, 3724, 3725) & !(Date %in% as.Date(c("2025-02-07", "2024-11-26"))))) %>% 
  mutate(status = if_else(Date %in% as.Date(c("2024-11-12", "2024-11-19", "2024-11-26", "2024-12-03")), 
                          "before", "after")) %>% 
  filter(!(ID %in% c(3727, 3728, 3729) & Date == as.Date("2025-02-20"))) %>% 
  select(ID, adiposity_index, status) 

#counting rows####
echo_data %>%
  group_by(ID) %>%
  summarize(status_count = n_distinct(status)) %>%
  filter(status_count == 2) %>%
  nrow() #So here we have 23 animals that contain the status before and after

id_counts <- echo_data %>%
  group_by(ID) %>%
  summarize(row_count = n())
id_counts  #double check the number of data per mouse , we should have 2 per animal

# Convert long format to wide format
echo_data_wide <- echo_data %>%
  pivot_wider(names_from = status, values_from = adiposity_index) 

t.test(echo_data_wide$after, echo_data_wide$before, paired = TRUE)

#plot####

echo_data <- echo_data %>%
  mutate(status = factor(status, levels = c("before", "after")))

echo_data %>%
    ggplot(aes(x = status, y = adiposity_index)) + 
    geom_point(aes(group = ID), alpha = 0.5) +
    geom_line(aes(group = ID), alpha = 0.5) +
    theme_classic() 
