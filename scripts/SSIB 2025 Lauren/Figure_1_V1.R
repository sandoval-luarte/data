# We aim to test if BW is statistically different after 12 weeks with LFD in NZO female mice

#Raw data /Users/carosandovalcaballero/Documents/GitHub/data/data/COHORT_3.csv
#         /Users/carosandovalcaballero/Documents/GitHub/data/data/COHORT_4.csv
#         /Users/carosandovalcaballero/Documents/GitHub/data/data/COHORT_5.csv

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)

#bw data import####

bw_data <- read_csv("SSIB 2025 Lauren/BW_compilate.csv")

#paired t test####
t.test(bw_data$starting_bw, bw_data$ending_bw, paired = TRUE)
nrow(bw_data) #number of individuals

# Convert to long format####
bw_long <- bw_data %>%
  pivot_longer(cols = c(starting_bw, ending_bw),
               names_to = "timepoint",
               values_to = "body_weight") %>% 
  mutate(timepoint = factor(timepoint, levels = c("starting_bw", "ending_bw")))

#plot 1 bw before and after 12 weeks of LFD ####
bw_long %>% 
  ggplot(aes(x = timepoint, y = body_weight)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  labs(
    x = "timepoint",
    y = "body weight (g)"
  ) +  # White background
  theme_classic()  

