# We aim to test the rate of change of FI over 12 weeks of LFD in NZO mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)

mdl_data_fi <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study 
  # filter(DIET == "LFD") %>%  #We exclude chow data I dont know why still appear the three 3 dates that correspond to chow
  group_by(ID) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTION PERIOD FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  drop_na(corrected_intake_kcal) %>% 
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal)) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  select(
    ID, cumsum_kcal, time, DIET
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols="cumsum_kcal") %>% 
  select(
    ID, time, name, value
  ) 

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_fi <- lmerTest::lmer(
  data = mdl_data_fi,
  value ~ time + (1+time|ID)
)
summary(lmer_fi)
view(coef(lmer_fi))
coef_fi <- coef(lmer_fi)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_fi= `(Intercept)`)

write_csv(coef_fi, "coef_fi.csv")  # Save the dataframe as an csv file