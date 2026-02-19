#This script aim to explore changes in BW after DREADD injection in orexin-cre mice


#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(car)

#BW ANALYSIS----
BW_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import

BW_data <- BW_data %>% 
  filter(COHORT %in% c(11, 17)) %>% 
  arrange(DATE) %>% 
  mutate(GROUP = case_when(
    ID %in% c(307, 316, 319, 522, 524) ~ "CONTROL_CERT", #5
    ID %in% c(297, 313, 318, 320,520,523) ~ "INH_DREADD_UNCERT", #6
    ID %in% c(321, 323, 325, 519, 525) ~ "CONTROL_UNCERT", #5
    ID %in% c(305, 306, 314, 322) ~ "WT" #4
  )) %>% 
  group_by(ID) %>% 
  mutate(
    surgery_date = min(
      DATE[COMMENTS %in% c(
        "BILATERAL_CONTROL_DREADD_INJECTION_SURGERY",
        "BILATERAL_INH_DREADD_INJECTION_SURGERY"
      )],
      na.rm = TRUE
    )
  ) %>% 
  filter(DATE >= surgery_date) %>% 
  drop_na(GROUP) %>% 
  mutate(day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))) %>% 
  filter(GROUP != "WT") 
BW_data  %>% 
  group_by(GROUP) %>%
  summarise(n_ID = n_distinct(ID)) 


##BW plot####
plot <- ggplot(BW_data, aes(x = day_rel, y = BW)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = ID), alpha = 0.7) +
  theme_minimal() +
  ylab("Body weight (g)") +
  xlab("Date") +
  facet_wrap(~ID)

plot

#BODY COMPOSITION ANALYSIS----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(11,17)) %>% 
  mutate(GROUP = case_when(
    ID %in% c(307, 316, 319, 522, 524) ~ "CONTROL_CERT", #5
    ID %in% c(297, 313, 318, 320,520,523) ~ "INH_DREADD_UNCERT", #6
    ID %in% c(321, 323, 325, 519, 525) ~ "CONTROL_UNCERT", #5
    ID %in% c(305, 306, 314, 322) ~ "WT" #4
  )) %>% 
  group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index,COHORT,DIET_FORMULA,n_measurement,GROUP) %>% 
  drop_na(GROUP) %>% 
  filter(GROUP != "WT") 
  
 
echoMRI_data  %>% 
  group_by(GROUP) %>%
  summarise(n_ID = n_distinct(ID))  #ok this is good

echoMRI_data_summary <- echoMRI_data  %>%
  ungroup() %>% 
  group_by(GROUP) %>%
  summarise(
    mean_ai = mean(adiposity_index, na.rm = TRUE),
    sem_ai  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )



ggplot() +
  geom_col(
    data = echoMRI_data_summary,
    aes(x = GROUP, y = mean_ai),
    width = 0.7
  ) +
  geom_errorbar(
    data = echoMRI_data_summary,
    aes(x = GROUP, ymin = mean_ai - sem_ai, ymax = mean_ai + sem_ai),
    width = 0.2
  ) +
  geom_jitter(
    data = echoMRI_data,
    aes(x = GROUP, y = adiposity_index),
    width = 0.1,
    size = 2
  ) +
  geom_text(
    data = echoMRI_data,
    aes(x = GROUP, y = adiposity_index, label = ID),
    vjust = -0.7,
    size = 3
  ) +
  labs(
    x = "Group",
    y = "Adiposity index (mean ± SEM)"
  ) +
  theme_classic()

shapiro.test(residuals(anova_model))

leveneTest(adiposity_index ~ GROUP, data = echoMRI_data)


#### STATS for group assignation ----
anova_model <- aov(adiposity_index ~ GROUP, data = echoMRI_data)
summary(anova_model)
# so groups are not different accordingly to adiposity index


