# This script aims to explore changes in body weight (BW) and adiposity index 
# (AI, fat/lean mass) in females and males C57BL6J fed with HFD (D12451i, research diets)

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)


 # BW over time data----

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(STRAIN=="C57BL/6J") %>% #JUST C57BL6J %>% 
  filter(DIET_FORMULA =="D12451i") %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) %>% 
  filter(day_rel <=100) %>% 
  mutate(
    age = factor(
      ifelse(COHORT == 7, "old-72 week old", "young-8 week old"),
      levels = c("young-8 week old", "old-72 week old")
    )
  ) %>% 
  ungroup() 

## summary ----

BW_binned <- BW_data %>%
  mutate(day_bin = floor(day_rel / 10) * 10) #showing data weekly 

BW_summary <- BW_binned %>%
  group_by(SEX, age, day_bin) %>%
  summarise(
    mean_BW_gain = mean(bw_rel),
    sem_BW_gain  = sd(bw_rel) / sqrt(n_distinct(ID)),
    mean_BW = mean(BW),
    sem_BW  = sd(BW) / sqrt(n_distinct(ID)),
    .groups = "drop"
  )

 # BW raw plot----

ggplot() +
  geom_line(
    data = BW_data,
    aes(x = day_rel, y = BW, group = ID),
    alpha = 0.25,
    color = "grey60"
  ) +
  geom_ribbon(
    data = BW_summary,
    aes(
      x = day_bin,
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW,
      fill = age
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = BW_summary,
    aes(x = day_bin, y = mean_BW, color = age),
    linewidth = 1.2
  ) +
  facet_wrap(~ SEX * age) +
  theme_classic() +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight (g)",
    color = "Age",
    fill = "Age"
  )+
  scale_x_continuous(limits = c(0, 80))


# BW % gain plot ----

ggplot() +
  geom_line(
    data = BW_data,
    aes(x = day_rel, y = bw_rel, group = ID),
    alpha = 0.25,
    color = "grey60"
  ) +
  geom_ribbon(
    data = BW_summary,
    aes(
      x = day_bin,
      ymin = mean_BW_gain - sem_BW_gain,
      ymax = mean_BW_gain + sem_BW_gain,
      fill = age
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = BW_summary,
    aes(x = day_bin, y = mean_BW_gain, color = age),
    linewidth = 1.2
  ) +
  facet_wrap(~ SEX * age) +
  theme_classic() +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight gain (%)",
    color = "Age",
    fill = "Age"
  )+
  scale_x_continuous(limits = c(0, 80))

# adiposity index over time ----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(STRAIN=="C57BL/6J") %>% #JUST C57BL6J %>% 
  filter(DIET_FORMULA =="D12451i") %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index,COHORT,SEX) %>%
  mutate(
    day_rel = Date - first(Date),
      ai_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
    fat_rel = 100 * (Fat - first(Fat)) / first(Fat)) %>% 
  filter(day_rel <= 100) %>% 
  mutate(
    age = factor(
      ifelse(COHORT == 7, "old-72 week old", "young-8 week old"),
      levels = c("young-8 week old", "old-72 week old")
    )
  ) %>% 
  ungroup() 

## summary ----

ai_binned <- echoMRI_data %>%
  mutate(day_bin = floor(day_rel / 7) * 7) #showing data weekly 

ai_summary <- ai_binned %>%
  group_by(SEX, age, day_bin) %>%
  summarise(
    mean_ai_gain = mean(ai_rel),
    sem_ai_gain  = sd(ai_rel) / sqrt(n_distinct(ID)),
    mean_ai = mean(adiposity_index),
    sem_ai  = sd(adiposity_index) / sqrt(n_distinct(ID)),
    .groups = "drop"
  )

# AI raw plot----

ggplot() +
  geom_line(
    data = echoMRI_data ,
    aes(x = day_rel, y = adiposity_index, group = ID),
    alpha = 0.25,
    color = "grey60"
  ) +
  geom_ribbon(
    data = ai_summary,
    aes(
      x = day_bin,
      ymin = mean_ai - sem_ai,
      ymax = mean_ai + sem_ai,
      fill = age
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = ai_summary,
    aes(x = day_bin, y = mean_ai, color = age),
    linewidth = 1.2
  ) +
  facet_wrap(~ SEX * age) +
  theme_classic() +
  labs(
    x = "Days relative to first measurement",
    y = "adiposity index (fat/lean mass)",
    color = "Age",
    fill = "Age"
  )+
  scale_x_continuous(limits = c(0, 90))

# AI % gain plot ----

ggplot() +
  geom_line(
    data = echoMRI_data,
    aes(x = day_rel, y = ai_rel, group = ID),
    alpha = 0.25,
    color = "grey60"
  ) +
  geom_ribbon(
    data = ai_summary,
    aes(
      x = day_bin,
      ymin = mean_ai_gain - sem_ai_gain,
      ymax = mean_ai_gain + sem_ai_gain,
      fill = age
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = ai_summary,
    aes(x = day_bin, y = mean_ai_gain, color = age),
    linewidth = 1.2
  ) +
  facet_wrap(~ SEX * age) +
  theme_classic() +
  labs(
    x = "Days relative to first measurement",
    y = "adiposity index gain (%)",
    color = "Age",
    fill = "Age"
  )+
  scale_x_continuous(limits = c(0, 90))


