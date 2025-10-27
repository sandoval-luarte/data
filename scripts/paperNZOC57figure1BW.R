# This script aims to explore changes in body weight in NZO and c57 after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#libraries----
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
library(patchwork)


#BW CSV data import RTIOXA 47 ----
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 
                7876, 7879, 7880, 7881,7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,
                7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871,
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ),
    day_rel = DATE - first(DATE)
  ) %>%
  mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 202 & GROUP == "restricted" ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 183 & GROUP == "restricted"~ "peak obesity",
      STRAIN == "C57BL6/J" & day_rel == 204 & GROUP == "ad lib" ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 241 & GROUP == "restricted"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 218 & GROUP == "restricted"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 239 & GROUP == "ad lib"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 218 & GROUP == "ad lib"~ "BW loss",
      STRAIN == "NZO/HlLtJ" & day_rel == 161 ~ "BW loss",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance", 
      STRAIN == "NZO/HlLtJ" &  COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    ) ) #%>% 
  #filter(STRAIN == "NZO/HlLtJ" | (STRAIN == "C57BL6/J" & DIET_FORMULA == "D12451i"))

BW_data %>% 
  group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO and 24 C57 

write_csv(BW_data, "../data/BW_data.csv") # Save as CSV


highlight_data <- BW_data %>%
  filter(STATUS == "peak obesity")

highlight_data_2 <- BW_data %>%
  filter(STATUS == "BW loss")

highlight_data_3 <- BW_data %>%
  filter(STATUS == "BW maintenance")

highlight_data_4 <- BW_data %>%
  filter(STATUS == "BW regain")

#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14))


#plot 1 description of BW over time per ID----
plot <- BW_data %>% 
  ggplot(aes(day_rel, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~STRAIN*GROUP) +
  #geom_smooth() +
  labs(
    x = "Days relative to baseline",
    y = "BW (grams)"
  ) +
  geom_vline(data = highlight_data, aes(xintercept = day_rel),
             linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(data = highlight_data_2, aes(xintercept = day_rel),
             linetype = "dashed", color = "green", alpha = 0.7) +
  geom_vline(data = highlight_data_3, aes(xintercept = day_rel),
             linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_vline(data = highlight_data_4, aes(xintercept = day_rel),
             linetype = "dashed", color = "orange", alpha = 0.7) +
  format.plot

plot

# Make STATUS an ordered factor
BW_data <- BW_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  filter(!is.na(STATUS))

BW_data %>%
  group_by(STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this is good


# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX, DIET_FORMULA
summary_data <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

#plot 2 BW in each STATUS collapsed by DIET_FORMULA----
ggplot() +
  # Individual points, smaller and more transparent
  geom_point(data = BW_data, 
             aes(x = STATUS, y = BW, color = SEX, group = ID),
             size = 1.2, alpha = 0.3) +
  geom_line(data = BW_data, 
            aes(x = STATUS, y = BW, color = DRUG, group = ID),
            alpha = 0.15) +
  
  # Mean ± SEM as shaded ribbon
  geom_ribbon(data = summary_data,
              aes(x = STATUS, ymin = mean_BW - sem_BW, ymax = mean_BW + sem_BW,
                  group = interaction(GROUP, SEX), fill = SEX),
              alpha = 0.2) +
  geom_line(data = summary_data,
            aes(x = STATUS, y = mean_BW, group = interaction(GROUP, SEX), color = SEX),
            size = 1) +
  
  # Facet by STRAIN (rows) × GROUP  (columns)
  facet_grid(~ STRAIN*GROUP, scales = "free_x") +
  
  scale_color_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  scale_fill_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  
  labs(
    x = "Status",
    y = "Body weight (grams)",
    color = "Sex",
    fill = "Sex"
  ) +
  format.plot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX using bw_rel
summary_data_rel <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP, SEX) %>%
  summarise(
    mean_bw_rel = mean(bw_rel, na.rm = TRUE),
    sem_bw_rel = sd(bw_rel, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

#plot 3 % BW gain in each STATUS collapsed by DIET_FORMULA----

ggplot() +
  # Individual points, smaller and more transparent
  geom_point(data = BW_data, 
             aes(x = STATUS, y = bw_rel, color = SEX, group = ID),
             size = 1.2, alpha = 0.3) +
  geom_line(data = BW_data, 
            aes(x = STATUS, y = bw_rel, group = ID),
            alpha = 0.15) +
  
  # Mean ± SEM as shaded ribbon
  geom_ribbon(data = summary_data_rel,
              aes(x = STATUS, ymin = mean_bw_rel - sem_bw_rel, ymax = mean_bw_rel + sem_bw_rel,
                  group = interaction(GROUP, SEX), fill = SEX),
              alpha = 0.2) +
  geom_line(data = summary_data_rel,
            aes(x = STATUS, y = mean_bw_rel, group = interaction(GROUP, SEX), color = SEX),
            size = 1) +
  
  # Facet by STRAIN (rows) × GROUP (columns)
  facet_grid(~ STRAIN*GROUP, scales = "free_x") +
  
  # Colors
  scale_color_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  scale_fill_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  
  # Labels and theme
  labs(
    x = "Status",
    y = "Body weight gain (%)",
    color = "Sex",
    fill = "Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  format.plot



