# This script aims to explore changes in FI in NZO and c57 after different stages of feeding:
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


#FI CSV data import RTIOXA 47 ----
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  group_by(ID) %>%
  arrange(DATE) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         date_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
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
    ) ) %>% 
filter(STRAIN == "NZO/HlLtJ") %>% 
  filter(ID==3708)

n_distinct(FI_data$ID) #9 NZO in total and 4 C57 in total

highlight_data <- FI_data %>%
  filter(STATUS == "peak obesity")

highlight_data_2 <- FI_data %>%
  filter(STATUS == "BW loss")

highlight_data_3 <- FI_data %>%
  filter(STATUS == "BW maintenance")

highlight_data_4 <- FI_data %>%
  filter(STATUS == "BW regain")

#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14))


#plot 1 description of FI over time per ID----
plot <- FI_data %>% 
  ggplot(aes(day_rel, FI_cum, group = ID)) +
  geom_point() +
  geom_line() +
  #facet_wrap(~STRAIN*GROUP) +
 # geom_smooth() +
  labs(
    x = "Days relative to baseline",
    y = "FI (cumulative kcal)"
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

