# This script aims to explore changes in cumulative food intake NZO and c57 within different stages of feeding:
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
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
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
    )) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) & DATE == as.Date("2025-04-04") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725) & DATE == as.Date("2025-04-18") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-04-25") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
      STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-06-12") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3715, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-21") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-25") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID == 3726 & DATE == as.Date("2025-06-27") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719, 3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
      STRAIN == "C57BL6/J"  & ID %in% c(7872, 7874, 7877, 7878) & DATE == as.Date("2025-03-05") ~ "peak obesity",
      STRAIN == "C57BL6/J"  & ID %in% c(7861, 7863, 7865, 7866) & DATE == as.Date("2025-03-24") ~ "peak obesity",
      STRAIN == "C57BL6/J"  & ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870,
                                        7871, 7873, 7875, 7876, 7879, 7880,
                                        7881, 7882, 7883) & DATE == as.Date("2025-03-24") ~ "peak obesity",
      STRAIN == "C57BL6/J"  & ID %in% c( 7873, 7874, 7875, 7876, 7877, 7878,
                                        7879, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-04-09") ~ "BW loss",
      STRAIN == "C57BL6/J"  & ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870,
                                        7871) & DATE == as.Date("2025-04-30") ~ "BW loss",
      STRAIN == "C57BL6/J" & ID == 7872 & DATE == as.Date("2025-04-10") ~ "BW loss",
      STRAIN == "C57BL6/J"  & ID %in% c(7861, 7863, 7865, 7866) & DATE == as.Date("2025-05-02") ~ "BW loss",
      STRAIN == "C57BL6/J"  & ID %in% c(7860, 7868, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-08-15") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7861, 7862, 7863, 7864) & DATE == as.Date("2025-08-01") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7865, 7878) & DATE == as.Date("2025-08-02") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7866, 7877) & DATE == as.Date("2025-08-10") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7867, 7872, 7873, 7874, 7875, 7876) & DATE == as.Date("2025-08-06") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7869, 7870, 7871, 7879) & DATE == as.Date("2025-08-08") ~ "BW maintenance",
      STRAIN == "C57BL6/J"  & ID %in% c(7860, 7868, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-09-10") ~ "BW regain",
      STRAIN == "C57BL6/J"  & ID %in% c(7861, 7862, 7863, 7864, 7865, 7878) & DATE == as.Date("2025-08-27") ~ "BW regain",
      STRAIN == "C57BL6/J"  & ID %in% c(7866, 7877, 7879) & DATE == as.Date("2025-09-05") ~ "BW regain",
      STRAIN == "C57BL6/J"  & ID %in% c(7867, 7872, 7873, 7874, 7875, 7876) & DATE == as.Date("2025-09-01") ~ "BW regain",
      STRAIN == "C57BL6/J"  & ID %in% c(7869, 7870, 7871) & DATE == as.Date("2025-09-03") ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()


FI_data %>% 
group_by(STRAIN,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this we have 22 NZO and 46 C57 per STATUS

# Summarize cumulative FI per ID and STATUS ----
FI_stage_summary <- FI_data %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(FI_cum_end = max(FI_cum, na.rm = TRUE), .groups = "drop") %>%
  # Reshape into wide format: one row per ID, columns = each STATUS
  pivot_wider(
    names_from = STATUS,
    values_from = FI_cum_end
  ) %>%
  # Calculate kcal consumed between stages
  mutate(
    kcal_baseline_to_peak = `peak obesity` - baseline,
    kcal_peak_to_loss = `BW loss` - `peak obesity`,
    kcal_loss_to_maint = `BW maintenance` - `BW loss`,
    kcal_maint_to_regain = `BW regain` - `BW maintenance`
  )

FI_stage_summary #now we can plot this


# Convert to long format for plotting ----
FI_stage_long <- FI_stage_summary %>%
  select(ID, STRAIN, DRUG, GROUP,
         kcal_baseline_to_peak,
         kcal_peak_to_loss,
         kcal_loss_to_maint,
         kcal_maint_to_regain) %>%
  pivot_longer(
    cols = starts_with("kcal_"),
    names_to = "transition",
    values_to = "kcal"
  ) %>%
  mutate(
    transition = factor(
      transition,
      levels = c("kcal_baseline_to_peak",
                 "kcal_peak_to_loss",
                 "kcal_loss_to_maint",
                 "kcal_maint_to_regain"),
      labels = c("Baseline → Peak Obesity",
                 "Peak Obesity → BW Loss",
                 "BW Loss → BW Maintenance",
                 "BW Maintenance → BW Regain")
    )
  )

# Plot ---- 
ggplot(FI_stage_long, aes(x = transition, y = kcal, fill = GROUP)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.3
  ) +
  facet_wrap(~ STRAIN*GROUP, scales = "free_y") +
  labs(
    x = "Transition",
    y = "Calories Consumed (kcal)",
    fill = "Feeding Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) 
#this is weird because it seems C57 ate more from baseline to peak of obesity and that is
#obvious because they were exposed in mean 100 days more than NZO to the diets during that phase
#we should evaluate a better approximation to standardize for time of exposure to the diets


#OK I will standardize per BW of each ID within each STATUS


# --- Step 1: Import FI and BW data ---
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO females
  filter(DIET_FORMULA.x !="2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet") %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  group_by(ID) %>%
  arrange(DATE) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(
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
    )) %>%
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         day_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         STATUS = case_when(
           day_rel == 0 ~ "baseline", 
           STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708, 3709, 3710, 3711, 3713, 3714) & DATE == as.Date("2025-04-04") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725) & DATE == as.Date("2025-04-18") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-04-25") ~ "BW loss",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-06-11") ~ "BW maintenance", 
           STRAIN == "NZO/HlLtJ" & ID == 3710 & DATE == as.Date("2025-06-12") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3711) & DATE == as.Date("2025-06-13") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3715, 3716) & DATE == as.Date("2025-06-18") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719) & DATE == as.Date("2025-06-20") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3720, 3721, 3722) & DATE == as.Date("2025-06-21") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-06-25") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3727, 3728, 3729) & DATE == as.Date("2025-06-26") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID == 3726 & DATE == as.Date("2025-06-27") ~ "BW maintenance",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3706, 3707, 3708) & DATE == as.Date("2025-07-07") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3709, 3710, 3711) & DATE == as.Date("2025-07-09") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3713, 3714, 3716) & DATE == as.Date("2025-07-14") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3717, 3718, 3719, 3720, 3721, 3722) & DATE == as.Date("2025-07-16") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3723, 3724, 3725) & DATE == as.Date("2025-07-21") ~ "BW regain",
           STRAIN == "NZO/HlLtJ" & ID %in% c(3726, 3727, 3728, 3729) & DATE == as.Date("2025-07-22") ~ "BW regain",
           STRAIN == "C57BL6/J"  & ID %in% c(7872, 7874, 7877, 7878) & DATE == as.Date("2025-03-05") ~ "peak obesity",
           STRAIN == "C57BL6/J"  & ID %in% c(7861, 7863, 7865, 7866) & DATE == as.Date("2025-03-24") ~ "peak obesity",
           STRAIN == "C57BL6/J"  & ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870,
                                             7871, 7873, 7875, 7876, 7879, 7880,
                                             7881, 7882, 7883) & DATE == as.Date("2025-03-24") ~ "peak obesity",
           STRAIN == "C57BL6/J"  & ID %in% c( 7873, 7874, 7875, 7876, 7877, 7878,
                                              7879, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-04-09") ~ "BW loss",
           STRAIN == "C57BL6/J"  & ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870,
                                             7871) & DATE == as.Date("2025-04-30") ~ "BW loss",
           STRAIN == "C57BL6/J" & ID == 7872 & DATE == as.Date("2025-04-10") ~ "BW loss",
           STRAIN == "C57BL6/J"  & ID %in% c(7861, 7863, 7865, 7866) & DATE == as.Date("2025-05-02") ~ "BW loss",
           STRAIN == "C57BL6/J"  & ID %in% c(7860, 7868, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-08-15") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7861, 7862, 7863, 7864) & DATE == as.Date("2025-08-01") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7865, 7878) & DATE == as.Date("2025-08-02") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7866, 7877) & DATE == as.Date("2025-08-10") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7867, 7872, 7873, 7874, 7875, 7876) & DATE == as.Date("2025-08-06") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7869, 7870, 7871, 7879) & DATE == as.Date("2025-08-08") ~ "BW maintenance",
           STRAIN == "C57BL6/J"  & ID %in% c(7860, 7868, 7880, 7881, 7882, 7883) & DATE == as.Date("2025-09-10") ~ "BW regain",
           STRAIN == "C57BL6/J"  & ID %in% c(7861, 7862, 7863, 7864, 7865, 7878) & DATE == as.Date("2025-08-27") ~ "BW regain",
           STRAIN == "C57BL6/J"  & ID %in% c(7866, 7877, 7879) & DATE == as.Date("2025-09-05") ~ "BW regain",
           STRAIN == "C57BL6/J"  & ID %in% c(7867, 7872, 7873, 7874, 7875, 7876) & DATE == as.Date("2025-09-01") ~ "BW regain",
           STRAIN == "C57BL6/J"  & ID %in% c(7869, 7870, 7871) & DATE == as.Date("2025-09-03") ~ "BW regain",
           TRUE ~ NA_character_
         )) %>% 
  filter(!is.na(STATUS)) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  ungroup()

BW_data <- read_csv("../data/BW_data.csv") %>%
  select(ID, STATUS, BW) %>% 
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))


# --- Step 2: Summarize FI per ID & STATUS ---
FI_stage <- FI_data %>%
  group_by(ID, STRAIN, GROUP, DRUG, STATUS) %>%
  summarise(
    FI_cum_end = max(FI_cum, na.rm = TRUE),
    day_end = max(day_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(BW_data, by = c("ID", "STATUS"))

# Create transition columns
df_transitions <- FI_stage %>%
  arrange(ID, day_end) %>%
  group_by(ID) %>%
  mutate(
    STATUS_TRANSITION = paste(STATUS, "->", lead(STATUS)),
    BW_change = lead(BW) - BW,
    delta_days = lead(day_end) - day_end,
    delta_FI_cum = lead(FI_cum_end) - FI_cum_end,
  ) %>%
  # Remove last row (no transition after final STATUS)
  filter(!is.na(BW_change)) %>%
  ungroup()

df_transitions %>%
  select(ID, STATUS_TRANSITION, BW_change, delta_days)

df_transitions <- df_transitions %>% 
  select(ID, STRAIN, GROUP, DRUG, STATUS_TRANSITION, BW_change, delta_days, delta_FI_cum)

df_transitions <- df_transitions %>%
  mutate(
    delta_days_num = as.numeric(gsub(" days", "", delta_days)),
    FI_kcal_per_day = delta_FI_cum / delta_days_num,
    FI_kcal_per_BWchange = ifelse(abs(BW_change) > 0.001, delta_FI_cum / BW_change, NA),
    energy_efficiency = ((BW_change /delta_FI_cum)/delta_days_num),
    DRUG = factor(DRUG)) 

# Define order for the transitions
status_order <- c(
  "baseline -> peak obesity",
  "peak obesity -> BW loss",
  "BW loss -> BW maintenance",
  "BW maintenance -> BW regain"
)

# Prepare data
df_transitions <- df_transitions %>%
  mutate(
    delta_days_num = as.numeric(gsub(" days", "",delta_days)),
    FI_kcal_per_day = delta_FI_cum / delta_days_num,
    facet_group = interaction(STRAIN,GROUP),
    STATUS_TRANSITION = factor(STATUS_TRANSITION, levels = status_order)) 

# Bar plot with fixed y-axis across all facets
ggplot(df_transitions, aes(x = STATUS_TRANSITION, y = FI_kcal_per_day, fill = STRAIN)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, 
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~ facet_group, scales = "fixed") +   # <-- fixed y-axis
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    x = "Status transition",
    y = "Food intake (kcal/day)",
    fill = "Strain",
  )

# Bar plot with FI_kcal_per_BWchange and free y-axis
ggplot(df_transitions, aes(x = STATUS_TRANSITION, y = FI_kcal_per_BWchange, fill = STRAIN)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", fun.data = mean_se, 
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~ facet_group, scales = "free_y") +   # free y-axis
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    x = "Status transition",
    y = "Food intake (kcal / g BW change)",
    fill = "Strain"
  )


# Bar plot with energy efficiency (g BW gained per kcal)
ggplot(df_transitions, aes(x = STATUS_TRANSITION, y = energy_efficiency, fill= STRAIN)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~ facet_group, scales = "fixed") +
  geom_hline(data = data.frame(yintercept = 0), 
             aes(yintercept = yintercept), 
             color = "black", linetype = "dashed") +  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    x = "Status transition",
    y = "Energy efficiency (g BW change per kcal consumed)",
    fill = "Strain"
  )


