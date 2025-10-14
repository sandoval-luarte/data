library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(ARTool)

# C57BL6J RTIOXA-43 × 5 days in a row ----

echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% # C57BL6J RTIOXA-43 × 5 days in a row 
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_donated",
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_Medchem",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    ),
    STATUS = case_when(
      Date == as.Date("2025-04-14") ~ "baseline",
      Date == as.Date("2025-04-21") ~ "end drug"
    )
  ) %>%
  select(ID, Fat, Lean, Date, SEX, STRAIN, adiposity_index, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  mutate(
    delta_fat  = Fat[STATUS == "end drug"] - Fat[STATUS == "baseline"],
    delta_lean = Lean[STATUS == "end drug"] - Lean[STATUS == "baseline"],
    delta_adiposity_index = delta_fat / delta_lean
  ) %>%
  ungroup() %>%
  # ensure one row per ID for plotting
  distinct(ID, DRUG, delta_adiposity_index) %>%
  # explicitly set factor order
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_donated", "RTIOXA_43_Medchem")))

# One-way ANOVA
anova_result <- aov(delta_adiposity_index ~ DRUG, 
                    data = distinct(echoMRI_data_43, ID, DRUG, delta_adiposity_index))
summary(anova_result)

# Post-hoc pairwise t-tests with Bonferroni correction
pairwise_result <- pairwise.t.test(
  x = distinct(echoMRI_data_43, ID, DRUG, delta_adiposity_index)$delta_adiposity_index,
  g = distinct(echoMRI_data_43, ID, DRUG, delta_adiposity_index)$DRUG,
  p.adjust.method = "bonferroni"
)
pairwise_result


ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(x = NULL, y = "ΔFat / ΔLean") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(
      c("vehicle", "RTIOXA_43_Medchem"),
      c("vehicle", "RTIOXA_43_donated")
    ),
    label = "p.signif"
  ) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_Medchem" = "darkgreen"
  ))

# C57BL6J and NZO RTIOXA-47 × 5 weeks ----
library(dplyr)
library(readr)
library(ARTool)
library(artTools)
library(emmeans)
library(ggplot2)
library(ggpubr)

# --- 1. Load and preprocess data ---
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 1 & COHORT < 6) %>% # NZO and C57BL6J females
  filter(!ID %in% c(3712, 3715)) %>%  # remove animals that died
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876,
                7879, 7880, 7881, 7882, 7883) ~ "ad lib",
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
    STATUS = case_when(
      STRAIN == "NZO/HlLtJ" & Date == as.Date("2025-05-27") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                                                  "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-06-05") ~ "BW maintenance",
      STRAIN == "C57BL6/J" & Date %in% as.Date(c("2025-09-11", "2025-09-10","2025-09-05","2025-09-04",
                                                 "2025-09-02","2025-09-01","2025-08-28","2025-08-27")) ~ "BW regain",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(STATUS)) %>%
  select(ID, Date, Fat, Lean, Weight, adiposity_index, GROUP, DRUG, SEX, STRAIN, DIET_FORMULA, STATUS)

# --- 2. Calculate delta values ---
echoMRI_delta <- echoMRI_data %>%
  group_by(ID, DRUG, DIET_FORMULA, STRAIN, SEX) %>%
  mutate(
    delta_fat  = Fat[STATUS == "BW regain"][1] - Fat[STATUS == "BW maintenance"][1],
    delta_lean = Lean[STATUS == "BW regain"][1] - Lean[STATUS == "BW maintenance"][1],
    delta_adiposity_index = delta_fat / delta_lean
  ) %>%
  ungroup() %>%
  distinct(ID, DRUG, delta_adiposity_index, STRAIN, SEX, DIET_FORMULA) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    SEX = factor(SEX),
    STRAIN = factor(STRAIN)
  ) %>%
  select(ID, DRUG, STRAIN, SEX, delta_adiposity_index)

# --- 3. Run ART + post-hoc per STRAIN × SEX ---
art_results <- echoMRI_delta %>%
  group_by(STRAIN, SEX) %>%
  group_map(~ {
    df <- .x
    keys <- .y
    if(length(unique(df$DRUG)) == 2) {
      art_model <- art(delta_adiposity_index ~ DRUG, data = df)
      lm_model <- artlm(art_model, "DRUG")
      emm <- emmeans(lm_model, ~ DRUG)
      pairwise <- pairs(emm, adjust = "bonferroni") %>% as.data.frame()
      pairwise$STRAIN <- keys$STRAIN
      pairwise$SEX <- keys$SEX
      return(pairwise)
    } else {
      return(NULL)
    }
  }) %>%
  bind_rows()


# --- 4. Prepare plotting positions ---
y_positions <- echoMRI_delta %>%
  group_by(STRAIN, SEX) %>%
  summarise(y_pos = max(delta_adiposity_index, na.rm = TRUE) * 1.05, .groups = "drop")

pairwise_df <- art_results %>%
  left_join(y_positions, by = c("STRAIN", "SEX")) %>%
  rename(y.position = y_pos) %>%
  mutate(
    group1 = "vehicle",
    group2 = "RTIOXA_47",
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

ggplot(echoMRI_delta, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(x = DRUG, y = delta_adiposity_index), 
             data = echoMRI_delta, 
             alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "ΔFat / ΔLean") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  facet_wrap(~STRAIN*SEX) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  stat_pvalue_manual(pairwise_df, label = "p.signif",
                     inherit.aes = FALSE)  # important


library(dplyr)
library(ggplot2)
library(ggpubr)

# --- 1. Add Experiment column ---
echoMRI_data_43_plot <- echoMRI_data_43 %>%
  mutate(
    Experiment = "RTIOXA-43 (C57BL6J, 5 days)",
    STRAIN = "C57BL6J",
    SEX = "M"  # set SEX since 43 dataset is all male
  )

echoMRI_delta_plot <- echoMRI_delta %>%
  mutate(Experiment = "RTIOXA-47 (C57BL6J & NZO, 5 weeks)")

# --- 2. Combine datasets ---
combined_data <- bind_rows(
  echoMRI_data_43_plot %>% select(ID, DRUG, STRAIN, SEX, delta_adiposity_index, Experiment),
  echoMRI_delta_plot %>% select(ID, DRUG, STRAIN, SEX, delta_adiposity_index, Experiment)
)

# --- 3. Optional: Combine pairwise results for plotting (simplified example) ---
# For RTIOXA-43 (vehicle vs each drug)
pairwise_43 <- data.frame(
  Experiment = "RTIOXA-43 (C57BL6J, 5 days)",
  STRAIN = "C57BL6J",
  SEX = "M",
  group1 = c("vehicle", "vehicle"),
  group2 = c("RTIOXA_43_donated", "RTIOXA_43_Medchem"),
  y.position = max(echoMRI_data_43$delta_adiposity_index) * 1.05,
  p.signif = c("*", "*") # or use actual t-test results
)

# For RTIOXA-47, reuse pairwise_df from your code
pairwise_47 <- pairwise_df %>%
  mutate(Experiment = "RTIOXA-47 (C57BL6J & NZO, 5 weeks)")

pairwise_combined <- bind_rows(pairwise_43, pairwise_47)

# --- 4. Plot combined ---
ggplot(combined_data, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(x = DRUG, y = delta_adiposity_index), 
             alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "ΔFat / ΔLean") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  facet_wrap(~Experiment + STRAIN + SEX, scales = "free_x") +
  scale_fill_manual(values = c(
    "vehicle" = "white", 
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_Medchem" = "darkgreen",
    "RTIOXA_47" = "orange"
  )) +
  stat_pvalue_manual(pairwise_combined, label = "p.signif", inherit.aes = FALSE)

