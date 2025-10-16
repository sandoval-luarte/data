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

# -----------------------------
# --- RTIOXA-43 × 5 days -----
# -----------------------------
echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>%
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_donated",
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_Medchem",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    ),
    STATUS = case_when(
      Date == as.Date("2025-04-14") ~ "baseline",
      Date == as.Date("2025-04-21") ~ "end_drug"
    )
  ) %>%
  select(ID, Fat, Lean, SEX, STRAIN, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    delta_fat = Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"],
    delta_lean = Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"],
    delta_adiposity_index = delta_fat / delta_lean,
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_donated", "RTIOXA_43_Medchem")))

# One-way ANOVA
anova_result <- aov(delta_fat ~ DRUG, data = echoMRI_data_43)
summary(anova_result)

# Post-hoc pairwise t-tests (Bonferroni)
pairwise_result <- pairwise.t.test(
  x = echoMRI_data_43$delta_fat,
  g = echoMRI_data_43$DRUG,
  p.adjust.method = "bonferroni"
)
pairwise_result

# -----------------------------
# --- RTIOXA-47 × 5 days -----
# -----------------------------
BW_data_47 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 12) %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079) ~ "vehicle"
    ),
    STRAIN = "C57BL6/J"
  ) %>%
  ungroup()

# Select baseline and end
BW_compare <- BW_data_47 %>%
  filter(n_measurement %in% c(1, 2)) %>%
  select(ID, DRUG, STRAIN, SEX, n_measurement, Fat, Lean)

# Pivot wider
BW_wide <- BW_compare %>%
  pivot_wider(
    id_cols = c(ID, DRUG, STRAIN, SEX),
    names_from = n_measurement,
    values_from = c(Fat, Lean),
    names_prefix = "n"
  ) %>%
  mutate(delta_fat = Fat_n2 - Fat_n1,
         DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")))

# One-way ANOVA
anova_result_47 <- aov(delta_fat ~ DRUG, data = BW_wide)
summary(anova_result_47)

# Pairwise t-test
pairwise_result_47 <- pairwise.t.test(
  x = BW_wide$delta_fat,
  g = BW_wide$DRUG,
  p.adjust.method = "bonferroni"
)
pairwise_result_47

# -----------------------------
# --- RTIOXA-47 × 5 weeks -----
# -----------------------------
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 1 & COHORT < 6) %>%
  filter(!ID %in% c(3712, 3715)) %>%
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

# Calculate delta values
echoMRI_delta <- echoMRI_data %>%
  group_by(ID, DRUG, DIET_FORMULA, STRAIN, SEX, GROUP) %>%
  summarise(
    delta_fat = Fat[STATUS == "BW regain"][1] - Fat[STATUS == "BW maintenance"][1],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")))

# ART + post-hoc
art_results <- echoMRI_delta %>%
  group_by(STRAIN, SEX) %>%
  group_map(~ {
    df <- .x
    keys <- .y
    if(length(unique(df$DRUG)) == 2) {
      art_model <- art(delta_fat ~ DRUG, data = df)
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
  bind_rows() %>%
  left_join(echoMRI_delta %>% select(STRAIN, SEX, GROUP) %>% distinct(), 
            by = c("STRAIN", "SEX"))

# Plotting positions for p-values
y_positions <- echoMRI_delta %>%
  group_by(STRAIN, SEX, GROUP) %>%
  summarise(y_pos = max(delta_fat, na.rm = TRUE) * 1.05, .groups = "drop")

pairwise_df <- art_results %>%
  left_join(y_positions, by = c("STRAIN", "SEX", "GROUP")) %>%
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

# -----------------------------
# --- PLOTS -------------------
# -----------------------------
# Y-axis limit for 5-day plots
y_max_5D <- max(c(echoMRI_data_43$delta_fat, BW_wide$delta_fat), na.rm = TRUE) * 1.1

plot_43_5D <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2) +
  scale_y_continuous(limits = c(0, y_max_5D)) +
  theme_minimal() +
  labs(title = "RTIOXA-43 (5 days)", x = NULL, y = "ΔFat") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("vehicle", "RTIOXA_43_Medchem"), c("vehicle", "RTIOXA_43_donated")),
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_donated" = "#66C2A5", "RTIOXA_43_Medchem" = "darkgreen"))

plot_47_5D <- ggplot(BW_wide, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2) +
  scale_y_continuous(limits = c(0, y_max_5D)) +
  theme_minimal() +
  labs(title = "RTIOXA-47 (5 days)", x = NULL, y = "ΔFat") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("vehicle", "RTIOXA_47")),
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange"))

plot_47_5W <- ggplot(echoMRI_delta, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(x = DRUG, y = delta_fat), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(title = "RTIOXA-47 (5 weeks)", x = NULL, y = "ΔFat") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  facet_grid(GROUP ~ STRAIN + SEX) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  stat_pvalue_manual(pairwise_df, label = "p.signif", inherit.aes = FALSE)

# -----------------------------
# --- Combine plots -----------
# -----------------------------
combined_plot <- ggarrange(
  plot_43_5D,
  plot_47_5D,
  plot_47_5W,
  ncol = 3,
  nrow = 1,
  labels = c("A", "B", "C"),
  align = "v"
)

combined_plot
