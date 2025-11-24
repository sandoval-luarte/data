
#deltafat/deltalean
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
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_Z",
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_M",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    ),
    STATUS = case_when(
      Date == as.Date("2025-04-14") ~ "baseline",
      Date == as.Date("2025-04-21") ~ "end_drug"
    )
  ) %>%
  select(ID, Fat, Lean, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    delta_fat = Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"],
    delta_lean = Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"],
    delta_adiposity_index = delta_fat / delta_lean,
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M")))

# -----------------------------
# --- ART model -----
# -----------------------------
art_model <- art(delta_adiposity_index ~ DRUG, data = echoMRI_data_43)

# ANOVA table
anova(art_model)

# -----------------------------
# --- Post-hoc pairwise comparisons (ART) -----
# -----------------------------
# Use art.con() for pairwise comparisons
pairwise_drug <- art.con(art_model, "DRUG")
pairwise_drug

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
  mutate(delta_adiposity_index = ((Fat_n2 - Fat_n1)/(Lean_n2 - Lean_n1)),
         DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  filter(ID != 8076) #this animal has the same Lean mass at the and and at the start so goes to infinite 

# Shapiro-Wilk test for all data
shapiro.test(BW_wide$delta_adiposity_index) #so this data is normal so we can use a ANOVA then

# One-way ANOVA
anova_result_47 <- aov(delta_adiposity_index ~ DRUG, data = BW_wide)
summary(anova_result_47)

# Pairwise t-test
pairwise_result_47 <- pairwise.t.test(
  x = BW_wide$delta_adiposity_index,
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
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki")) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index, GROUP, DRUG, SEX, STRAIN, STATUS)

# Calculate delta values
echoMRI_delta <- echoMRI_data %>%
  group_by(ID, DRUG, STRAIN, SEX, GROUP) %>%
  summarise(
    delta_adiposity_index = ((Fat[STATUS == "BW regain"][1] - Fat[STATUS == "BW maintenance"][1])/(Lean[STATUS == "BW regain"][1] - Lean[STATUS == "BW maintenance"][1])),
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")))

# Shapiro-Wilk test for all data
shapiro.test(echoMRI_delta$delta_adiposity_index) #so this data is normal so we can use a ANOVA then

# ART + post-hoc
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
  bind_rows() %>%
  left_join(echoMRI_delta %>% select(STRAIN, SEX, GROUP) %>% distinct(), 
            by = c("STRAIN", "SEX"))

# Plotting positions for p-values
y_positions <- echoMRI_delta %>%
  group_by(STRAIN, SEX, GROUP) %>%
  summarise(y_pos = max(delta_adiposity_index, na.rm = TRUE) * 1.05, .groups = "drop")

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
y_max_5D <- max(echoMRI_data_43$delta_adiposity_index, na.rm = TRUE) * 1.1

plot_43_5D <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "ΔFat/ΔLean") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("vehicle", "RTIOXA_43_M"), c("vehicle", "RTIOXA_43_Z")),
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D

plot_47_5D <- ggplot(BW_wide, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +  # jitter to see overlapping points
  theme_minimal() +
  labs(x = NULL, y = "ΔFat/ΔLean") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("vehicle", "RTIOXA_47")),
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("vehicle" = "gray70", "RTIOXA_47" = "#8DA0CB")) +
  theme(legend.position = "none")
plot_47_5D 


plot_47_5W <- ggplot(echoMRI_delta, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(x = DRUG, y = delta_adiposity_index), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "ΔFat/ΔLean") +
  facet_grid(GROUP ~ STRAIN + SEX, scales = "free_y") +  # <-- free y-axis per facet
  scale_fill_manual(values = c("vehicle" = "gray70", "RTIOXA_47" = "#8DA0CB")) +
  stat_pvalue_manual(pairwise_df, label = "p.signif", inherit.aes = FALSE)
plot_47_5W 

# -----------------------------
# --- Combine plots -----------
# -----------------------------
combined_plot <- ((plot_43_5D | plot_47_5D) / plot_47_5W) +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A")+
  theme(legend.position = "none")

combined_plot

