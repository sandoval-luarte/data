
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
library(car)
library(patchwork)


#format plot
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

# RTIOXA-43 × 5 days -----

# first we used the adiposity index as a criteria to assignt the groups

echoMRI_data_43_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_Z",
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_M",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M")))

plot_43_5D_assignation <- ggplot(echoMRI_data_43_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_assignation

#pairwise comparison to evaluate if the adiposity index is the same among the groups at the baseline

#to check if variances differ between groups,
echoMRI_data_43_assignation <- echoMRI_data_43_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_43_assignation) #p = 0.6489 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation$adiposity_index,
  g = echoMRI_data_43_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between any of the three groups

#Now after we injected RTIOXA-43 × 5 days

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
  select(ID, Fat, Lean,adiposity_index,Weight, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    delta_adiposity_index = adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"],
    delta_lean = Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"],
    delta_fat = Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"],
    delta_bw = Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M")))

#plot bw changes

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ body weight") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_bw

#plot adiposity index changes

plot_43_5D_ai <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ adiposity index") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_ai

#plot lean mass changes

plot_43_5D_lean <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ lean mass") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_lean

#plot fat mass changes

plot_43_5D_fat <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ fat mass") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_fat

# combined plot for RTIOXA-43 x 5 days

combined_plot <- (plot_43_5D_bw | plot_43_5D_ai) /
  (plot_43_5D_fat | plot_43_5D_lean) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

combined_plot

#data analysis
 #1 does rtioxa-43 decreases bw after 5 days of injection?
#to check if variances differ between groups,
echoMRI_data_43 <- echoMRI_data_43 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_43) #p = 0.75 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_43) #p = 0.63 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_43) #p = 0.18 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_43) #p = 0.4 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_43))
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_43))
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_43))
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_43))

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, 
# or adiposity after 5 days. Probably because data are underpowered (n=3 per group).

# RTIOXA-47 × 5 days -----

echoMRI_data_47_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 12) %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8075, 8077, 8078) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079) ~ "vehicle"
    )
  ) %>%
  filter(n_measurement == 1) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_47")))


plot_47_5D_assignation <- ggplot(echoMRI_data_47_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange"))+
  theme(legend.position = "none")

plot_47_5D_assignation

#pairwise comparison to evaluate if the adiposity index is the same between groups at the baseline

#to check if variances differ between groups,
echoMRI_data_47_assignation <- echoMRI_data_47_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_47_assignation) #p = 0.9> 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_47_assignation$adiposity_index,
  g = echoMRI_data_47_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups


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

