#We aim to explore changes in body composition (bw, fat mass, lean mass and adiposity index)
#in 3 adult males C57BL6J after 5 days of IP injection of RTIOXA-47

#libraries

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

# group assignation -----

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

# RTIOXA-47 × 5 days -----

echoMRI_data_47 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 12) %>%
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079) ~ "vehicle"
    )
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
         STATUS = case_when(
           n_measurement == 1 ~ "baseline",
          n_measurement == 2 ~ "end_drug"
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
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")))

#plot bw changes

plot_47_5D_bw <- ggplot(echoMRI_data_47, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ body weight") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_bw

#plot adiposity index changes

plot_47_5D_ai <- ggplot(echoMRI_data_47, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ adiposity index") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_ai

#plot lean mass changes

plot_47_5D_lean <- ggplot(echoMRI_data_47, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ lean mass") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_lean

#plot fat mass changes

plot_47_5D_fat <- ggplot(echoMRI_data_47, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Δ fat mass") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_fat

# combined plot for RTIOXA-47 x 5 days

combined_plot <- (plot_47_5D_bw | plot_47_5D_ai) /
  (plot_47_5D_fat | plot_47_5D_lean) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

combined_plot

#data analysis
#1 does rtioxa-47 decreases bw after 5 days of injection?
#to check if variances differ between groups,
echoMRI_data_47 <- echoMRI_data_47 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_47) #p = 0.99 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_47) #p = 0.32 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_47) #p = 0.91 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_47) #p = 0.85 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_47))
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_47))
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_47))
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_47))

#conclusion> There is no indication that RTIOXA-47 decreases body weight, fat mass, 
# or adiposity after 5 days. Probably because data are underpowered (n=3 per group).
