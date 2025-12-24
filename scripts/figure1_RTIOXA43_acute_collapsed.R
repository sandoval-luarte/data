
#We aim to explore changes in body composition (bw, fat mass, lean mass and adiposity index)
#in 10 adult males C57BL6J after 5 days of IP injection of RTIOXA-43

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
library(car)
library(patchwork)
library(ggrepel)


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

# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(9, 14)) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8096, 8099, 8102, 41006, 41007, 41008, 41015, 41001, 41018, 41021) ~ "RTIOXA_43", 
      ID %in% c(8095, 8100, 8103, 41003, 41005, 41012, 41010, 41016, 41014, 41020) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

plot_43_5D_assignation <- ggplot(echoMRI_data_43_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_assignation

#pairwise comparison to evaluate if the adiposity index is the same among the groups at the baseline

#to check if variances differ between groups,
echoMRI_data_43_assignation <- echoMRI_data_43_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_43_assignation) #p = 0.47 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation$adiposity_index,
  g = echoMRI_data_43_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups

# RTIOXA-43 × 5 days ----

echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(9, 14)) %>% 
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21", "2025-12-18","2025-12-23"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102, 41006, 41007,
                41008, 41015, 41001, 41018, 41021) ~ "RTIOXA_43", 
      ID %in% c(8095, 8100, 8103, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle"),
    STATUS = case_when(
      COHORT == 9 & Date == as.Date("2025-04-14") ~ "baseline",
      COHORT == 9 & Date == as.Date("2025-04-21") ~ "end_drug",
      COHORT == 14 & Date == as.Date("2025-12-18") ~ "baseline",
      COHORT == 14 & Date == as.Date("2025-12-23") ~ "end_drug")
  ) %>%
  select(ID, Fat, Lean,adiposity_index,Weight, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% of body weight change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_bw

#plot adiposity index changes

plot_43_5D_ai <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% adiposity index change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_ai

#plot lean mass changes

plot_43_5D_lean <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% lean mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_lean

#plot fat mass changes

plot_43_5D_fat <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% fat mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
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

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_43) #p = 0.72 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_43) #p = 0.65 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_43) #p = 0.14 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_43) #p = 0.42 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_43))  #p = 0.4
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_43))  #p = 0.18
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_43))  #p = 0.59
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_43))  #p = 0.16

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, 
# or adiposity after 5 days. 
