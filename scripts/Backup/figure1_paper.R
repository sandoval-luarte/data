
#We aim to explore changes in body composition (bw, fat mass, lean mass and adiposity index)
#in 3 adult males C57BL6J after 5 days of IP injection of RTIOXA-43

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
## Differentiated analysis for RTIOXA 43 obtained from Yanan and from medchem----

# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8096,8099,  8102) ~ "RTIOXA_43_Z", 
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_M",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M"))) %>% 
  drop_na()

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

# RTIOXA-43 × 5 days ----
## Differentiated analysis for RTIOXA 43 obtained from Yanan and from medchem----


echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>%
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096,8099,   8102) ~ "RTIOXA_43_Z", #8099 eliminated for now
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
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% of body weight change") +
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
  labs(x = NULL, y = "% adiposity index change") +
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
  labs(x = NULL, y = "% lean mass change") +
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
  labs(x = NULL, y = "% fat mass change") +
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
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_43) #p = 0.85 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_43) #p = 0.18 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_43) #p = 0.46 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_43))  #p = 0.61
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_43))  #p = 0.15
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_43))  #p = 0.15
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_43))  #p = 0.22

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, 
# or adiposity after 5 days. Probably because data are underpowered (n=3 per group).

# group assignation -----
## Collapsed RTIOXAs 43 ----

# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation2 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION2 = case_when(
      ID %in% c(8096, 8102, 8097, 8098, 8101, 8099) ~ "RTIOXA_43",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION2 = factor(DRUG_ASSIGNATION2, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

plot_43_5D_assignation2 <- ggplot(echoMRI_data_43_assignation2, aes(x = DRUG_ASSIGNATION2, y = adiposity_index, fill = DRUG_ASSIGNATION2)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_assignation2

#pairwise comparison to evaluate if the adiposity index is the same among the groups at the baseline

#to check if variances differ between groups,
echoMRI_data_43_assignation2 <- echoMRI_data_43_assignation2 %>%
  mutate(DRUG_ASSIGNATION2 = factor(DRUG_ASSIGNATION2))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION2, data = echoMRI_data_43_assignation2) #p = 0.29 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation2$adiposity_index,
  g = echoMRI_data_43_assignation2$DRUG_ASSIGNATION2,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups

# RTIOXA-43 × 5 days (collapsed) ---- 

echoMRI_data_432 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>%
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096,  8102, 8097, 8098, 8101, 8099) ~ "RTIOXA_43",
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
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% body weight change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_bw2

#plot adiposity index changes

plot_43_5D_ai2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% adiposity index change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_ai2

#plot lean mass changes

plot_43_5D_lean2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% lean mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_lean2

#plot fat mass changes

plot_43_5D_fat2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% fat mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_fat2

# combined plot for RTIOXA-43 x 5 days collapsed

combined_plot2 <- (plot_43_5D_bw2 | plot_43_5D_ai2) /
  (plot_43_5D_fat2 | plot_43_5D_lean2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

combined_plot2

#data analysis
#1 does rtioxa-43 decreases bw after 5 days of injection? (collapsed)
#to check if variances differ between groups,
echoMRI_data_432 <- echoMRI_data_432 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_432) #p = 0.88 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_432) #p = 0.72 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_432) #p = 0.13 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_432) #p = 0.02 < 0.05, so variances are not considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_432)) #p = 0.64
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_432)) #p = 0.08 it seems the significance is given per ID 8098
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_432)) #p = 0.71
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_432)) #p = 0.10

#conclusion> There is no indication that RTIOXA-43 (collapsed) decreases body weight, fat mass, 
# or adiposity after 5 days. Probably because data are underpowered (n=3 per vehicle group and n=6 per RTIOXA 43 group).
# It seems the potential significance of RTIOXA=43 in fat mass is because ID 8098
#OUTLAYER ANALYSIS####

#lets check if that is true and if we have outlayers in delta fat subgroup of data

subdata <- echoMRI_data_432 %>% 
  filter(DRUG =="RTIOXA_43") %>% 
  select(ID,DRUG,delta_fat,delta_adiposity_index)

ggplot(subdata, aes(y = delta_fat, x = "")) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
  geom_point(aes(label = ID), position = position_jitter(width = 0.1)) +
  geom_text_repel(aes(label = ID)) +
  labs(y = "% fat mass change", x = "") +
  theme_minimal() #IT SEEMS LIKE ID 8098 IS AN OUTLIER

Q1 <- quantile(subdata$delta_fat, 0.25)
Q3 <- quantile(subdata$delta_fat, 0.75)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

subdata %>%
  mutate(outlier_IQR = delta_fat < lower_bound | delta_fat > upper_bound) 

#conclusion 
#8098 IS an outlayer for %change in fat mass and adiposity index


#so lets run the analysis without ID 8098
# group assignation without ID 8098 -----
## Differentiated analysis for RTIOXA 43 obtained from Yanan and from medchem----

# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8096, 8102,8099) ~ "RTIOXA_43_Z", 
      ID %in% c(8097, 8101) ~ "RTIOXA_43_M",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M"))) %>% 
  drop_na()

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

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_43_assignation) #p = 0.22 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation$adiposity_index,
  g = echoMRI_data_43_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between any of the three groups

# RTIOXA-43 × 5 days ----
## Differentiated analysis for RTIOXA 43 obtained from Yanan and from medchem----


echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>%
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_Z",
      ID %in% c(8097, 8101) ~ "RTIOXA_43_M", #8098 bye for now
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
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_Z", "RTIOXA_43_M"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% body weight change") +
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
  labs(x = NULL, y = "% adiposity index change") +
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
  labs(x = NULL, y = "% lean mass change") +
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
  labs(x = NULL, y = "% fat mass change") +
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
#1 does rtioxa-43 decreases bw after 5 days of injection? WITHOUT 8098
#to check if variances differ between groups,
echoMRI_data_43 <- echoMRI_data_43 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_43) #p = 0.65 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_43) #p = 0.46 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_43) #p = 0.07 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_43) #p = 0.74 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_43)) #p = 0.1
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_43))#p = 0.31
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_43))#p = 0.28 
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_43))#p = 0.33

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, 
# or adiposity after 5 days. Probably because data are underpowered (n=3 per group).

# group assignation -----
## Collapsed RTIOXAs 43 WITHOUT 8098----

# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation2 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION2 = case_when(
      ID %in% c(8096, 8102, 8097, 8099, 8101) ~ "RTIOXA_43", #8098 eliminated for now
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION2 = factor(DRUG_ASSIGNATION2, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

plot_43_5D_assignation2 <- ggplot(echoMRI_data_43_assignation2, aes(x = DRUG_ASSIGNATION2, y = adiposity_index, fill = DRUG_ASSIGNATION2)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43_Z" = "#66C2A5", "RTIOXA_43_M" = "darkgreen")) +
  theme(legend.position = "none")

plot_43_5D_assignation2

#pairwise comparison to evaluate if the adiposity index is the same among the groups at the baseline

#to check if variances differ between groups,
echoMRI_data_43_assignation2 <- echoMRI_data_43_assignation2 %>%
  mutate(DRUG_ASSIGNATION2 = factor(DRUG_ASSIGNATION2))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION2, data = echoMRI_data_43_assignation2) #p = 0.27 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation2$adiposity_index,
  g = echoMRI_data_43_assignation2$DRUG_ASSIGNATION2,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups

# RTIOXA-43 × 5 days (collapsed) ---- 

echoMRI_data_432 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 9) %>%
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096,  8102, 8097, 8099, 8101) ~ "RTIOXA_43", #8098 eliminated for now
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
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% body weight change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")#+
  #geom_text_repel(aes(label = ID)) 

plot_43_5D_bw2

#plot adiposity index changes

plot_43_5D_ai2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% adiposity index change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none") 

plot_43_5D_ai2

#plot lean mass changes

plot_43_5D_lean2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% lean mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_lean2

#plot fat mass changes

plot_43_5D_fat2 <- ggplot(echoMRI_data_432, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% fat mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_fat2

# combined plot for RTIOXA-43 x 5 days collapsed

combined_plot2 <- (plot_43_5D_bw2 | plot_43_5D_ai2) /
  (plot_43_5D_fat2 | plot_43_5D_lean2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

combined_plot2

#data analysis
#1 does rtioxa-43 decreases bw after 5 days of injection? (collapsed)
#to check if variances differ between groups,
echoMRI_data_432 <- echoMRI_data_432 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_432) #p = 0.34 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_432) #p = 0.34 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_432) #p = 0.23 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_432) #p = 0.19 < 0.05, so variances are not considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_432)) #p = 0.89
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_432)) #p = 0.12
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_432)) #p = 0.97
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_432)) #p = 0.11

#conclusion> the potential reduction in fat mass change in RTIOXA-43 group was per ID 8098 that actually is an outlayer

#Food intake analysis for RTIOXA 43 x 5 days no data----
