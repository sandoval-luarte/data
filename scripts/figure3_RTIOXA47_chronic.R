#We aim to explore changes in body composition (bw, fat mass, lean mass and adiposity index)
#and food intake 
#in adult females C57BL6J and NZO mice and males C57BL6J after 4 weeks of IP injection of RTIOXA-47

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
library(tidyverse)
library(ggrepel)

#format plot
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "top",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

# group assignation -----

echoMRI_data_47_chronic_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 1 & COHORT < 6) %>%
  filter(!ID %in% c(3712, 3715)) %>% #died over the experiment
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
    DRUG_ASSIGNATION = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index, GROUP, DRUG_ASSIGNATION, SEX, STRAIN,n_measurement,DIET_FORMULA) %>% 
mutate(STATUS = case_when(
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
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_47"))) %>% 
  filter(STATUS =="BW maintenance")

echoMRI_data_47_chronic_assignation %>%
  group_by(STRAIN,SEX,DIET_FORMULA,STATUS,GROUP) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

plot_47_chronic_assignation <- ggplot(echoMRI_data_47_chronic_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index prior injections") + 
 # format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange"))+
  facet_wrap(~STRAIN*GROUP)+
  theme(legend.position = "none")

plot_47_chronic_assignation

#pairwise comparison to evaluate if the adiposity index is the same among groups prior to injections

#to check if variances differ between groups,
echoMRI_data_47_chronic_assignation <- echoMRI_data_47_chronic_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

levene_results <- echoMRI_data_47_chronic_assignation %>%
  group_by(STRAIN,GROUP) %>%
  do(tidy(leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = .))) %>%
  ungroup()

levene_results #This means none of the STRAIN X GROUP groups have different variance between vehicle and RTIOXA_47.
#So diet and sex here are not important (I collapsed males and females C57)

pairwise_results <- echoMRI_data_47_chronic_assignation %>%
  group_by(STRAIN,GROUP) %>%
  group_modify(~ {
    
    # run Levene test for this subgroup
    lev <- car::leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = .x)
    equal_var <- lev$`Pr(>F)`[1] > 0.05   # TRUE = equal variances
    
    # run appropriate pairwise t-test
    test <- pairwise.t.test(
      x = .x$adiposity_index,
      g = .x$DRUG_ASSIGNATION,
      p.adjust.method = "bonferroni",
      pool.sd = equal_var
    )
    
    broom::tidy(test)
  }) %>%
  ungroup()
pairwise_results #there is no significant differences in adiposity index within each group prior to injections, which is good

#  RTIOXA-47 × 5 weeks -----

echoMRI_data_47_chronic<- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT > 1 & COHORT < 6) %>%
  filter(!ID %in% c(3712, 3715)) %>% #died over the experiment
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
  select(ID, Date, Fat, Lean, Weight, adiposity_index, GROUP, DRUG, SEX, STRAIN, STATUS, n_measurement,DIET_FORMULA) %>% 
  group_by(ID, DRUG, SEX, STRAIN, DIET_FORMULA,GROUP) %>% 
  summarise(
    delta_adiposity_index = (adiposity_index[STATUS == "BW regain"] - adiposity_index[STATUS == "BW maintenance"])/adiposity_index[STATUS == "BW maintenance"],
    delta_lean = (Lean[STATUS == "BW regain"] - Lean[STATUS == "BW maintenance"])/Lean[STATUS == "BW maintenance"],
    delta_fat = (Fat[STATUS == "BW regain"] - Fat[STATUS == "BW maintenance"])/Fat[STATUS == "BW maintenance"],
    delta_bw = (Weight[STATUS == "BW regain"] - Weight[STATUS == "BW maintenance"])/Weight[STATUS == "BW maintenance"]) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki"))

echoMRI_data_47_chronic %>%
  group_by(STRAIN,GROUP,DRUG,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) #this is good, we collapsed SEX in C57BL6J

#plot bw changes

plot_47_chronic_bw <- ggplot(echoMRI_data_47_chronic, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% body weight gain") +
#  format.plot+
  facet_wrap(~STRAIN*GROUP)+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_chronic_bw

#plot adiposity index changes

plot_47_chronic_ai <- ggplot(echoMRI_data_47_chronic, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% adiposity index gain") +
 # format.plot+
  facet_wrap(~STRAIN*GROUP)+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_chronic_ai

#plot lean mass changes

plot_47_chronic_lean <- ggplot(echoMRI_data_47_chronic, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% lean mass gain") +
  #format.plot+
  facet_wrap(~STRAIN*GROUP)+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_chronic_lean

#plot fat mass changes

plot_47_chronic_fat <- ggplot(echoMRI_data_47_chronic, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% fat mass gain") +
 # format.plot+
  facet_wrap(~STRAIN*GROUP)+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_chronic_fat

#data analysis
#1 does rtioxa-47 decreases bw after 5 weeks of injection within each group?
#to check if variances differ between groups,
echoMRI_data_47_chronic <- echoMRI_data_47_chronic %>%
  mutate(DRUG = factor(DRUG))

delta_vars <- c("delta_adiposity_index", "delta_fat", "delta_lean", "delta_bw")

levene_results <- echoMRI_data_47_chronic %>%
  group_by(STRAIN, GROUP) %>%
  group_modify(~{
    map_dfr(delta_vars, function(varname) {
      test <- leveneTest(reformulate("DRUG", varname), data = .x)
      tidy(test) %>%
        mutate(variable = varname)
    })
  }) %>%
  ungroup()

levene_significant <- levene_results %>%
  filter(p.value < 0.05)
View(levene_significant) #good so ANOVA assumptions are met


anova_results <- echoMRI_data_47_chronic %>%
  group_by(STRAIN, GROUP) %>%
  group_modify(~{
    map_dfr(
      c("delta_bw", "delta_fat", "delta_lean", "delta_adiposity_index"),
      function(varname){
        test <- aov(reformulate("DRUG", varname), data = .x)
        tidy(test) %>%
          mutate(variable = varname)
      }
    )
  }) %>%
  ungroup()

anova_significant <- anova_results %>%
  filter(p.value < 0.05)
anova_significant


#conclusion> There is no indication that RTIOXA-47 decreases body weight, fat mass, 
# or adiposity after 5 weeks of injection in any strain.
#a trend to increased delta bw, delta fat and delta adiposity index ocurred
#that goes in opposite direction of what we want.

#food intake analysis ----

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>%
  filter(!ID %in% c(3712, 3715)) %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876,
                7879, 7880, 7881, 7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728,
                3729, 7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )
  ) %>% 
  mutate(
    inj_day = min(DATE[str_detect(COMMENTS, regex("DAY_1_INJECTIONS", ignore_case = TRUE))], na.rm = TRUE)
  ) %>%
  # remove IDs that never had DAY_1_INJECTIONS
  filter(!is.infinite(inj_day)) %>%
  filter(DATE >= inj_day) %>%
  ungroup() %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki"))%>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
select(ID,BW,DATE,COMMENTS,SEX,STRAIN,DIET_FORMULA,GROUP,DRUG )

BW_data2 <- BW_data %>%
  filter(COMMENTS %in% c("DAY_1_INJECTIONS", "DAY_4_SABLE_AND_SAC")) %>% 
  group_by(ID, DRUG, SEX, STRAIN, DIET_FORMULA,GROUP) %>% 
  summarise(
    delta_BW = (BW[COMMENTS == "DAY_4_SABLE_AND_SAC"] - BW[COMMENTS == "DAY_1_INJECTIONS"])/BW[COMMENTS == "DAY_1_INJECTIONS"])
  
BW_data2 %>%
  group_by(STRAIN,SEX,GROUP,DRUG,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

plot_47_chronic_bw2 <- ggplot(BW_data2, aes(x = DRUG, y = delta_BW, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% body weight gain2") +
  #  format.plot+
  facet_wrap(~STRAIN*GROUP)+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_chronic_bw2 #both graphs showed the same trend

#FI data ----
FI_data <- read_csv("../data/FI.csv") %>% 
filter(COHORT > 1 & COHORT < 6) %>%
  filter(!ID %in% c(3712, 3715)) %>% #died over the experiment
   group_by(ID) %>% 
   arrange(DATE) %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876,
                7879, 7880, 7881, 7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728,
                3729, 7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )
  ) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>%
  mutate(
    inj_day = min(DATE[str_detect(COMMENTS, regex("DAY_1_INJECTIONS", ignore_case = TRUE))], na.rm = TRUE)
  ) %>%
  # remove IDs that never had DAY_1_INJECTIONS
  filter(!is.infinite(inj_day)) %>%
  filter(DATE >= inj_day) %>%
  ungroup() %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki"))%>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  select(ID,DIET_FORMULA,DATE,COMMENTS,SEX,corrected_intake_kcal, STRAIN,GROUP,DRUG )

# left join FI_data with BW_data by ID and COMMENTS
 FI_BW_joined <- FI_data %>%
   left_join(BW_data %>% select(ID, DATE, BW), by = c("ID", "DATE")) %>% 
   drop_na(BW) %>% 
   ungroup()
   
 FI_BW_joined %>%
   count(ID) %>% 
   print(n = Inf)
 
 FI_BW_joined%>%
   group_by(STRAIN
            ) %>%
   summarise(n_ID = n_distinct(ID)) #this is good

# # standarization of daily food intake per daily BW
 FI_BW_joined <- FI_BW_joined %>%
   group_by(ID) %>% 
   mutate(kcal_per_gBW = corrected_intake_kcal / BW)
 
# # cumulative sum per id
 FI_BW_joined_sum <- FI_BW_joined %>%
   group_by(ID, DRUG,STRAIN,GROUP) %>%   # group by ID and DRUG
   summarise(
     total_kcal = sum(kcal_per_gBW, na.rm = TRUE),
     .groups = "drop"
   )
 
 FI_BW_joined_sum %>%
   group_by(STRAIN) %>%
   summarise(n_ID = n_distinct(ID)) #this is good
 
 
plot_47_foodintake <- ggplot(FI_BW_joined_sum, aes(x = DRUG, y = total_kcal, fill = DRUG)) +
   stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
   geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
   theme_minimal() +
   labs(x = NULL, y = "Total kcal per g body weight over 5 weeks") +
   scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
   theme(legend.position = "none")+
  facet_wrap(~STRAIN)+
   geom_text_repel(aes(label = ID)) 
 plot_47_foodintake

#it seems there is no significant effects of RTIOXA 47 in food intake after 4 weeks of treatment



