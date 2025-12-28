# This script aims to explore changes in body weight (BW) in females and males CD1 dams
#exposed to perinatal BPA 5 ug/kg and fed with HFD (D12451i, research diets)
#or LFD (D12451Hi, research diets)

#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(lmerTest)
library(emmeans)
library(pracma)
library(lubridate)
library(broom)
library(stringr)
library(forcats)
library(patchwork)


 # BW over time data----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(15, 16)) %>% 
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) %>% 
  left_join(METABPA, by= "ID") %>% 
  select(
    -SEX.y,
    -DIET_FORMULA.y
  ) %>% 
  rename(
    SEX = SEX.x,
    DIET_FORMULA = DIET_FORMULA.x
  )

BW_data %>% 
  group_by(SEX,BPA_EXPOSURE, DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 


BW_summary <- BW_data %>%
  group_by(day_rel,BPA_EXPOSURE,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#plot 1

plot_bw_sex <- ggplot(
  BW_summary,
  aes(x = day_rel,
      y = mean_BW,
      color = BPA_EXPOSURE,
      fill  = BPA_EXPOSURE)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = mean_BW - sem_BW,
        ymax = mean_BW + sem_BW),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

plot_bw_sex

#stats LME----
# Females
female_data <- BW_data %>% filter(SEX == "F")
lme_female <- lmer(BW ~ BPA_EXPOSURE * day_rel + (1|ID), data = female_data)

# Males
male_data <- BW_data %>% filter(SEX == "M")
lme_male <- lmer(BW ~ BPA_EXPOSURE * day_rel + (1|ID), data = male_data)

# Females
# Round day_rel to nearest measured day
emmeans_f <- emmeans(lme_female, ~ BPA_EXPOSURE | day_rel, at = list(day_rel = unique(female_data$day_rel)))
# Males
# Round day_rel to nearest measured day
emmeans_m <- emmeans(lme_male, ~ BPA_EXPOSURE | day_rel, at = list(day_rel = unique(male_data$day_rel)))

contrast_f <- contrast(emmeans_f, method = "pairwise") %>% as.data.frame()
contrast_m <- contrast(emmeans_m, method = "pairwise") %>% as.data.frame()

# Make sure all contrasts are YES - NO females
contrast_f <- contrast_f %>%
  mutate(
    estimate = ifelse(contrast == "NO - YES", -estimate, estimate),
    t.ratio = ifelse(contrast == "NO - YES", -t.ratio, t.ratio),
    contrast = "YES - NO"
  )
sig_daysf <- contrast_f %>% filter(p.value < 0.05)

# Make sure all contrasts are YES - NO males
contrast_m <- contrast_m %>%
  mutate(
    estimate = ifelse(contrast == "NO - YES", -estimate, estimate),
    t.ratio = ifelse(contrast == "NO - YES", -t.ratio, t.ratio),
    contrast = "YES - NO"
  )
sig_daysm <- contrast_m %>% filter(p.value < 0.05)

first_sig_dayf <- sig_daysf %>% slice_min(day_rel)
first_sig_dayf$day_rel

first_sig_daym <- sig_daysm %>% slice_min(day_rel)
first_sig_daym$day_rel


# Define first significant day for females only
first_sig_day_f <- tibble(
  SEX = "F",
  day_sig = 56
)

# Compute ymin and ymax for females only, keeping SEX
sig_rect <- first_sig_day_f %>%
  left_join(
    BW_summary %>%
      filter(SEX == "F") %>%
      group_by(SEX) %>%
      summarise(ymin = min(mean_BW - sem_BW),
                ymax = max(mean_BW + sem_BW),
                .groups = "drop"),   # keep SEX
    by = "SEX"
  ) %>%
  mutate(
    xmin = day_sig,
    xmax = max(BW_summary$day_rel)
  )

# Plot 2----
plot_bw_sig_sex <- ggplot(
  BW_summary,
  aes(x = day_rel, y = mean_BW, color = BPA_EXPOSURE, fill = BPA_EXPOSURE)
) +
  geom_rect(
    data = sig_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "red", alpha = 0.1
  ) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = mean_BW - sem_BW, ymax = mean_BW + sem_BW),
    alpha = 0.25, color = NA
  ) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)
plot_bw_sig_sex 

##-- BW over time with diet formula just for F ----

BW_summary_diet <- BW_data %>%
  filter(SEX=="F") %>% 
  group_by(day_rel,BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#stats LMER for BW----
# Females
female_data <- BW_data %>% filter(SEX == "F")
lme_femalediet <- lmer(BW ~ BPA_EXPOSURE * DIET_FORMULA * day_rel + (1|ID), data = female_data)
# Females
# Round day_rel to nearest measured day
emmeans_fdiet <- emmeans(lme_femalediet, ~ BPA_EXPOSURE | day_rel * DIET_FORMULA)
contrast_fdiet <- contrast(emmeans_fdiet, method = "pairwise") %>% as.data.frame()

# Make sure all contrasts are YES - NO females
contrast_fdiet <- contrast_fdiet %>%
  mutate(
    estimate = ifelse(contrast == "NO - YES", -estimate, estimate),
    t.ratio = ifelse(contrast == "NO - YES", -t.ratio, t.ratio),
    contrast = "YES - NO"
  )
sig_daysfdiet <- contrast_fdiet %>% filter(p.value < 0.05)


# BW over time summary for females by diet
BW_summary_diet <- BW_data %>%
  filter(SEX=="F") %>% 
  group_by(day_rel, BPA_EXPOSURE, DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# First significant day per diet from your LME contrasts
first_sig_dayf_diet <- sig_daysfdiet %>%
  group_by(DIET_FORMULA) %>%
  slice_min(day_rel) %>%
  ungroup()

# Create shaded rectangle data per diet
sig_rectdiet <- BW_summary_diet %>%
  group_by(DIET_FORMULA) %>%
  summarise(
    ymin = min(mean_BW - sem_BW),
    ymax = max(mean_BW + sem_BW),
    .groups = "drop"
  ) %>%
  left_join(first_sig_dayf_diet %>% select(DIET_FORMULA, day_sig = day_rel),
            by = "DIET_FORMULA") %>%
  mutate(
    xmin = day_sig,
    xmax = max(BW_summary_diet$day_rel)
  )

# Plot 3----
plot_bw_female_diet <- ggplot(
  BW_summary_diet,
  aes(x = day_rel, y = mean_BW, color = BPA_EXPOSURE, fill = BPA_EXPOSURE)
) +
  geom_rect(
    data = sig_rectdiet,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "red", alpha = 0.1
  ) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = mean_BW - sem_BW, ymax = mean_BW + sem_BW),
    alpha = 0.25, color = NA
  ) +
  facet_wrap(~ DIET_FORMULA) +
  labs(
    title = "                Females",
    x = "Days relative to first measurement",
    y = "Body weight (g)",
    color = "BPA exposure",
    fill = "BPA exposure"
  ) +
  theme_classic(base_size = 14) 

plot_bw_female_diet 

# combined plots----
combined_plot_bw <- (plot_bw_sig_sex / plot_bw_female_diet) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

combined_plot_bw #figure 2

# check if BW of the four groups females with or without HFD and with or without BPA exposure had the same BW at day 0

bw_day0_female <- BW_data %>%
  filter(
    SEX == "F",
    day_rel == 0
  )

bw_day0_female %>%
  group_by(DIET_FORMULA, BPA_EXPOSURE) %>%
  summarise(n = n(), .groups = "drop")

lm_day0 <- lm(
  BW ~ BPA_EXPOSURE * DIET_FORMULA,
  data = bw_day0_female
)

anova(lm_day0)

emmeans_day0 <- emmeans(lm_day0, ~ BPA_EXPOSURE | DIET_FORMULA)
contrast(emmeans_day0, method = "pairwise")

# BW gain over time data----

BW_gain <- BW_data %>% 
  filter(SEX =="F") %>% 
  select(ID,day_rel,bw_rel,BPA_EXPOSURE,DIET_FORMULA,SEX)

BW_gain %>% 
  group_by(BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

BW_gainsummary <- BW_gain  %>%
  group_by(day_rel,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(
    mean_BWgain = mean(bw_rel, na.rm = TRUE),
    sem_BWgain  = sd(bw_rel, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ggplot(BW_gainsummary,
       aes(x = day_rel,
           y = mean_BWgain,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_BWgain - sem_BWgain,
                  ymax = mean_BWgain + sem_BWgain),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~DIET_FORMULA) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight gain",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

#stats LME----
# Females
lme_femalegain <- lmer(bw_rel ~ BPA_EXPOSURE * day_rel + (1|ID), data = BW_gain)
# Females
# Round day_rel to nearest measured day
emmeans_fgain <- emmeans(lme_femalegain, ~ BPA_EXPOSURE | day_rel, at = list(day_rel = unique(BW_gain$day_rel)))
contrast_fgain <- contrast(emmeans_fgain, method = "pairwise") %>% as.data.frame()

# Make sure all contrasts are YES - NO females
contrast_fgain <- contrast_fgain %>%
  mutate(
    estimate = ifelse(contrast == "NO - YES", -estimate, estimate),
    t.ratio = ifelse(contrast == "NO - YES", -t.ratio, t.ratio),
    contrast = "YES - NO"
  )
sig_daysfgain <- contrast_fgain %>% filter(p.value < 0.05)

first_sig_dayfgain <- sig_daysfgain %>% slice_min(day_rel)
first_sig_dayfgain$day_rel #there is no day in which BW gain of BPA females are significant higher than non BPA female in both diets

# Speed = rate of change of body weight over time data----

BW_speed <- BW_data %>%
  group_by(ID, SEX, BPA_EXPOSURE) %>%
  do(tidy(lm(BW ~ day_rel, data = .))) %>%
  filter(term == "day_rel") %>%
  rename(
    speed_g_per_day = estimate,
    se = std.error
  )

BW_speed %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 

# plot 1 BPA within sex ----
p1 <- ggplot(BW_speed,
             aes(x = BPA_EXPOSURE,
                 y = speed_g_per_day,
                 fill = BPA_EXPOSURE)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_grid(~ SEX) +
  labs(
    y = "BW gain speed (g/day)",
    x = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

## stats speed----

lmer_speed <- lmer(
  BW ~ day_rel * BPA_EXPOSURE * SEX + (day_rel | ID),
  data = BW_data
)
summary(lmer_speed)


#Male mice exhibited a significantly faster rate of body
#weight gain than females (≈0.08 g/day greater; t = 3.54).
#BPA exposure did not significantly alter growth rate in 
#either sex, and no significant BPA × sex interaction on
#growth speed was observed.

#BPA-exposed females are estimated to gain slightly 
#faster than BPA-free females. 
# BUT t = 1.36 so this is not statistically significant


# ---- Plot 2: Sex difference ----
p2 <- ggplot(BW_speed,
             aes(x = SEX,
                 y = speed_g_per_day)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = BPA_EXPOSURE),
              width = 0.15, size = 2, alpha = 0.7) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("F", "M")),
    label = "p.format"
  ) +
  labs(
    y = "BW gain speed (g/day)",
    x = "Sex",
    color = "BPA exposure"
  ) +
  theme_classic(base_size = 14)


combined_plot <- p1 | p2
combined_plot

# Body comp over time ----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16)) %>% 
  filter(!ID %in% c(9354, 9367, 9368, 9372, 9414)) %>%  # 9367 and 9368 have confused data from 8/1/25 echoMRI 
  # 9354 and 9372 lack of complete schedule of measurements in echoMRI
  # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
  group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index,COHORT,DIET_FORMULA) %>% 
  left_join(METABPA, by= "ID") %>% 
  ungroup() %>% 
  select(
    -DIET_FORMULA.y) %>% 
  rename(
    DIET_FORMULA = DIET_FORMULA.x)
  

#cohort 15 are 24 animals originally but if I eliminate 4 animals so total animals are 20 per date (9367, 9368, 9354, 9372)
#cohort 16 are 15 animals originally but if I eliminate 1 animal so total animals are 14 per date (9414)

echoMRI_data %>% 
  group_by(SEX) %>%
  summarise(n_ID = n_distinct(ID)) 


echoMRI_data_comparisons <- echoMRI_data %>% 
  mutate(
    n_measurement= case_when(
      COHORT == 16 & Date == "2025-08-04" ~ "0 days with O.D",
      COHORT == 16 & Date == "2025-09-09" ~ "35 days with O.D",
      COHORT == 16 & Date == "2025-10-07" ~ "63 days with O.D",
      COHORT == 15 & Date == "2025-08-01" ~ "98 days with O.D",
      COHORT == 16 & Date == "2025-11-17" ~ "98 days with O.D",
      COHORT == 15 & Date == "2025-09-09" ~ "133 days with O.D",
      COHORT == 16 & Date == "2025-12-18" ~ "133 days with O.D",
      COHORT == 15 & Date == "2025-10-07" ~ "150 days with O.D"))

bodycomp_summary <- echoMRI_data_comparisons %>%
  group_by(n_measurement, BPA_EXPOSURE,SEX) %>%
  summarise(
    mean_AI = mean(adiposity_index, na.rm = TRUE),
    sem_AI  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    mean_fat = mean(Fat, na.rm = TRUE),
    sem_fat  = sd(Fat, na.rm = TRUE) / sqrt(n()),
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean  = sd(Lean, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>% 
  mutate(n_measurement = factor(n_measurement, levels = c("0 days with O.D",
                                                          "35 days with O.D",
                                                          "63 days with O.D",
                                                          "98 days with O.D",
                                                          "133 days with O.D",
                                                          "150 days with O.D"))) %>% 
  ungroup()


## plot 1 males and females ----
p1 <- ggplot(bodycomp_summary,
             aes(x = n_measurement,
                 y = mean_AI,
                 color = BPA_EXPOSURE,
                 fill  = BPA_EXPOSURE,
                 group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_AI - sem_AI,
                  ymax = mean_AI + sem_AI),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days of measurement",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

# LME model: AI ~ time * BPA_EXPOSURE * SEX + random intercept for ID
lmer_AI <- lmer(adiposity_index ~ n_measurement * BPA_EXPOSURE * SEX +
                  (1 | ID),
                data = echoMRI_data_comparisons)
emm_AI <- emmeans(lmer_AI, ~ BPA_EXPOSURE | SEX * n_measurement)
AI_contrasts <- contrast(emm_AI, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  select(SEX, n_measurement, contrast, estimate, SE, df, t.ratio, p.value)
AI_first_sig <- AI_contrasts %>%
  filter(p.value < 0.05) %>%
  group_by(SEX) %>%
  slice_min(n_measurement)  # first day with significant difference
AI_contrasts

## AI with diet_formula FEMALES----

bodycomp_summarydiet <- echoMRI_data_comparisons %>%
  filter(SEX=="F") %>% 
  group_by(n_measurement, BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(
    mean_AI = mean(adiposity_index, na.rm = TRUE),
    sem_AI  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    mean_fat = mean(Fat, na.rm = TRUE),
    sem_fat  = sd(Fat, na.rm = TRUE) / sqrt(n()),
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean  = sd(Lean, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>% 
  mutate(n_measurement = factor(n_measurement, levels = c("0 days with O.D",
                                                          "35 days with O.D",
                                                          "63 days with O.D",
                                                          "98 days with O.D",
                                                          "133 days with O.D",
                                                          "150 days with O.D"))) %>% 
  filter(n_measurement %in% c("98 days with O.D",
                              "133 days with O.D",
                              "150 days with O.D")) %>% 
  ungroup()

## plot 2 females per diet formula ----
p2 <- ggplot(bodycomp_summarydiet,
             aes(x = n_measurement,
                 y = mean_AI,
                 color = BPA_EXPOSURE,
                 fill  = BPA_EXPOSURE,
                 group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_AI - sem_AI,
                  ymax = mean_AI + sem_AI),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ DIET_FORMULA) +
  labs(
    title = "                      Females",
    x = "Days of measurement",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots vertically
combined_plot <- p1 / p2
combined_plot

## LMER females AI from day 98----

# Filter for females and days >= 98

echoMRI_females_post98 <- echoMRI_data_comparisons %>%
  filter(SEX == "F") %>%
  mutate(n_measurement_days = case_when(
    n_measurement == "0 days with O.D"   ~ 0,
    n_measurement == "35 days with O.D"  ~ 35,
    n_measurement == "63 days with O.D"  ~ 63,
    n_measurement == "98 days with O.D"  ~ 98,
    n_measurement == "133 days with O.D" ~ 133,
    n_measurement == "150 days with O.D" ~ 150
  )) %>%
  filter(n_measurement_days >= 98)

# Fit LME model including diet
lmer_AI_females_post98 <- lmer(adiposity_index ~ n_measurement * BPA_EXPOSURE + DIET_FORMULA + (1 | ID),
                               data = echoMRI_females_post98)

emm_AI_females_post98 <- emmeans(lmer_AI_females_post98, ~ BPA_EXPOSURE | n_measurement * DIET_FORMULA)
AI_contrasts_females_post98 <- contrast(emm_AI_females_post98, method = "pairwise", adjust = "none") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  select(n_measurement, DIET_FORMULA, contrast, estimate, SE, df, t.ratio, p.value)
AI_contrasts_females_post98

## fat mass ----
## plot 1 males and females ----
p1fat <- ggplot(bodycomp_summary,
             aes(x = n_measurement,
                 y = mean_fat,
                 color = BPA_EXPOSURE,
                 fill  = BPA_EXPOSURE,
                 group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_fat - sem_fat,
                  ymax = mean_fat + sem_fat),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days of measurement",
    y = "fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1fat

# LME model: AI ~ time * BPA_EXPOSURE * SEX + random intercept for ID
lmer_fat <- lmer(Fat ~ n_measurement * BPA_EXPOSURE * SEX +
                  (1 | ID),
                data = echoMRI_data_comparisons)
emm_fat <- emmeans(lmer_fat, ~ BPA_EXPOSURE | SEX * n_measurement)
fat_contrasts <- contrast(emm_fat, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  select(SEX, n_measurement, contrast, estimate, SE, df, t.ratio, p.value)
fat_first_sig <- fat_contrasts %>%
  filter(p.value < 0.05) %>%
  group_by(SEX) %>%
  slice_min(n_measurement)  # first day with significant difference
fat_contrasts

## lean mass ----
## plot 1 males and females ----
p1lean <- ggplot(bodycomp_summary,
             aes(x = n_measurement,
                 y = mean_lean,
                 color = BPA_EXPOSURE,
                 fill  = BPA_EXPOSURE,
                 group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_lean - sem_lean,
                  ymax = mean_lean + sem_lean),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days of measurement",
    y = "lean mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1lean

# fat and lean plot combined----

# Combine plots vertically
combined_plot <- p1fat / p1lean
combined_plot


# LME model: AI ~ time * BPA_EXPOSURE * SEX + random intercept for ID
lmer_lean <- lmer(Lean ~ n_measurement * BPA_EXPOSURE * SEX +
                   (1 | ID),
                 data = echoMRI_data_comparisons)
emm_lean <- emmeans(lmer_lean, ~ BPA_EXPOSURE | SEX * n_measurement)
lean_contrasts <- contrast(emm_lean, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  select(SEX, n_measurement, contrast, estimate, SE, df, t.ratio, p.value)
lean_first_sig <- lean_contrasts %>%
  filter(p.value < 0.05) %>%
  group_by(SEX) %>%
  slice_min(n_measurement)  # first day with significant difference
lean_contrasts


# OGTT----

#OGTT WAS DONE JUST FOR COHORT 15 WHEN ALL ANIMALS WERE 27 WEEK OLD

OGTT <- read_csv("~/Documents/GitHub/data/data/OGTT_CD1.csv") %>% 
  mutate(DATE = mdy(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID)%>% 
  left_join(METABPA, by= "ID")

ogtt_long <- OGTT  %>% 
  pivot_longer(
    cols = starts_with("time_"),
    names_to = "time_min",
    values_to = "glucose"
  ) %>% 
  mutate(
    time_min = gsub("time_", "", time_min),  # remove "time_"
    time_min = as.numeric(time_min)           # convert to numeric
  )

ogtt_long %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 


auc_df <- ogtt_long %>% 
  arrange(ID, time_min) %>% 
  group_by(ID, BPA_EXPOSURE,SEX) %>% 
  summarise(
    AUC = trapz(time_min, glucose),
    .groups = "drop"
  )

ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
 facet_wrap(~ SEX) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()

## two way anova for AUC (SEX x BPA)----

auc_df$BPA_EXPOSURE <- factor(auc_df$BPA_EXPOSURE)
auc_df$SEX <- factor(auc_df$SEX)

model <- aov(AUC ~ BPA_EXPOSURE * SEX , data = auc_df)
summary(model)

#Males clearly have higher AUCs than females (F) overall
#There is no statistically supported difference between
#BPA-exposed and non-exposed females (p=0.52) non-significant BPA × SEX interaction
#there is a trend (p=0.06) for HFD to increased AUC


# FI over time data----

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(15,16)) %>% 
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) %>% 
  left_join(METABPA, by= "ID") %>% 
  select(
    -SEX.y,
    -DIET_FORMULA.x,
    -DIET_FORMULA.y
  ) %>% 
  rename(
    SEX = SEX.x
  ) %>% 
  mutate(
   FIcumulative = cumsum(corrected_intake_kcal)) %>% 
  ungroup()

FI_data  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID))

FI_summary <- FI_data %>%
  group_by(day_rel, SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_FI = mean(FIcumulative, na.rm = TRUE),
    sem_FI  = sd(FIcumulative, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# cumulative food intake plot----

FI_plot <- ggplot(FI_summary, aes(x = day_rel, y = mean_FI, color = BPA_EXPOSURE, fill = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_FI - sem_FI, ymax = mean_FI + sem_FI), alpha = 0.25, color = NA) +
  facet_wrap(~ SEX) +
  labs(x = "Days", y = "Cumulative food intake (kcal)", color = "BPA exposure", fill = "BPA exposure") +
  theme_classic(base_size = 14)
FI_plot

mod_FI <- lmer(
  FIcumulative ~ day_rel * BPA_EXPOSURE * SEX +
    (1 | ID),
  data = FI_data
)
summary(mod_FI)
anova(mod_FI)

lmer(FIcumulative ~ day_rel * BPA_EXPOSURE + (1|ID),
     data = filter(FI_data, SEX == "F"))

lmer(FIcumulative ~ day_rel * BPA_EXPOSURE + (1|ID),
     data = filter(FI_data, SEX == "M"))

summary(lmer(FIcumulative ~ day_rel * BPA_EXPOSURE + (1|ID),
             data = filter(FI_data, SEX == "F")))

summary(lmer(FIcumulative ~ day_rel * BPA_EXPOSURE + (1|ID),
             data = filter(FI_data, SEX == "M")))

mod_F_F <- lmer(
  FIcumulative ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = filter(FI_data, SEX == "F")
)

emm_F <- emmeans(
  mod_F_F,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = c(30, 60, 90))
)

contrast(emm_F, method = "pairwise")

mod_F_M <- lmer(
  FIcumulative ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = filter(FI_data, SEX == "M")
)

emm_M <- emmeans(
  mod_F_M,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = c(30, 60, 90))
)

contrast(emm_M, method = "pairwise")


emm_F_56 <- emmeans(
  mod_F_F,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = 56)
)

contrast(emm_F_56, method = "pairwise")

emm_M_56 <- emmeans(
  mod_F_M,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = 56)
)

contrast(emm_M_56, method = "pairwise")

# energy eficiency analysis ----
EE_data <- BW_data %>%
  select(ID, SEX, BPA_EXPOSURE, day_rel, BW) %>%
  left_join(FI_data %>% select(ID, day_rel, FIcumulative), by = c("ID", "day_rel"))

# Check that IDs are still valid
unique(EE_data$ID)

intervals <- c(0, 30, 60, 90, 120, 150)  # days

library(dplyr)
library(purrr)
library(tibble)

EE_per_interval <- EE_data %>%
  group_by(ID) %>%
  group_split() %>%  # split into list by ID
  map_dfr(function(df) {
    res <- tibble()
    intervals <- c(0, 30, 60, 90, 120, 150)
    
    for(i in 2:length(intervals)){
      start_day <- intervals[i-1]
      end_day   <- intervals[i]
      
      # Closest rows
      row_start <- df[which.min(abs(df$day_rel - start_day)), ]
      row_end   <- df[which.min(abs(df$day_rel - end_day)), ]
      
      delta_BW <- row_end$BW - row_start$BW
      delta_FI <- row_end$FIcumulative - row_start$FIcumulative
      EE_val   <- delta_BW / delta_FI
      
      res <- bind_rows(res, tibble(
        ID = df$ID[1],
        SEX = as.character(df$SEX[1]),
        BPA_EXPOSURE = as.character(df$BPA_EXPOSURE[1]),
        interval = paste0(start_day, "-", end_day),
        delta_BW = delta_BW,
        delta_FI = delta_FI,
        EE = EE_val
      ))
    }
    res
  })

unique(EE_per_interval$ID)
nrow(EE_per_interval)
head(EE_per_interval)

EE_per_interval <- EE_per_interval %>%
  mutate(
    SEX = factor(SEX, levels = c("F","M")),
    BPA_EXPOSURE = factor(BPA_EXPOSURE, levels = c("NO","YES")),
    interval = factor(interval, levels = c("0-30","30-60","60-90","90-120","120-150"))
  )

EE_summary <- EE_per_interval %>%
  group_by(interval, SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_EE = mean(EE, na.rm = TRUE),
    sem_EE  = sd(EE, na.rm = TRUE)/sqrt(n()),
    n = n(),
    .groups = "drop"
  )

EE_summary

library(lme4)
library(emmeans)

EE_lme <- lmer(EE ~ BPA_EXPOSURE * SEX + interval + (1|ID), data = EE_per_interval)
summary(EE_lme)

# Post-hoc comparisons (optional)
EE_emm <- emmeans(EE_lme, ~ BPA_EXPOSURE | SEX * interval)
pairs(EE_emm)


# Energy Efficiency plot ----
EE_plot <- ggplot(EE_summary, aes(x = interval, y = mean_EE, group = BPA_EXPOSURE, color = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_EE - sem_EE, ymax = mean_EE + sem_EE), width = 0.2) +
  facet_wrap(~SEX) +
  labs(
    x = "Interval (days)",
    y = "Energy Efficiency (g BW gain / kcal food)",
    color = "BPA exposure"  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

EE_plot


# Combine vertically
combined_plot <- FI_plot / EE_plot + plot_layout(heights = c(1, 1))

# Display
combined_plot





# indirect calorimetry columbus data ----

##counts####
ical_data_counts_coh15 <- read_csv("~/Documents/GitHub/data/data/iCal_Kotz_082425_counts.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("93", ID))) %>% 
  select(-`8/27/25 11:00`, #we want to keep just 24 hours
         -`8/27/25 12:00`,
         -`8/27/25 13:00`,
         -`8/27/25 14:00`,
         -`8/27/25 15:00`,
         -`8/27/25 16:00`,
         -`8/27/25 17:00`,
         -`8/27/25 18:00`,
         -`8/27/25 19:00`,
         -`8/28/25 20:00`,
         -`8/28/25 21:00`,
         -`8/28/25 22:00`,
         -`8/28/25 23:00`,
         -`8/29/25 0:00`,
         -`8/29/25 1:00`,
         -`8/29/25 2:00`,
         -`8/29/25 3:00`,
         -`8/29/25 4:00`,
         -`8/29/25 5:00`,
         -`8/29/25 6:00`,
         -`8/29/25 7:00`,
         -`8/29/25 8:00`,
         -`8/29/25 9:00`,
         -`8/29/25 10:00`)


ical_data_counts_coh16 <- read_csv("~/Documents/GitHub/data/data/ical analysis 11262025_counts.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("9", ID))) %>% 
  select(-`11/24/25 9:00`, #we want to keep just 24 hours
         -`11/24/25 10:00`,
         -`11/24/25 11:00`,
         -`11/24/25 12:00`,
         -`11/24/25 13:00`,
         -`11/24/25 14:00`,
         -`11/24/25 15:00`,
         -`11/24/25 16:00`,
         -`11/24/25 17:00`,
         -`11/24/25 18:00`,
         -`11/24/25 19:00`,
         -`11/25/25 20:00`,
         -`11/25/25 21:00`,
         -`11/25/25 22:00`,
         -`11/25/25 23:00`,
         -`11/26/25 0:00`,
         -`11/26/25 1:00`,
         -`11/26/25 2:00`,
         -`11/26/25 3:00`,
         -`11/26/25 4:00`,
         -`11/26/25 5:00`,
         -`11/26/25 6:00`,
         -`11/26/25 7:00`,
         -`11/26/25 8:00`,
         -`11/26/25 9:00`)

ical_long15 <- ical_data_counts_coh15 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "count") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16 <- ical_data_counts_coh16 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "count") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16 <- ical_long16 %>% 
  group_by(ID) %>%
  mutate(
    day = as.integer(as.Date(datetime) - min(as.Date(datetime))) + 1
  ) %>%
  ungroup() %>% 
  mutate(
    datetime_day = paste0(
      "day ", day, " ",
      format(datetime, "%H:%M")
    )
  )

ical_long15 <- ical_long15 %>%
  mutate(cohort = 15)

ical_long16 <- ical_long16 %>%
  mutate(cohort = 16)

common_cols <- intersect(names(ical_long15), names(ical_long16)) 

ical_long15 <- ical_long15 %>% select(all_of(common_cols))
ical_long16 <- ical_long16 %>% select(all_of(common_cols))

ical_long_all <- bind_rows(ical_long15, ical_long16) #this is key, here we combined

ical_long_all %>% 
  group_by(cohort) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have 24 animals for cohort 15 and 15 animals for cohort 16
#so in total 39 animals

ical_long_all <- ical_long_all %>%
  mutate(
    # shift so experimental day starts at 20:00
    datetime_shifted = datetime - hours(20),
    
    # experimental day (20:00 → 19:59)
    exp_day = as.integer(
      as.Date(datetime_shifted) -
        min(as.Date(datetime_shifted))
    ) + 1,
    
    # clock hour
    hour = hour(datetime),
    
    # ordering variable: 20 → 43 (20–23, then 0–19)
    hour_ordered = if_else(hour < 20, hour + 24, hour)
  ) %>%
  arrange(ID, exp_day, hour_ordered) %>%
  group_by(ID, exp_day) %>%
  mutate(count_total = cumsum(count)) %>%
  ungroup() %>% 
  filter(!ID == 9406) %>% #9406 has a very weird pattern
  group_by(ID, exp_day) %>%
   mutate(
    relative_total_count = count_total - first(count_total))

ical_long_all %>%
  filter(ID == min(ID)) %>%
  select(datetime, exp_day, hour, hour_ordered, count_total) %>%
  arrange(datetime)

ical_long_all <- ical_long_all %>%
  mutate(
    hour_label = factor(
      hour_ordered,
      levels = 20:43,
      labels = c(20:23, 0:19)
    )
  ) %>% 
  mutate(
    lights = if_else(
      hour >= 20 | hour < 6,
      "OFF",
      "ON"
    )
  )


ggplot(
  ical_long_all,
  aes(
    x = hour_label,
    y = relative_total_count,
    group = interaction(ID, exp_day)
  )
) +
  geom_line(size = 1) +
  facet_wrap(~ ID) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "relative cumulative Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 8)
  )

ical_long_allgrouped <- ical_long_all %>% 
  group_by(hour_label, SEX, BPA_EXPOSURE) %>% 
  summarise(
    mean_counts = mean(relative_total_count, na.rm = TRUE),
    sem_counts  = sd(relative_total_count, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

ggplot(
  ical_long_allgrouped,
  aes(
    x = hour_label,
    y = mean_counts,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_counts - sem_counts,
      ymax = mean_counts + sem_counts
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(size = 1) +
  facet_wrap(~ SEX) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "Relative cumulative count (mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 10)
  )


ical_long_all <- ical_long_all %>%
  group_by(ID, exp_day) %>%
  mutate(
    relative_total_counts_lights = case_when(
      
      # Lights OFF: 20:00–05:59 → start at 0 at 20:00
      lights == "OFF" ~
        relative_total_count - min(relative_total_count[lights == "OFF"], na.rm = TRUE),
      
      # Lights ON: 06:00–19:59 → start at 0 at 06:00
      lights == "ON" ~
        relative_total_count - min(relative_total_count[lights == "ON"], na.rm = TRUE)
    )
  ) %>%
  ungroup() 

  ical_long_all <- ical_long_all %>%
  mutate(
    hour_phase = case_when(
      lights == "OFF" ~ factor(
        hour,
        levels = c(20:23, 0:5)
      ),
      lights == "ON" ~ factor(
        hour,
        levels = 6:19
      )
    )
  )


ical_long_allgroupedlights <- ical_long_all %>% 
  group_by(hour_label, SEX, BPA_EXPOSURE,lights) %>% 
  summarise(
    mean_counts = mean(relative_total_counts_lights, na.rm = TRUE),
    sem_counts  = sd(relative_total_counts_lights, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

ical_long_allgroupedlights <- ical_long_allgroupedlights %>%
  mutate(
    hour_phase = case_when(
      lights == "OFF" ~ factor(hour_label, levels = c(20:23, 0:5)),
      lights == "ON"  ~ factor(hour_label, levels = 6:19)
    )
  )

ggplot(
  ical_long_allgroupedlights,
  aes(
    x = hour_phase,
    y = mean_counts,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE,
    group = BPA_EXPOSURE
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_counts - sem_counts,
      ymax = mean_counts + sem_counts
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(size = 1) +
  facet_grid(SEX ~ lights, scales = "free_x") +
  labs(
    x = "Time of day",
    y = "Relative cumulative count (mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 11)
  )

ical_auc_data <- ical_long_all %>%
  mutate(
    time_in_phase = case_when(
      lights == "OFF" ~ match(hour, c(20:23, 0:5)) - 1,
      lights == "ON"  ~ match(hour, 6:19) - 1
    )
  )
auc_per_ID <- ical_auc_data %>%
  group_by(ID, SEX, BPA_EXPOSURE, lights) %>%
  arrange(time_in_phase, .by_group = TRUE) %>%
  summarise(
    AUC = trapz(time_in_phase, relative_total_counts_lights),
    .groups = "drop"
  )

auc_wide <- auc_per_ID %>%
  tidyr::pivot_wider(
    names_from  = lights,
    values_from = AUC,
    names_prefix = "AUC_"
  )

auc_per_ID %>% count(ID, lights)

summary(auc_wide$AUC_OFF)
summary(auc_wide$AUC_ON)

auc_per_ID$ID <- factor(auc_per_ID$ID)

model_auc <- aov(
  AUC ~ SEX * BPA_EXPOSURE * lights + Error(ID / lights),
  data = auc_per_ID
)

summary(model_auc)

##heat####
ical_data_heat_coh15 <- read_csv("~/Documents/GitHub/data/data/iCal_Kotz_082425_heat.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("93", ID))) %>% 
  select(-`8/27/25 11:00`, #we want to keep just 24 hours
         -`8/27/25 12:00`,
         -`8/27/25 13:00`,
         -`8/27/25 14:00`,
         -`8/27/25 15:00`,
         -`8/27/25 16:00`,
         -`8/27/25 17:00`,
         -`8/27/25 18:00`,
         -`8/27/25 19:00`,
         -`8/28/25 20:00`,
         -`8/28/25 21:00`,
         -`8/28/25 22:00`,
         -`8/28/25 23:00`,
         -`8/29/25 0:00`,
         -`8/29/25 1:00`,
         -`8/29/25 2:00`,
         -`8/29/25 3:00`,
         -`8/29/25 4:00`,
         -`8/29/25 5:00`,
         -`8/29/25 6:00`,
         -`8/29/25 7:00`,
         -`8/29/25 8:00`,
         -`8/29/25 9:00`,
         -`8/29/25 10:00`)

ical_data_heat_coh16 <- read_csv("~/Documents/GitHub/data/data/ical analysis 11262025_heat.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("9", ID))) %>% 
  select(-`11/24/25 9:00`, #we want to keep just 24 hours
         -`11/24/25 10:00`,
         -`11/24/25 11:00`,
         -`11/24/25 12:00`,
         -`11/24/25 13:00`,
         -`11/24/25 14:00`,
         -`11/24/25 15:00`,
         -`11/24/25 16:00`,
         -`11/24/25 17:00`,
         -`11/24/25 18:00`,
         -`11/24/25 19:00`,
         -`11/25/25 20:00`,
         -`11/25/25 21:00`,
         -`11/25/25 22:00`,
         -`11/25/25 23:00`,
         -`11/26/25 0:00`,
         -`11/26/25 1:00`,
         -`11/26/25 2:00`,
         -`11/26/25 3:00`,
         -`11/26/25 4:00`,
         -`11/26/25 5:00`,
         -`11/26/25 6:00`,
         -`11/26/25 7:00`,
         -`11/26/25 8:00`,
         -`11/26/25 9:00`)

ical_long15heat <- ical_data_heat_coh15 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "kcal_hr") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16heat <- ical_data_heat_coh16 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "kcal_hr") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16heat <- ical_long16heat %>% 
  group_by(ID) %>%
  mutate(
    day = as.integer(as.Date(datetime) - min(as.Date(datetime))) + 1
  ) %>%
  ungroup() %>% 
  mutate(
    datetime_day = paste0(
      "day ", day, " ",
      format(datetime, "%H:%M")
    )
  )

ical_long15heat <- ical_long15heat %>%
  mutate(cohort = 15)

ical_long16heat <- ical_long16heat %>%
  mutate(cohort = 16)

common_cols <- intersect(names(ical_long15heat), names(ical_long16heat)) 

ical_long15heat <- ical_long15heat %>% select(all_of(common_cols))
ical_long16heat <- ical_long16heat %>% select(all_of(common_cols))

ical_long_allheat <- bind_rows(ical_long15heat, ical_long16heat) #this is key, here we combined

ical_long_allheat %>% 
  group_by(cohort) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have 24 animals for cohort 15 and 15 animals for cohort 16
#so in total 39 animals

ical_long_allheat <- ical_long_allheat %>%
  mutate(
    # shift so experimental day starts at 20:00
    datetime_shifted = datetime - hours(20),
    
    # experimental day (20:00 → 19:59)
    exp_day = as.integer(
      as.Date(datetime_shifted) -
        min(as.Date(datetime_shifted))
    ) + 1,
    
    # clock hour
    hour = hour(datetime),
    
    # ordering variable: 20 → 43 (20–23, then 0–19)
    hour_ordered = if_else(hour < 20, hour + 24, hour)
  ) %>%
  arrange(ID, exp_day, hour_ordered) %>%
  group_by(ID, exp_day) %>%
  mutate(kcal_hr_total = cumsum(kcal_hr)) %>%
  ungroup() %>% 
  #filter(!ID == 9406) %>% #9406 has a very weird pattern
  group_by(ID, exp_day) %>%
  mutate(
    relative_total_kcal_hr = kcal_hr_total - first(kcal_hr_total))

ical_long_allheat %>%
  filter(ID == min(ID)) %>%
  select(datetime, exp_day, hour, hour_ordered, relative_total_kcal_hr, kcal_hr_total) %>%
  arrange(datetime)

ical_long_allheat <- ical_long_allheat %>%
  mutate(
    hour_label = factor(
      hour_ordered,
      levels = 20:43,
      labels = c(20:23, 0:19)
    )
  ) %>% 
  mutate(
    lights = if_else(
      hour >= 20 | hour < 6,
      "OFF",
      "ON"
    )
  )


ggplot(
  ical_long_allheat,
  aes(
    x = hour_label,
    y = relative_total_kcal_hr,
    group = interaction(ID, exp_day)
  )
) +
  geom_line(size = 1) +
  facet_wrap(~ ID) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "relative cumulative TEE (kcal_hr)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 8)
  )

ical_long_allgroupedheat <- ical_long_allheat %>% 
  group_by(hour_label, SEX, BPA_EXPOSURE) %>% 
  summarise(
    mean_kcal_hr = mean(relative_total_kcal_hr, na.rm = TRUE),
    sem_kcal_hr  = sd(relative_total_kcal_hr, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

ggplot(
  ical_long_allgroupedheat,
  aes(
    x = hour_label,
    y = mean_kcal_hr,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_kcal_hr - sem_kcal_hr,
      ymax = mean_kcal_hr + sem_kcal_hr
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(size = 1) +
  facet_wrap(~ SEX) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "Relative cumulative TEE (kcal per h mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 10)
  )


ical_long_allheat <- ical_long_allheat %>%
  group_by(ID, exp_day) %>%
  mutate(
    relative_total_kcal_hr_lights = case_when(
      
      # Lights OFF: 20:00–05:59 → start at 0 at 20:00
      lights == "OFF" ~
        relative_total_kcal_hr - min(relative_total_kcal_hr[lights == "OFF"], na.rm = TRUE),
      
      # Lights ON: 06:00–19:59 → start at 0 at 06:00
      lights == "ON" ~
        relative_total_kcal_hr - min(relative_total_kcal_hr[lights == "ON"], na.rm = TRUE)
    )
  ) %>%
  ungroup() 

ical_long_allheat <- ical_long_allheat %>%
  mutate(
    hour_phase = case_when(
      lights == "OFF" ~ factor(
        hour,
        levels = c(20:23, 0:5)
      ),
      lights == "ON" ~ factor(
        hour,
        levels = 6:19
      )
    )
  )


ical_long_allgroupedlightsheat <- ical_long_allheat %>% 
  group_by(hour_label, SEX, BPA_EXPOSURE,lights) %>% 
  summarise(
    mean_kcal_hr = mean(relative_total_kcal_hr_lights, na.rm = TRUE),
    sem_kcal_hr  = sd(relative_total_kcal_hr_lights, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

ical_long_allgroupedlightsheat <- ical_long_allgroupedlightsheat %>%
  mutate(
    hour_phase = case_when(
      lights == "OFF" ~ factor(hour_label, levels = c(20:23, 0:5)),
      lights == "ON"  ~ factor(hour_label, levels = 6:19)
    )
  )

ggplot(
  ical_long_allgroupedlightsheat,
  aes(
    x = hour_phase,
    y = mean_kcal_hr,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE,
    group = BPA_EXPOSURE
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_kcal_hr - sem_kcal_hr,
      ymax = mean_kcal_hr + sem_kcal_hr
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(size = 1) +
  facet_grid(SEX ~ lights, scales = "free_x") +
  labs(
    x = "Time of day",
    y = "Relative TEE (kcal per hr mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 11)
  )

##RER####
ical_data_RER_coh15 <- read_csv("~/Documents/GitHub/data/data/iCal_Kotz_082425_RER.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("93", ID))) %>% 
  select(-`8/27/25 11:00`, #we want to keep just 24 hours
         -`8/27/25 12:00`,
         -`8/27/25 13:00`,
         -`8/27/25 14:00`,
         -`8/27/25 15:00`,
         -`8/27/25 16:00`,
         -`8/27/25 17:00`,
         -`8/27/25 18:00`,
         -`8/27/25 19:00`,
         -`8/28/25 20:00`,
         -`8/28/25 21:00`,
         -`8/28/25 22:00`,
         -`8/28/25 23:00`,
         -`8/29/25 0:00`,
         -`8/29/25 1:00`,
         -`8/29/25 2:00`,
         -`8/29/25 3:00`,
         -`8/29/25 4:00`,
         -`8/29/25 5:00`,
         -`8/29/25 6:00`,
         -`8/29/25 7:00`,
         -`8/29/25 8:00`,
         -`8/29/25 9:00`,
         -`8/29/25 10:00`)

ical_data_RER_coh16 <- read_csv("~/Documents/GitHub/data/data/ical analysis 11262025_RER.csv") %>% 
  rename(ID = Subject) %>% 
  mutate(ID = as.numeric(paste0("9", ID))) %>% 
  select(-`11/24/25 9:00`, #we want to keep just 24 hours
         -`11/24/25 10:00`,
         -`11/24/25 11:00`,
         -`11/24/25 12:00`,
         -`11/24/25 13:00`,
         -`11/24/25 14:00`,
         -`11/24/25 15:00`,
         -`11/24/25 16:00`,
         -`11/24/25 17:00`,
         -`11/24/25 18:00`,
         -`11/24/25 19:00`,
         -`11/25/25 20:00`,
         -`11/25/25 21:00`,
         -`11/25/25 22:00`,
         -`11/25/25 23:00`,
         -`11/26/25 0:00`,
         -`11/26/25 1:00`,
         -`11/26/25 2:00`,
         -`11/26/25 3:00`,
         -`11/26/25 4:00`,
         -`11/26/25 5:00`,
         -`11/26/25 6:00`,
         -`11/26/25 7:00`,
         -`11/26/25 8:00`,
         -`11/26/25 9:00`)

ical_long15RER <- ical_data_RER_coh15 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "RER") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16RER <- ical_data_RER_coh16 %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "RER") %>% 
  left_join(METABPA, by = "ID") %>% 
  mutate(datetime = mdy_hm(datetime_raw))

ical_long16RER <- ical_long16RER %>% 
  group_by(ID) %>%
  mutate(
    day = as.integer(as.Date(datetime) - min(as.Date(datetime))) + 1
  ) %>%
  ungroup() %>% 
  mutate(
    datetime_day = paste0(
      "day ", day, " ",
      format(datetime, "%H:%M")
    )
  )

ical_long15RER <- ical_long15RER %>%
  mutate(cohort = 15)

ical_long16RER <- ical_long16RER %>%
  mutate(cohort = 16)

common_cols <- intersect(names(ical_long15RER), names(ical_long16RER)) 

ical_long15RER <- ical_long15RER %>% select(all_of(common_cols))
ical_long16RER <- ical_long16RER %>% select(all_of(common_cols))

ical_long_allRER <- bind_rows(ical_long15RER, ical_long16RER) #this is key, here we combined

ical_long_allRER %>% 
  group_by(cohort) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have 24 animals for cohort 15 and 15 animals for cohort 16
#so in total 39 animals

ical_long_allRER <- ical_long_allRER %>%
  mutate(
    # shift so experimental day starts at 20:00
    datetime_shifted = datetime - hours(20),
    
    # experimental day (20:00 → 19:59)
    exp_day = as.integer(
      as.Date(datetime_shifted) -
        min(as.Date(datetime_shifted))
    ) + 1,
    
    # clock hour
    hour = hour(datetime),
    
    # ordering variable: 20 → 43 (20–23, then 0–19)
    hour_ordered = if_else(hour < 20, hour + 24, hour)
  ) %>%
  arrange(ID, exp_day, hour_ordered) %>%
  group_by(ID, exp_day)

ical_long_allRER <- ical_long_allRER %>%
  mutate(
    hour_label = factor(
      hour_ordered,
      levels = 20:43,
      labels = c(20:23, 0:19)
    )
  ) %>% 
  mutate(
    lights = if_else(
      hour >= 20 | hour < 6,
      "OFF",
      "ON"
    )
  )


ggplot(
  ical_long_allRER,
  aes(
    x = hour_label,
    y = RER,
    group = interaction(ID, exp_day)
  )
) +
  geom_line(size = 1) +
  facet_wrap(~ ID) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "RER (VCO2/VO2)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 8)
  )

ical_long_allgroupedRER <- ical_long_allRER %>% 
  group_by(hour_label, SEX, BPA_EXPOSURE) %>% 
  summarise(
    mean_RER = mean(RER, na.rm = TRUE),
    sem_RER  = sd(RER, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

ggplot(
  ical_long_allgroupedRER,
  aes(
    x = hour_label,
    y = mean_RER,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_RER - sem_RER,
      ymax = mean_RER + sem_RER
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(size = 1) +
  facet_wrap(~ SEX) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "RER (VCO2/VO2 mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 10)
  )

#RER ≈ 0.7: Mostly fat oxidation
#RER ≈ 0.85: Mix of fat and carbohydrate oxidation
#RER ≈ 1.0: Mostly carbohydrate oxidation
