# This script aims to explore changes in body composition and behavior in females and males CD1 offspring
#exposed to perinatal BPA 50 ug/kg and fed with HFD (D12451i, research diets) or HCD (D12451Hi, research diets)

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
library(lmerTest)
library(emmeans)
library(pracma)
library(lubridate)
library(broom)
library(stringr)
library(forcats)
library(patchwork)
library(ggpattern)
library(car)
library(broom) 

# BPA effects on dams----
dams_data <- read_csv("../data/DAMSBPAINFO.csv") %>% 
  mutate(SEX = ifelse(SEX == FALSE, "F",
                      ifelse(SEX == TRUE, "M", as.character(SEX)))) %>% 
  group_by(ID) %>% 
  mutate(F_M_sexratio = female_pups_number / male_pups_number) %>% 
  ungroup()

# Add female/male sex ratio
dams_data <- dams_data %>%
  mutate(F_M_sexratio = female_pups_number / male_pups_number)

####histogram check data distribution for litter size-----

dams_data  %>% 
  group_by(BPA_EXPOSURE,Liter_size,COHORT) %>%
  summarise(n_ID = n_distinct(ID)) 

ggplot(dams_data , aes(x = Liter_size)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic() #clearly here we have an "outlayer" but the issue is we are underpowered
#he histogram is mostly decorative rather than diagnostically powerful

ggplot(dams_data, aes(sample = Liter_size)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(dams_data$Liter_size) #W = 0.82607, p-value = 0.04034

#With such a small sample, distributional plots (histogram, QQ plot, Shapiro) are inherently unstable.
#Given the small sample size per group (n = 4–5) and the presence of an extreme value
#data were summarized using median and IQR and compared using a non-parametric Wilcoxon test

#####STATS for Liter size----

wilcox_test <- wilcox.test(Liter_size ~ BPA_EXPOSURE,
                           data = dams_data,
                           exact = FALSE)

wilcox_test

# get the median
dams_data %>%
  group_by(BPA_EXPOSURE) %>%
  summarise(
    median = median(Liter_size),
    IQR_low = quantile(Liter_size, 0.25),
    IQR_high = quantile(Liter_size, 0.75)
  )

ggplot(dams_data, aes(x = BPA_EXPOSURE, y = Liter_size)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.08, size = 3) +
  stat_summary(fun = median,
               geom = "point",
               size = 6,
               shape = 18,
               color = "red") +
  theme_classic() +
  labs(y = "Litter size",
       x = "BPA Exposure")

####histogram check data distribution for females-----

ggplot(dams_data , aes(x = female_pups_number)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(dams_data, aes(sample = female_pups_number)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(dams_data$female_pups_number) #W = 0.83547, p-value = 0.05142
#The data is not normally distributed so analyzing the mean is not correct
#I should say let see what what happen with the median

####histogram check data distribution for males-----

ggplot(dams_data , aes(x = male_pups_number)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(dams_data, aes(sample = male_pups_number)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(dams_data$male_pups_number) #W = 0.89631, p-value = 0.2314
#The data is normally distributed so analyzing the mean is  correct
#I should say let see what what happen with the median


summary_stats <- dams_data %>%
  group_by(BPA_EXPOSURE) %>%
  summarise(
    mean_Liter_size = mean(Liter_size, na.rm = TRUE),
    sem_Liter_size  = sd(Liter_size, na.rm = TRUE) / sqrt(n()),
    n_Liter_size    = n(),
    
    mean_gestation_days = mean(gestation_lenght_days, na.rm = TRUE),
    sem_gestation_days  = sd(gestation_lenght_days, na.rm = TRUE) / sqrt(n()),
    n_gestation_days    = n(),
    
    mean_female_pups = mean(female_pups_number, na.rm = TRUE),
    sem_female_pups  = sd(female_pups_number, na.rm = TRUE) / sqrt(n()),
    n_female_pups    = n(),
    
    mean_male_pups = mean(male_pups_number, na.rm = TRUE),
    sem_male_pups  = sd(male_pups_number, na.rm = TRUE) / sqrt(n()),
    n_male_pups    = n(),
    
    mean_F_M_sexratio = mean(F_M_sexratio, na.rm = TRUE),
    sem_F_M_sexratio  = sd(F_M_sexratio, na.rm = TRUE) / sqrt(n()),
    n_F_M_sexratio    = n()
  )


# First, define fill manually: YES -> grey pattern, NO -> white
summary_stats <- summary_stats %>%
  mutate(fill_color = ifelse(BPA_EXPOSURE == "YES", "pattern", "black"))

# Panel A example
pA <- ggplot(summary_stats, aes(x = BPA_EXPOSURE, y = mean_gestation_days)) +
  # draw NO and YES bars separately
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="NO"),
           aes(y = mean_gestation_days),
           fill = "black", color = "black", width = 0.6) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="YES"),
           aes(y = mean_gestation_days),
           fill = "gray", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_gestation_days - sem_gestation_days,
                    ymax = mean_gestation_days + sem_gestation_days),
                width = 0.2, color = "black", size = 0.8) +
  labs(y = "Gestational days", x = "") +
  theme_classic(base_size = 14)

pA

# Panel B: Liter size
pB <- ggplot(summary_stats, aes(x = BPA_EXPOSURE, y = mean_Liter_size)) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="NO"),
           fill = "black", color = "black", width = 0.6) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="YES"),
           fill = "gray", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_Liter_size - sem_Liter_size,
                    ymax = mean_Liter_size + sem_Liter_size),
                width = 0.2, color = "black", size = 0.8) +
  labs(y = "Liter size", x = "") +
  theme_classic(base_size = 14)

# Panel C: Female pups
pC <- ggplot(summary_stats, aes(x = BPA_EXPOSURE, y = mean_female_pups)) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="NO"),
           fill = "black", color = "black", width = 0.6) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="YES"),
           fill = "gray", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_female_pups - sem_female_pups,
                    ymax = mean_female_pups + sem_female_pups),
                width = 0.2, color = "black", size = 0.8) +
  labs(y = "Female pups", x = "") +
  theme_classic(base_size = 14)

# Panel D: Male pups
pD <- ggplot(summary_stats, aes(x = BPA_EXPOSURE, y = mean_male_pups)) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="NO"),
           fill = "black", color = "black", width = 0.6) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="YES"),
           fill = "gray", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_male_pups - sem_male_pups,
                    ymax = mean_male_pups + sem_male_pups),
                width = 0.2, color = "black", size = 0.8) +
  labs(y = "Male pups", x = "") +
  theme_classic(base_size = 14)

# Panel E: Female/male sex ratio
pE <- ggplot(summary_stats, aes(x = BPA_EXPOSURE, y = mean_F_M_sexratio)) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="NO"),
           fill = "black", color = "black", width = 0.6) +
  geom_col(data = summary_stats %>% filter(BPA_EXPOSURE=="YES"),
           fill = "gray", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_F_M_sexratio - sem_F_M_sexratio,
                    ymax = mean_F_M_sexratio + sem_F_M_sexratio),
                width = 0.2, color = "black", size = 0.8) +
  labs(y = "Female/male sex ratio", x = "BPA EXPOSURE") +
  theme_classic(base_size = 14)

##plot---- 
(pA | pB) /
  (pC | pD) /
  pE +
  plot_annotation(tag_levels = 'A')

## stats of BPA on dams----
# List of parameters
parameters <- c("gestation_lenght_days", "Liter_size",
                "female_pups_number", "male_pups_number", "F_M_sexratio")

#Function to run appropriate test comparing YES vs NO
run_comparison <- function(param) {
  
  # Check normality per group
  normality <- dams_data %>%
    group_by(BPA_EXPOSURE) %>%
    summarise(
      shapiro_p = if(length(unique(get(param))) > 1) {
        shapiro.test(get(param))$p.value
      } else {
        NA  # cannot test if all identical
      },
      .groups = "drop"
    )
  
  # Decide which test to use
  use_wilcox <- any(normality$shapiro_p < 0.05, na.rm = TRUE) | any(is.na(normality$shapiro_p))
  
  if (use_wilcox) {
    test_res <- wilcox.test(get(param) ~ BPA_EXPOSURE, data = dams_data)
    test_name <- "Wilcoxon"
  } else {
    test_res <- t.test(get(param) ~ BPA_EXPOSURE, data = dams_data)
    test_name <- "t-test"
  }
  
  data.frame(
    parameter = param,
    test = test_name,
    p_value = test_res$p.value
  )
}

# Run for all parameters and combine results
comparison_results <- lapply(parameters, run_comparison) %>% bind_rows()

comparison_results

#BODY WEIGHT (BW) ANALYSIS----
## BW separated by diet (HFD or HCD) ----
METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
filter(!grepl("-", ID)) %>%  #I eliminate from metadata all animals that were measured for NORT
  mutate(ID = as.numeric(ID)) 

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(15, 16, 18)) %>% 
  filter(!ID ==9406) %>%  #9406 has a  weird pattern in locomotion
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
    -DIET_FORMULA.y,
    -COHORT.y
  ) %>% 
  rename(
    SEX = SEX.x,
    DIET_FORMULA = DIET_FORMULA.x,
    COHORT = COHORT.x
  ) %>% 
  mutate(
    week_rel = day_rel / 7
  ) %>% 
  filter(week_rel<=19) #the last week of measurement for cohort 15 is 21, for cohort 16 is 19 so 19 is the common end

BW_data  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

BW_summary <- BW_data %>%
  group_by(week_rel,BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS females in HCD----
bw_fem_hcd <- BW_data %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_fem_hcd)   # should be > 0
table(bw_fem_hcd$BPA_EXPOSURE)

fem_bw_hcd <- lmer(
  BW ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_fem_hcd
)

anova(fem_bw_hcd)

emm_week_fem_hcd <- emmeans(
  fem_bw_hcd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_fem_hcd$week_rel)))
)

week_contrasts_fem_hcd <- contrast(
  emm_week_fem_hcd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

week_stats_fem_hcd <- as.data.frame(week_contrasts_fem_hcd) 
week_stats_fem_hcd
#so for females in HCD after week 9 BPA females are heavier than controls
#this difference in BW is extended until the end of the study (week 19)

### STATS females in HFD----

bw_fem_hfd <- BW_data %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_fem_hfd)   # should be > 0
table(bw_fem_hfd$BPA_EXPOSURE)

fem_bw_hfd <- lmer(
  BW ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_fem_hfd
)

anova(fem_bw_hfd)

emm_week_fem_hfd <- emmeans(
  fem_bw_hfd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_fem_hfd$week_rel)))
)

week_contrasts_fem_hfd <- contrast(
  emm_week_fem_hfd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

week_stats_fem_hfd <- as.data.frame(week_contrasts_fem_hfd) 
week_stats_fem_hfd 
#so for females in HFD they started with different BWs at week 0 so we can say nothing

### STATS males in HCD----

bw_male_hcd <- BW_data %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_male_hcd)   # should be > 0
table(bw_male_hcd$BPA_EXPOSURE)

male_bw_hcd <- lmer(
  BW ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_male_hcd
)

anova(male_bw_hcd)

emm_week_male_hcd <- emmeans(
  male_bw_hcd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_male_hcd$week_rel)))
)

week_contrasts_male_hcd <- contrast(
  emm_week_male_hcd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

week_stats_male_hcd <- as.data.frame(week_contrasts_male_hcd) 
week_stats_male_hcd
#so for there is no weeks in which BPA males are heavier than controls

### STATS males in HFD----

bw_male_hfd <- BW_data %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_male_hfd)   # should be > 0
table(bw_male_hfd$BPA_EXPOSURE)

male_bw_hfd <- lmer(
  BW ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_male_hfd
)

anova(male_bw_hfd)

emm_week_male_hfd <- emmeans(
  male_bw_hfd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_male_hfd$week_rel)))
)

week_contrasts_male_hfd <- contrast(
  emm_week_male_hfd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

week_stats_male_hfd <- as.data.frame(week_contrasts_male_hfd) 
week_stats_male_hfd 
#so for males they started with different BWs (BPA vs no BPA ones) but they have the same BW after 17 weeks with HFD

### overal mix model ----
model <-lmer(
  BW ~ week_rel * BPA_EXPOSURE * SEX * DIET_FORMULA +
    (1 | ID),
  data = BW_data
)
emmeans(model, pairwise ~ BPA_EXPOSURE | SEX * DIET_FORMULA * week_rel)

#for females HCD BPA changes the trajectory of BW over time, BPA animals become heavier starting around week 9–10
#for females HFD BPA animals are overall heavier, and trend toward different trajectory, BPA females are heavier from baseline, and remain heavier throughout.
#for males HCF No BPA effect on BW
#for males HFD BPA animals are consistently heavier, but slope is similar, Baseline difference present

#### STATS with baseline as a covariate (adjusts for BW starting differences)----

#Do BPA animals gain weight differently over time, independent of where they started?

baseline_bw <- BW_data %>%
  filter(week_rel == 0) %>%   # baseline timepoint
  select(ID, baseline_BW = BW)
BW_data2 <- BW_data %>%
  left_join(baseline_bw, by = "ID")
BW_data2 <- BW_data2 %>%
  select(ID, week_rel, BW, baseline_BW,SEX,BPA_EXPOSURE,DIET_FORMULA) 
BW_data2_no0 <- BW_data2 %>%
  filter(week_rel > 0)
#“Do BPA animals gain weight differently over time, controlling for baseline?
lmer(
  BW ~ week_rel * BPA_EXPOSURE * SEX * DIET_FORMULA +
    baseline_BW +
    (1 | ID),
  data = BW_data2_no0
)
emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "week_rel")
#After adjusting for baseline body weight,BPA increases body weight gain rate in females (stronger under HCD), but not in males.

#### plot A: BW over time separated by diet----

BW_summary <- BW_data %>%
  group_by(week_rel, BPA_EXPOSURE, SEX, DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

plot_bw_sex <- ggplot(
  BW_summary,
  aes(
    x = week_rel,
    y = mean_BW,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW
    ),
    alpha = 0.15,
    color = NA
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    x = "Weeks",
    y = "BW (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)#+
 # scale_x_continuous(
 #   limits = c(0, 19),
  #  breaks = seq(0, 19, by = 2)
#  )

plot_bw_sex

####plot B: slopes of BW over time (Rate of BW gain (g/week)) ----

slopes <- as.data.frame(
  emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "week_rel")
)
pairs(
  emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "week_rel")
)

plot_bw_slopes <- ggplot(
  slopes,
  aes(x = BPA_EXPOSURE,
      y = week_rel.trend,
      color = BPA_EXPOSURE,
      shape = BPA_EXPOSURE)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c("D12450Hi" = "HCD",
                       "D12451i"  = "HFD")
    )
  ) +
  labs(
    y = "Rate of BW gain (g/week)",
    x = "BPA exposure"
  ) +
  theme_classic(base_size = 14)+
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))
plot_bw_slopes

#### alternative plot: delta BW ( week 19 - week 0) ----
BW_week_delta <- BW_data %>% 
  filter(week_rel %in% c(0, 19)) %>% 
  filter(COHORT %in% c(15, 16)) %>%   # Cohort 18 is still in progress
  select(ID, SEX, DIET_FORMULA, COHORT, BPA_EXPOSURE, week_rel, BW) %>% 
  tidyr::pivot_wider(
    names_from = week_rel,
    values_from = BW,
    names_prefix = "wk_"
  ) %>% 
  mutate(
    delta_BW = wk_19 - wk_0
  )
#### alternative STATS delta BW ( week 19 - week 0) --------
BW_week_delta %>%
  group_by(SEX, DIET_FORMULA) %>%
  group_modify(~ tidy(t.test(delta_BW ~ BPA_EXPOSURE, data = .x)))

delta_bwplot <- ggplot(BW_week_delta,
       aes(x = BPA_EXPOSURE, y = delta_BW, fill = BPA_EXPOSURE)) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    x = "BPA exposure",
    y = expression(Delta*" BW (g, week 19–0)")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

delta_bwplot

#### BW at baseline (week 0) separated by diet----

BW_week0_raw <- BW_data %>% 
  filter(week_rel == 0)

BW_week0_sum <- BW_week0_raw %>% 
  group_by(BPA_EXPOSURE,SEX,DIET_FORMULA) %>% 
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
#### STATS BW at baseline (week 0) separated by diet----
bw_week0HCF <- BW_week0_raw %>% 
  filter(DIET_FORMULA=="D12450Hi")
bw_week0HFD <- BW_week0_raw %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week0HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week0HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week0HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week0HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 
#No significant differences in BW at baseline between BPA-exposed and control groups within any sex and diet condition.
#### plot C: BW at baseline (week 0)----
plot_wk_0 <- ggplot(BW_week0_sum,
                    aes(x = BPA_EXPOSURE, y = mean_BW, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_BW - sem_BW,
        ymax = mean_BW + sem_BW),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = BW_week0_raw,
    aes(x = BPA_EXPOSURE, y = BW),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "BW (g) at baseline",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_wk_0 

#Baseline body weight did not differ significantly between BPA-exposed and control groups
#within each sex and diet condition, although small numerical differences were observed in some groups.

#### BW at week 19 separated by diet----

BW_week19_raw <- BW_data %>% 
  filter(week_rel == 19)

BW_week19_sum <- BW_week19_raw %>% 
  group_by(BPA_EXPOSURE,SEX,DIET_FORMULA) %>% 
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
#### STATS BW at week 19 separated by diet----
bw_week19HCF <- BW_week19_raw %>% 
  filter(DIET_FORMULA=="D12450Hi")
bw_week19HFD <- BW_week19_raw %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week19HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week19HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week19HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_week19HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

# significant differences in BW at week 19 between BPA-exposed and control females in both OD. Nothing in males

#### plot D: BW at week 19----
plot_wk_19<- ggplot(BW_week19_sum,
                    aes(x = BPA_EXPOSURE, y = mean_BW, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_BW - sem_BW,
        ymax = mean_BW + sem_BW),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = BW_week19_raw,
    aes(x = BPA_EXPOSURE, y = BW),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "BW (g) at week 19",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_wk_19

# FIGURE 1 BW  ----

y_min <- min(c(BW_week0_raw$BW, BW_week19_raw$BW), na.rm = TRUE)
y_max <- max(c(BW_week0_raw$BW, BW_week19_raw$BW), na.rm = TRUE)

plot_bw_sex <- plot_bw_sex + labs(tag = "A")+
   scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray30"))
#delta_bwplot <- delta_bwplot + labs(tag = "B")+ 
 # scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  #scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  #theme(legend.position = "none")+
plot_wk_0<- plot_wk_0 + labs(tag = "C")+ 
 scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
theme(legend.position = "none") 
plot_wk_19 <- plot_wk_19 + labs(tag = "D")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")
plot_bw_slopes <- plot_bw_slopes + labs(tag = "B")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
 theme(legend.position = "none")

plot_wk_0 <- plot_wk_0 +
  coord_cartesian(ylim = c(y_min, y_max))

plot_wk_19 <- plot_wk_19 +
  coord_cartesian(ylim = c(y_min, y_max))

combined_plot <-combined_plot <- (plot_bw_sex | plot_bw_slopes) / (plot_wk_0 | plot_wk_19  )
combined_plot

# BODY COMPOSITION ANALYSIS----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")%>% 
  filter(!grepl("-", ID)) %>%  #I eliminate from metadata all animals that were measured for NORT
  mutate(ID = as.numeric(ID)) 

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16,18)) %>% 
  filter(!ID %in% c(9406,9354)) %>%  #9406 has a  weird pattern in locomotion so we exclude this animal from all the analyses
  #9354 last measurement were done after 26 wks with the OD 
group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, adiposity_index,COHORT,DIET_FORMULA) %>% 
  left_join(METABPA, by= "ID") %>% 
  ungroup() %>% 
  select(
    -DIET_FORMULA.y,
    -COHORT.y) %>% 
  rename(
    DIET_FORMULA = DIET_FORMULA.x,
    COHORT = COHORT.x) %>% 
  filter(!(DIET_FORMULA == "2918_teklad_Irradiated_Global_18%_Protein_Rodent_Diet"))

## Date assignation for longitudinal analysis ----

echoMRI_data %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

echoMRI_data_comparisons_collapsed <- echoMRI_data %>% 
  mutate(
    n_measurement= case_when(
      COHORT == 15 & Date == "2025-04-29" ~ "0",
      COHORT == 16 & Date == "2025-08-04" ~ "0",
      COHORT == 18 & Date %in% c("2026-01-21", "2026-01-23") ~ "0",
      COHORT == 15 & Date == "2025-05-28" ~ "4",
      COHORT == 16 & Date == "2025-09-09" ~ "4", #really this is 5 wks
     COHORT == 18 & Date == "2026-02-23" ~ "4", #really this is 4.5 wks
      COHORT == 15 & Date == "2025-07-07" ~ "10",
      COHORT == 16 & Date == "2025-10-07" ~ "10", #really 9 wks
     COHORT == 18 & Date == "2026-04-01" ~  "10",
      COHORT == 15 & Date == "2025-08-01" ~ "13",
      COHORT == 16 & Date == "2025-11-17" ~ "13", #really this is 15 wks
  #   COHORT == 18 & Date == "2026-04-29" ~ "14",  
      COHORT == 15 & Date == "2025-09-09" ~ "19",
      COHORT == 16 & Date == "2025-12-18" ~ "19",
  #   COHORT == 18 & Date == "2026-06-03" ~ "19",
      COHORT == 15 & Date == "2025-10-07" ~ "23",
      COHORT == 16 & Date == "2026-01-09" ~ "23" #really this is 22 wks ,
  #   COHORT == 18 & Date == "2026-07-08" ~ "24"
    )) %>% 
  mutate(
    n_measurement = as.numeric(
      n_measurement,
      levels = c("0", "4", "10", "13", "19", "23")
    )
  ) %>% 
  filter(!n_measurement ==23)   # for being consistent with BW data we will consider end of the study as week 19

### ADIPOSITY INDEX----
### plot A: Adiposity index over time separated by diet ----

AI_summary <- echoMRI_data_comparisons_collapsed %>%
  group_by(n_measurement, BPA_EXPOSURE, SEX, DIET_FORMULA) %>%
  summarise(
    mean_ai = mean(adiposity_index, na.rm = TRUE),
    sem_ai  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

plot_ai_sex <- ggplot(
  AI_summary,
  aes(
    x = n_measurement,
    y = mean_ai,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE,
    group = BPA_EXPOSURE
  )
) +
  geom_line(linewidth = 1) +
  #geom_point(size = 2) +
  geom_ribbon(
    aes(
      ymin = mean_ai - sem_ai,
      ymax = mean_ai + sem_ai
    ),
    alpha = 0.15,
    color = NA
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    x = "Weeks",
    y = "Adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray40")) +
  theme_classic(base_size = 14) 

plot_ai_sex

### STATS Adiposity index with baseline as a covariate (adjusts for AI starting differences)----
baseline_AI <- echoMRI_data_comparisons_collapsed %>%
  filter(n_measurement ==0) %>%   # baseline timepoint
  select(ID, baseline_ai = adiposity_index)
AI_data2 <- echoMRI_data_comparisons_collapsed %>%
  left_join(baseline_AI, by = "ID")
AI_data2 <- AI_data2 %>%
  select(ID, n_measurement, adiposity_index, baseline_ai,SEX,BPA_EXPOSURE,DIET_FORMULA) 

#Do BPA animals change adiposity index differently over time, after accounting for where they started, within each sex and diet?
AI_data2_no0 <- AI_data2 %>%
  filter(n_measurement >0)

#“Do BPA animals gain weight differently over time, controlling for baseline?
model <- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE * SEX * DIET_FORMULA +
    baseline_ai +
    (1 | ID),
  data = AI_data2_no0
)

emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")

pairs(emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement"))

###plot B: slopes of adiposity index over time (Rate of change in adiposity index (per week))----

slopes_ai <- as.data.frame(
  emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")
)

plot_ai_slopes <- ggplot(
  slopes_ai,
  aes(
    x = BPA_EXPOSURE,
    y = n_measurement.trend,
    color = BPA_EXPOSURE,
    shape = BPA_EXPOSURE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    y = "Rate of change in adiposity index (per week)",
    x = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

plot_ai_slopes

### Adiposity index at baseline (week 0) separated by diet----
AI_week0<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 0) 

AI_week0 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

AI_week0_sum <- AI_week0  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_AI_0 = mean(adiposity_index, na.rm = TRUE),
    sem_AI_0  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for adiposity index at baseline (week 0)----
ai_week0HCF <- AI_week0 %>% 
  filter(DIET_FORMULA=="D12450Hi")
ai_week0HFD <-  AI_week0 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  adiposity_index ~ BPA_EXPOSURE,
  data = ai_week0HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week0HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week0HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week0HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot C: adiposity index at baseline (week 0)----
plot_ai_0 <- ggplot(AI_week0_sum ,
                    aes(x = BPA_EXPOSURE, y = mean_AI_0, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_AI_0 - sem_AI_0,
        ymax = mean_AI_0 + sem_AI_0),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = AI_week0,
    aes(x = BPA_EXPOSURE, y = adiposity_index),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Adiposity index at baseline",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_ai_0 

### Adiposity index at week 19 separated by diet----
AI_week19<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 19) 

AI_week19 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

AI_week19_sum <- AI_week19  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_AI_19 = mean(adiposity_index, na.rm = TRUE),
    sem_AI_19  = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for adiposity index at week 19----
ai_week19HCF <- AI_week19 %>% 
  filter(DIET_FORMULA=="D12450Hi")
ai_week19HFD <-  AI_week19 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  adiposity_index ~ BPA_EXPOSURE,
  data = ai_week19HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week19HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week19HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = ai_week19HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot D: adiposity index at week 19----
plot_ai_19 <- ggplot(AI_week19_sum ,
                    aes(x = BPA_EXPOSURE, y = mean_AI_19, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_AI_19 - sem_AI_19,
        ymax = mean_AI_19 + sem_AI_19),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = AI_week19,
    aes(x = BPA_EXPOSURE, y = adiposity_index),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Adiposity index at week 19",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_ai_19

# FIGURE 2 ADIPOSITY INDEX  ----

y_min <- min(c(AI_week19$adiposity_index, AI_week19$adiposity_index), na.rm = TRUE)
y_max <- max(c(AI_week19$adiposity_index, AI_week19$adiposity_index), na.rm = TRUE)

plot_ai_sex <- plot_ai_sex + labs(tag = "A")+
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray30"))
plot_ai_0<- plot_ai_0 + labs(tag = "C")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none") 
plot_ai_19 <- plot_ai_19 + labs(tag = "D")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")
plot_ai_slopes <- plot_ai_slopes + labs(tag = "B")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")

plot_ai_0 <- plot_ai_0 +
  coord_cartesian(ylim = c(y_min, y_max))

plot_ai_19 <- plot_ai_19 +
  coord_cartesian(ylim = c(y_min, y_max))

combined_plot <-combined_plot <- (plot_ai_sex | plot_ai_slopes) / (plot_ai_0 | plot_ai_19  )
combined_plot

#ALTERNATIVE ANALYSIS ----
#delta Adiposity index, lean mass and fat mass
## data preparation

echoMRI_data_comparisons_collapsed_delta <- echoMRI_data_comparisons_collapsed %>% 
  filter(n_measurement %in% c(0, 19)) %>% 
  filter(COHORT %in% c(15, 16)) %>%   # Cohort 18 is still in progress
  select(ID, SEX, DIET_FORMULA, COHORT, BPA_EXPOSURE, n_measurement, Fat,Lean,adiposity_index) %>% 
  tidyr::pivot_wider(
    names_from = n_measurement,
    values_from = c(Fat, Lean, adiposity_index),
    names_prefix = "wk_"
  ) %>% 
  mutate(
    delta_Fat = Fat_wk_19 - Fat_wk_0,
    delta_Lean = Lean_wk_19 - Lean_wk_0,
    delta_adiposity_index = adiposity_index_wk_19 - adiposity_index_wk_0
  )

#### plot B: delta adiposity index ( week 19 - week 0) 
delta_adiposity_indexplot <- ggplot(echoMRI_data_comparisons_collapsed_delta,
                                    aes(x = BPA_EXPOSURE, y = delta_adiposity_index, fill = BPA_EXPOSURE)) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    x = "BPA exposure",
    y = expression(Delta*" adiposity index (week 19–0)")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

delta_adiposity_indexplot

### STATS delta adiposity index ( week 19 - week 0) 
echoMRI_data_comparisons_collapsed_delta %>%
  group_by(SEX, DIET_FORMULA) %>%
  group_modify(~ tidy(t.test(delta_adiposity_index ~ BPA_EXPOSURE, data = .x)))

##Figure alternative ADIPOSITY INDEX
plot_ai_sex  <- plot_ai_sex  + labs(tag = "A")+
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray40"))
delta_adiposity_indexplot <- delta_adiposity_indexplot + labs(tag = "B") + 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) 
  plot_ai_slopes<- plot_ai_slopes + labs(tag = "C")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")

combined_plot <- plot_ai_sex | delta_adiposity_indexplot| plot_ai_slopes
combined_plot

#### STATS delta Fat ( week 19 - week 0) 
echoMRI_data_comparisons_collapsed_delta %>%
  group_by(SEX, DIET_FORMULA) %>%
  group_modify(~ tidy(t.test(delta_Fat ~ BPA_EXPOSURE, data = .x)))
##### plot alternative: delta fat mass ( week 19 - week 0) 
delta_fatplot <- ggplot(echoMRI_data_comparisons_collapsed_delta,
                        aes(x = BPA_EXPOSURE, y = delta_Fat, fill = BPA_EXPOSURE)) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    x = "BPA exposure",
    y = expression(Delta*" Fat mass (g, week 19–0)")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

delta_fatplot

#### STATS delta lean ( week 19 - week 0) 
echoMRI_data_comparisons_collapsed_delta %>%
  group_by(SEX, DIET_FORMULA) %>%
  group_modify(~ tidy(t.test(delta_Lean ~ BPA_EXPOSURE, data = .x)))

#BPA increases lean mass gain in females under HFD
# BPA does not significantly affect fat mass gain
# BPA effects in males are minimal or absent
# Δ analyses show trends consistent with longitudinal models but are less sensitive


##### plot alternative: delta lean mass ( week 19 - week 0) 
delta_leanplot <- ggplot(echoMRI_data_comparisons_collapsed_delta,
                        aes(x = BPA_EXPOSURE, y = delta_Lean, fill = BPA_EXPOSURE)) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    x = "BPA exposure",
    y = expression(Delta*" Lean mass (g, week 19–0)")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

delta_leanplot

### FAT MASS----
## plot A: fat mass separated by diet over time----

fat_summary <- echoMRI_data_comparisons_collapsed %>%
  group_by(n_measurement, BPA_EXPOSURE, SEX, DIET_FORMULA) %>%
  summarise(
    mean_fat = mean(Fat, na.rm = TRUE),
    sem_fat  = sd(Fat, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

plot_fat_sex <- ggplot(
  fat_summary,
  aes(
    x = n_measurement,
    y = mean_fat,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE,
    group = BPA_EXPOSURE
  )
) +
  geom_line(linewidth = 1) +
  #geom_point(size = 2) +
  geom_ribbon(
    aes(
      ymin = mean_fat - sem_fat,
      ymax = mean_fat + sem_fat
    ),
    alpha = 0.15,
    color = NA
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    x = "Weeks",
    y = "Fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "gray40")) +
  theme_classic(base_size = 14)

plot_fat_sex

## STATS with baseline as a covariate (adjusts for fat mass starting differences)----

baseline_fat <- echoMRI_data_comparisons_collapsed %>%
  filter(n_measurement ==0) %>%   # baseline timepoint
  select(ID, baseline_fat = Fat)
fat_data2 <- echoMRI_data_comparisons_collapsed %>%
  left_join(baseline_fat, by = "ID")
fat_data2 <- fat_data2 %>%
  select(ID, n_measurement, Fat, baseline_fat,SEX,BPA_EXPOSURE,DIET_FORMULA) 

#Do BPA animals change fat mass differently over time, after accounting for where they started, within each sex and diet?
fat_data2_no0 <- fat_data2 %>%
  filter(n_measurement >0)

#“Do BPA animals gain fat mass differently over time, controlling for baseline?
model <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE * SEX * DIET_FORMULA +
    baseline_fat +
    (1 | ID),
  data = fat_data2_no0
)

emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")

pairs(emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement"))

####plot B: slopes of fat mass over time (Rate of fat mass gain (g/week)) ----

slopes_fat <- as.data.frame(
  emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")
)

plot_fat_slopes <- ggplot(
  slopes_fat,
  aes(
    x = BPA_EXPOSURE,
    y = n_measurement.trend,
    color = BPA_EXPOSURE,
    shape = BPA_EXPOSURE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    y = "Rate of fat mass gain (g/week)",
    x = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

plot_fat_slopes

### Fat mass at baseline (week 0) separated by diet----
fat_week0<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 0) 

fat_week0 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

fat_week0_sum <- fat_week0  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_fat_0 = mean(Fat, na.rm = TRUE),
    sem_fat_0  = sd(Fat, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for fat mass at baseline (week 0)----
fat_week0HCF <- fat_week0 %>% 
  filter(DIET_FORMULA=="D12450Hi")
fat_week0HFD <-  fat_week0 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  Fat ~ BPA_EXPOSURE,
  data = fat_week0HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  Fat  ~ BPA_EXPOSURE,
  data =fat_week0HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  Fat  ~ BPA_EXPOSURE,
  data = fat_week0HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  adiposity_index  ~ BPA_EXPOSURE,
  data = fat_week0HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot C: adiposity index at baseline (week 0)----
plot_fat_0 <- ggplot(fat_week0_sum ,
                    aes(x = BPA_EXPOSURE, y = mean_fat_0, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_fat_0 - sem_fat_0,
        ymax = mean_fat_0 + sem_fat_0),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = fat_week0,
    aes(x = BPA_EXPOSURE, y = Fat),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Fat mass (g) at baseline",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_fat_0 

### Fat mass at week 19 separated by diet----
fat_week19<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 19) 

fat_week19 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

fat_week19_sum <- fat_week19  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_fat_19 = mean(Fat, na.rm = TRUE),
    sem_fat_19  = sd(Fat, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for fat mass at week 19----
fat_week19HCF <- fat_week19 %>% 
  filter(DIET_FORMULA=="D12450Hi")
fat_week19HFD <-  fat_week19 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  Fat ~ BPA_EXPOSURE,
  data = fat_week19HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  Fat  ~ BPA_EXPOSURE,
  data = fat_week19HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  Fat  ~ BPA_EXPOSURE,
  data = fat_week19HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  Fat  ~ BPA_EXPOSURE,
  data = fat_week19HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot D: fat mass at week 19----
plot_fat_19 <- ggplot(fat_week19_sum ,
                     aes(x = BPA_EXPOSURE, y = mean_fat_19, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_fat_19 - sem_fat_19,
        ymax = mean_fat_19 + sem_fat_19),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = fat_week19,
    aes(x = BPA_EXPOSURE, y = Fat),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Fat mass (g) at week 19",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_fat_19

# FIGURE 3 FAT MASS  ----

y_min <- min(c(fat_week19$Fat, fat_week19$Fat), na.rm = TRUE)
y_max <- max(c(fat_week19$Fat, fat_week19$Fat), na.rm = TRUE)

plot_fat_sex <- plot_fat_sex + labs(tag = "A")+
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray30"))
plot_fat_0<- plot_fat_0 + labs(tag = "C")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none") 
plot_fat_19 <- plot_fat_19 + labs(tag = "D")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")
plot_fat_slopes <- plot_fat_slopes + labs(tag = "B")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")

plot_fat_0 <- plot_fat_0 +
  coord_cartesian(ylim = c(y_min, y_max))

plot_fat_19 <- plot_fat_19 +
  coord_cartesian(ylim = c(y_min, y_max))

combined_plot <-combined_plot <- (plot_fat_sex | plot_fat_slopes) / (plot_fat_0 | plot_fat_19  )
combined_plot

### LEAN MASS----
## plot A: lean mass separated by diet over time----
lean_summary <- echoMRI_data_comparisons_collapsed %>%
  group_by(n_measurement, BPA_EXPOSURE, SEX, DIET_FORMULA) %>%
  summarise(
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean  = sd(Lean, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

plot_lean_sex <- ggplot(
  lean_summary,
  aes(
    x = n_measurement,
    y = mean_lean,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE,
    group = BPA_EXPOSURE
  )
) +
  geom_line(linewidth = 1) +
  #geom_point(size = 2) +
  geom_ribbon(
    aes(
      ymin = mean_lean - sem_lean,
      ymax = mean_lean + sem_lean
    ),
    alpha = 0.15,
    color = NA
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    x = "Weeks",
    y = "Lean mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "gray40")) +
  theme_classic(base_size = 14)

plot_lean_sex

## STATS with baseline as a covariate (adjusts for lean mass starting differences)----

baseline_lean <- echoMRI_data_comparisons_collapsed %>%
  filter(n_measurement ==0) %>%   # baseline timepoint
  select(ID, baseline_lean = Lean)
lean_data2 <- echoMRI_data_comparisons_collapsed %>%
  left_join(baseline_lean, by = "ID")
lean_data2 <- lean_data2 %>%
  select(ID, n_measurement, Lean, baseline_lean,SEX,BPA_EXPOSURE,DIET_FORMULA) 

#Do BPA animals change lean mass differently over time, after accounting for where they started, within each sex and diet?
lean_data2_no0 <- lean_data2 %>%
  filter(n_measurement >0)

#“Do BPA animals lean mass differently over time, controlling for baseline?
model <- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE * SEX * DIET_FORMULA +
    baseline_lean +
    (1 | ID),
  data = lean_data2_no0
)

emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")
pairs(emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement"))

####plot B: slopes of lean mass over time (Rate of lean mass gain (g/week)) ----

slopes_lean <- as.data.frame(
  emtrends(model, ~ BPA_EXPOSURE | SEX * DIET_FORMULA, var = "n_measurement")
)
plot_lean_slopes <- ggplot(
  slopes_lean,
  aes(
    x = BPA_EXPOSURE,
    y = n_measurement.trend,
    color = BPA_EXPOSURE,
    shape = BPA_EXPOSURE
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    y = "Rate of lean mass gain (g/week)",
    x = "BPA exposure"
  ) +
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
plot_lean_slopes

### lean mass at baseline (week 0) separated by diet----
lean_week0<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 0) 

lean_week0 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

lean_week0_sum <- lean_week0  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_lean_0 = mean(Lean, na.rm = TRUE),
    sem_lean_0  = sd(Lean, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for lean mass at baseline (week 0)----
lean_week0HCF <- lean_week0 %>% 
  filter(DIET_FORMULA=="D12450Hi")
lean_week0HFD <-  lean_week0 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  Lean ~ BPA_EXPOSURE,
  data = lean_week0HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  Lean  ~ BPA_EXPOSURE,
  data =lean_week0HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
 Lean  ~ BPA_EXPOSURE,
  data = lean_week0HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
  Lean  ~ BPA_EXPOSURE,
  data = lean_week0HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot C: lean mass at baseline (week 0)----
plot_lean_0 <- ggplot(lean_week0_sum ,
                     aes(x = BPA_EXPOSURE, y = mean_lean_0, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_lean_0 - sem_lean_0,
        ymax = mean_lean_0 + sem_lean_0),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = lean_week0,
    aes(x = BPA_EXPOSURE, y = Lean),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Lean mass (g) at baseline",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_lean_0 

### lean mass at week 19 separated by diet----
lean_week19<- echoMRI_data_comparisons_collapsed  %>% 
  filter(n_measurement == 19) 

lean_week19 %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

lean_week19_sum <- lean_week19  %>% 
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>% 
  summarise(
    mean_lean_19 = mean(Lean, na.rm = TRUE),
    sem_lean_19  = sd(Lean, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### STATS for lean mass at week 19----
lean_week19HCF <- lean_week19 %>% 
  filter(DIET_FORMULA=="D12450Hi")
lean_week19HFD <-  lean_week19 %>% 
  filter(DIET_FORMULA=="D12451i")
# Females in HCD
t_femaleHCD <- t.test(
  Lean ~ BPA_EXPOSURE,
  data = lean_week19HCF  %>% filter(SEX == "F"))
# Females in HFD
t_femaleHFD <- t.test(
  Lean  ~ BPA_EXPOSURE,
  data = lean_week19HFD  %>% filter(SEX == "F")
)
# Males in HCD
t_maleHCD <- t.test(
  Lean  ~ BPA_EXPOSURE,
  data = lean_week19HCF  %>% filter(SEX == "M"))
# Males in HFD
t_maleHFD <- t.test(
 Lean ~ BPA_EXPOSURE,
  data = lean_week19HFD  %>% filter(SEX == "M")
)
t_femaleHCD 
t_femaleHFD 
t_maleHCD
t_maleHFD 

### plot D: lean mass at week 19----
plot_lean_19 <- ggplot(lean_week19_sum ,
                      aes(x = BPA_EXPOSURE, y = mean_lean_19, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_lean_19 - sem_lean_19,
        ymax = mean_lean_19 + sem_lean_19),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = lean_week19,
    aes(x = BPA_EXPOSURE, y = Lean),
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  labs(
    y = "Lean mass (g) at week 19",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("NO" = "gray80", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "black"))
plot_lean_19

# FIGURE 4 LEAN MASS  ----

y_min <- min(c(lean_week0$Lean, lean_week0$Lean), na.rm = TRUE)
y_max <- max(c(lean_week19$Lean, lean_week19$Lean), na.rm = TRUE)

plot_lean_sex <- plot_lean_sex + labs(tag = "A")+
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_fill_manual(values = c("NO" = "gray50", "YES" = "gray30"))
plot_lean_0<- plot_lean_0 + labs(tag = "C")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none") 
plot_lean_19 <- plot_lean_19 + labs(tag = "D")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")
plot_lean_slopes <- plot_lean_slopes + labs(tag = "B")+ 
  scale_color_manual(values = c("NO" = "gray50", "YES" = "black")) +
  scale_shape_manual(values = c("NO" = 16, "YES" = 17))+
  theme(legend.position = "none")

plot_lean_0 <- plot_lean_0 +
  coord_cartesian(ylim = c(y_min, y_max))

plot_lean_19 <- plot_lean_19 +
  coord_cartesian(ylim = c(y_min, y_max))

combined_plot <-combined_plot <- (plot_lean_sex | plot_lean_slopes) / (plot_lean_0 | plot_lean_19  )
combined_plot

#LENGTH ANALYSIS----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  filter(!grepl("-", ID)) %>%  #I eliminate from metadata all animals that were measured for NORT
  mutate(ID = as.numeric(ID))

length_data <- read_csv("../data/METABPA_LENGTH.csv") %>% 
  left_join(METABPA, by= "ID") %>% 
  group_by(ID) %>% 
  mutate(
    DATE = mdy(DATE),
    DOB  = mdy(DOB),
    age_weeks = as.numeric(difftime(DATE, DOB, units = "weeks")))%>%
  mutate(age_week_round = round(age_weeks)) %>% 
  ungroup()

length_data %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

####histogram check data distribution-----

ggplot(length_data, aes(x = LENGTH_CM)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(length_data, aes(sample = LENGTH_CM)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(length_data$LENGTH_CM) #The data are normally distributed so analysing the mean is correct

length_summary <- length_data  %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_length = mean(LENGTH_CM, na.rm = TRUE),
    sem_length  = sd(LENGTH_CM, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n = n_distinct(ID),
    .groups = "drop"
  ) 

##STATS----
# Female groups
female_data <- length_data %>% filter(SEX == "F")
shapiro.test(female_data$LENGTH_CM[female_data$BPA_EXPOSURE == "YES"])
shapiro.test(female_data$LENGTH_CM[female_data$BPA_EXPOSURE == "NO"])

# Male groups
male_data <- length_data %>% filter(SEX == "M")
shapiro.test(male_data$LENGTH_CM[male_data$BPA_EXPOSURE == "YES"])
shapiro.test(male_data$LENGTH_CM[male_data$BPA_EXPOSURE == "NO"])

# Female
leveneTest(LENGTH_CM ~ BPA_EXPOSURE, data = female_data)

# Male
leveneTest(LENGTH_CM ~ BPA_EXPOSURE, data = male_data)

# Welch t-test for females
t.test(LENGTH_CM ~ BPA_EXPOSURE, data = female_data, var.equal = FALSE)
# Standard t-test for males
t.test(LENGTH_CM ~ BPA_EXPOSURE, data = male_data, var.equal = TRUE)

### Length plot----
length_plot <- ggplot(
  length_data,
  aes(x = BPA_EXPOSURE, y = LENGTH_CM, fill = BPA_EXPOSURE)
) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ SEX) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    x = "BPA exposure",
    y = "Length (cm)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

length_plot

# FOOD INTAKE ANALYSIS----

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(15,16)) %>% #COHORT 18 IS STILL IN PROGRESS
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) %>% 
  left_join(METABPA, by= "ID") %>% 
  select(
    -SEX.y,
    -DIET_FORMULA.y,
    -DIET_FORMULA.x,
    -COHORT.y
  ) %>% 
  rename(
    SEX = SEX.x,
    COHORT = COHORT.x
  ) %>% 
  mutate(
    FIcumulative = cumsum(corrected_intake_kcal)) %>% 
  ungroup() %>% 
  filter(!ID ==9406) %>% #9406 has a  weird pattern in locomotion
  mutate(
    week_rel = day_rel / 7
  ) %>% 
  filter(week_rel<=19) #the last week of measurement for cohort 15 is 21, for cohort 16 is 19 so 19 is the common end

FI_data   %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

####histogram check data distribution-----

ggplot(FI_data, aes(x = FIcumulative)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(FI_data, aes(sample = FIcumulative)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(FI_data$FIcumulative) #The data are normally distributed so analyzing the mean is correct

####FI separated by diet----

FI_data  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID))

FI_final <- FI_data %>%
  arrange(ID, DATE) %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup()

FI_plotA <- ggplot(
  FI_final,
  aes(x = BPA_EXPOSURE, y = FIcumulative, fill = BPA_EXPOSURE)
) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_fill_manual(values = c("NO" = "gray80", "YES" = "black")) +
  labs(
    y = "Cumulative food intake (kcal)",
    x = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

FI_plotA

##### Do males eat more than females within each diet, collapsing BPA?

FI_final %>%
  group_by(DIET_FORMULA) %>%
  group_modify(~ tidy(t.test(FIcumulative ~ SEX, data = .x)))

#Males do NOT eat more than females in your dataset
#Differences in BW are likely not driven by food intake

FI_plotB <-ggplot(
  FI_final,
  aes(x = SEX, y = FIcumulative, fill = SEX)
) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.6,
    color = "black"
  ) +
  
  # SEM
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8
  ) +
  
  # Individual animals
  geom_jitter(
    width = 0.15,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  
  # Separate by diet
  facet_wrap(
    ~ DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  
  # Colors
  scale_fill_manual(values = c("F" = "gray70", "M" = "black")) +
  
  # Labels
  labs(
    x = "Sex",
    y = "Cumulative food intake (kcal)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

FI_plotB 

# Figure 5 (Cumulative food intake at week 19) ----
FI_plotA <- FI_plotA + labs(tag = "A")
FI_plotB <- FI_plotB + labs(tag = "B")

combined_plot <- FI_plotA | FI_plotB
combined_plot

# OGTT----

#OGTT WAS DONE FOR COHORT 15 AND 16 WHEN ALL ANIMALS WERE 27 WEEK OLD

#cohort 15: D.O.B 4/10/2025 - 10/13/2025 = 26.6 wks
#cohort 16: D.O.B 7/13/2025 - 1/19/2026 = 27.1 wks
#cohort 18: D.O.B 12/31/2025 and 1/2/2026 - 07/03/2026 = 26.4 wks

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  mutate(ID= as.numeric(ID))

OGTT <- read_csv("~/Documents/GitHub/data/data/OGTT_CD1.csv") %>% 
  mutate(DATE = mdy(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID)%>% 
  left_join(METABPA, by= "ID") %>% 
  filter(!ID ==9406)  #9406 has a  weird pattern in locomotion
  
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

####histogram check data distribution-----

ggplot(auc_df, aes(x = AUC)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(auc_df, aes(sample = AUC)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic()

shapiro.test(auc_df$AUC) #The data are normally distributed so analyzing the mean is correct

#mean plot

ogtta <- ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
 facet_wrap(~ SEX*DIET_FORMULA) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()
ogtta

#median
# Compute median and 25/75th percentiles for each group
summary_df <- auc_df %>%
  group_by(BPA_EXPOSURE, SEX,DIET_FORMULA) %>%
  summarise(
    med = median(AUC),
    ymin = quantile(AUC, 0.25),
    ymax = quantile(AUC, 0.75)
  )

# Plot median
ogtta <- ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_point(data = summary_df, aes(x = BPA_EXPOSURE, y = med), size = 3) +
  geom_errorbar(data = summary_df, aes(x = BPA_EXPOSURE, y = med, ymin = ymin, ymax = ymax), width = 0.2) +
  facet_wrap(~ SEX*DIET_FORMULA) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()
ogtta

##### STATS for females COLLAPSED FOR DIET----

ogtt_fem<- auc_df %>%
  filter(
    SEX == "F",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_fem)
table(ogtt_fem$BPA_EXPOSURE)

by(
  ogtt_fem$AUC,
  ogtt_fem$BPA_EXPOSURE,
  shapiro.test
)

leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_fem)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_fem
)


##### STATS for males----

ogtt_m<- auc_df %>%
  filter(
    SEX == "M",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_m)
table(ogtt_m$BPA_EXPOSURE)

by(
  ogtt_m$AUC,
  ogtt_m$BPA_EXPOSURE,
  shapiro.test
)

leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_m)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_m
)

####OGTT collapsed by BPA exposure----

ogtt_long %>% 
  group_by(SEX) %>%
  summarise(n_ID = n_distinct(ID)) 


auc_df <- ogtt_long %>% 
  arrange(ID, time_min) %>% 
  group_by(ID,SEX) %>% 
  summarise(
    AUC = trapz(time_min, glucose),
    .groups = "drop"
  )

ogttb <- ggplot(auc_df, aes(x = SEX, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(
    x = "SEX",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()
ogttb 

##### STATS females vs males----

ogtt<- auc_df 
nrow(ogtt)
table(ogtt$SEX)

by(
  ogtt$AUC,
  ogtt$SEX,
  shapiro.test
)

leveneTest(AUC ~ SEX, data = ogtt)

t.test(
  AUC ~ SEX,
  data = ogtt
)

#this means that females manage better the glucose captation than males independent of BPA exposure
  
####OGTT separated by diet----

ogtt_long %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 


auc_df <- ogtt_long %>% 
  arrange(ID, time_min) %>% 
  group_by(ID, BPA_EXPOSURE,SEX,DIET_FORMULA) %>% 
  summarise(
    AUC = trapz(time_min, glucose),
    .groups = "drop"
  )

ogttc <- ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ SEX*DIET_FORMULA,
             labeller = labeller(
               DIET_FORMULA = c(
                 "D12450Hi" = "HCD",
                 "D12451i"  = "HFD" ))) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()
ogttc

#####STATS for females in HCD----

ogtt_fem_hcd <- auc_df %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_fem_hcd)   # should be > 0
table(ogtt_fem_hcd$BPA_EXPOSURE)

by(
  ogtt_fem_hcd$AUC,
  ogtt_fem_hcd$BPA_EXPOSURE,
  shapiro.test
)

leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_fem_hcd)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_fem_hcd
)

##### STATS for females in HFD----

ogtt_fem_hfd <- auc_df %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_fem_hfd)   # should be > 0
table(ogtt_fem_hfd$BPA_EXPOSURE)

by(
  ogtt_fem_hfd$AUC,
  ogtt_fem_hfd$BPA_EXPOSURE,
  shapiro.test
)


leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_fem_hfd)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_fem_hfd
)

#####STATS for males in HCD----

ogtt_m_hcd <- auc_df %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_m_hcd)   # should be > 0
table(ogtt_m_hcd$BPA_EXPOSURE)

by(
  ogtt_m_hcd$AUC,
  ogtt_m_hcd$BPA_EXPOSURE,
  shapiro.test
)

leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_m_hcd)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_m_hcd
)

##### STATS for males in HFD----

ogtt_m_hfd <- auc_df %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(ogtt_m_hfd)   # should be > 0
table(ogtt_m_hfd$BPA_EXPOSURE)

by(
  ogtt_m_hfd$AUC,
  ogtt_m_hfd$BPA_EXPOSURE,
  shapiro.test
)


leveneTest(AUC ~ BPA_EXPOSURE, data = ogtt_m_hfd)

t.test(
  AUC ~ BPA_EXPOSURE,
  data = ogtt_m_hfd
)

# supplementary figure 7 (ogtta + ogttb + ogttc)----

sf7b <- ogttc + labs(tag = "A")
sf7c <- ogttb + labs(tag = "B")

combined_plot_ogtt <- sf7b | sf7c 
combined_plot_ogtt

# INDIRECT CALORIMETRY / COLUMBUS DATA ANALYSIS ----

#cohort 15: D.O.B 4/10/2025 - 08/24/2025 = 19.4 wks
#cohort 16: D.O.B 7/13/2025 - 11/24/2025 = 19.1 wks
#cohort 18: D.O.B 12/31/2025 and 1/2/2026 - 05/16/2026 = 19.4 wks

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  mutate(ID= as.numeric(ID))

##counts /locomotion ----
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

ical_long_all <- bind_rows(ical_long15, ical_long16) %>% #this is key, here we combined
  #filter(!ID %in% c(9367, 9366, 9404, 9363,9406)) #these animals are responsible for a skew behavior of the normal curve in locomotion data
  filter(!ID ==9406)  #9406 has a  weird pattern in locomotion
  
ical_long_all %>% 
filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 

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
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
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

#### individual locomotion data----

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
filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(hour_label, SEX, BPA_EXPOSURE,DIET_FORMULA) %>% 
  summarise(
    mean_counts = mean(relative_total_count, na.rm = TRUE),
    sem_counts  = sd(relative_total_count, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

dark_phase <- data.frame(
  xmin = 1,
  xmax = 11,
  ymin = -Inf,
  ymax = Inf
)

#### cumulative locomotion separated by diet----

ggplot(
  ical_long_allgrouped,
  aes(
    x = as.numeric(hour_label),
    y = mean_counts,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  # DARK PHASE BACKGROUND
  geom_rect(
    data = dark_phase,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.5
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
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_x_continuous(
    breaks = 1:24,
    labels = levels(ical_long_allgrouped$hour_label)
  ) +
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

# Summarize relative_total_count at hour == 19 including SEX and BPA_EXPOSURE so this means locomotion over 24h
relative_count_at_19 <- ical_long_all %>%
  filter(hour == 19) %>%                # keep only rows at hour 19
  group_by(ID) %>%
  summarise(
    relative_total_count_19 = first(relative_total_count),  # value at 19:00
    cohort = first(cohort),
    SEX = first(SEX),
    BPA_EXPOSURE = first(BPA_EXPOSURE),
    .groups = "drop"
  ) %>% 
 # filter(!ID %in% c(9367, 9366, 9404, 9363,9406))
filter(!ID ==9406)   #9406 has a  weird pattern in locomotion
  
  relative_count_at_19 %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 

####histogram check data distribution-----

ggplot(relative_count_at_19 , aes(x = relative_total_count_19 )) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic()

ggplot(relative_count_at_19, aes(sample = relative_total_count_19)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(SEX ~ BPA_EXPOSURE) +
  theme_classic()

# STATS locomotion Run unpaired t-test within each sex ----

# Check variance within each sex
variance_check <- relative_count_at_19 %>%
  group_by(SEX) %>%
  summarise(
    var_test = list(
      var.test(
        relative_total_count_19 ~ BPA_EXPOSURE, 
        data = cur_data()
      )
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    f_statistic = var_test$statistic,
    df1 = var_test$parameter[1],
    df2 = var_test$parameter[2],
    p_value = var_test$p.value
  ) %>%
  select(SEX, f_statistic, df1, df2, p_value)

variance_check #so because the assumption of equal variance is violated, so Welch’s t-test is justified.

t_test_results <- relative_count_at_19 %>%
  group_by(SEX) %>%
  summarise(
    t_test = list(t.test(
      relative_total_count_19 ~ BPA_EXPOSURE,  # formula: group comparison
      data = cur_data(),
      var.equal = FALSE                        # Welch's t-test (does not assume equal variance)
    )),
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),
    df = map_dbl(t_test, ~ .x$parameter),
    p_value = map_dbl(t_test, ~ .x$p.value),
    mean_NO = map_dbl(t_test, ~ .x$estimate[1]),   # mean for BPA_EXPOSURE == NO
    mean_YES = map_dbl(t_test, ~ .x$estimate[2])   # mean for BPA_EXPOSURE == YES
  ) %>%
  select(SEX, mean_NO, mean_YES, t_statistic, df, p_value)

t_test_results

# Compute mean and SEM per SEX × BPA_EXPOSURE
summary_19 <- relative_count_at_19 %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_count = mean(relative_total_count_19),
    sem_count = sd(relative_total_count_19) / sqrt(n()),
    .groups = "drop"
  )


# T-test results
t_test_summary <- relative_count_at_19 %>%
  group_by(SEX) %>%
  do(tidy(t.test(relative_total_count_19 ~ BPA_EXPOSURE, data = .))) %>%
  ungroup() %>%
  mutate(
    # Positions to display p-values on plot
    y_position = max(relative_count_at_19$relative_total_count_19) * 1.05
  ) %>%
  select(SEX, estimate1, estimate2, statistic, p.value, y_position)

t_test_summary

## supplementary figure 8A (SF8A) locomotion over 24h ----
sf8a <- ggplot(summary_19, aes(x = BPA_EXPOSURE, y = mean_count, fill = BPA_EXPOSURE)) +
  geom_col(alpha = 0.7, width = 0.6) +  # bars
  geom_errorbar(aes(ymin = mean_count - sem_count, ymax = mean_count + sem_count),
                width = 0.2) +           # SEM
  geom_jitter(
    data = relative_count_at_19,
    aes(x = BPA_EXPOSURE, y = relative_total_count_19, color = BPA_EXPOSURE),
    width = 0.15,
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ SEX) +
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  # Add p-values from t-test
  geom_text(
    data = t_test_summary,
    aes(x = 1.5,  # midpoint between NO (1) and YES (2)
        y = y_position,
        label = paste0("p = ", signif(p.value, 3))),
    inherit.aes = FALSE
  ) +
  labs(
    x = "BPA Exposure",
    y = "Counts over 24h",
    fill = "BPA Exposure",
    color = "BPA Exposure"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )
sf8a


# Summarize relative_total_count at hour == 19 including SEX and BPA_EXPOSURE and DIET so this means locomotion over 24h
relative_count_at_19_diet <- ical_long_all %>%
  filter(hour == 19) %>%
  group_by(ID) %>%
  summarise(
    relative_total_count_19 = first(relative_total_count),
    cohort = first(cohort),
    SEX = first(SEX),
    BPA_EXPOSURE = first(BPA_EXPOSURE),
    DIET_FORMULA = first(DIET_FORMULA),
    .groups = "drop"
  ) %>%
  mutate(
    DIET_FORMULA = dplyr::recode(DIET_FORMULA,
                                 "D12450Hi" = "HCD",
                                 "D12451i"  = "HFD")
  )


# females
model_f <- lm(relative_total_count_19 ~ BPA_EXPOSURE * DIET_FORMULA,
              data = relative_count_at_19_diet %>% filter(SEX == "F"))

anova(model_f)
summary(model_f)

# males
model_m <- lm(relative_total_count_19 ~ BPA_EXPOSURE * DIET_FORMULA,
              data = relative_count_at_19_diet %>% filter(SEX == "M"))

anova(model_m)
summary(model_m)


# females
emmeans(model_f, pairwise ~ BPA_EXPOSURE | DIET_FORMULA)

# males
emmeans(model_m, pairwise ~ BPA_EXPOSURE | DIET_FORMULA)

## supplementary figure 8B (SF8B) locomotion over 24h separated by diet ----
# Mean and SEM per SEX × BPA × DIET
summary_19_diet <- relative_count_at_19_diet %>%
  group_by(SEX, DIET_FORMULA, BPA_EXPOSURE) %>%
  summarise(
    mean_count = mean(relative_total_count_19),
    sem_count = sd(relative_total_count_19) / sqrt(n()),
    .groups = "drop"
  )

sf8b <- ggplot(summary_19_diet,
                    aes(x = DIET_FORMULA,
                        y = mean_count,
                        fill = BPA_EXPOSURE)) +
  
  geom_col(position = position_dodge(width = 0.7),
           alpha = 0.7,
           width = 0.6) +
  
  geom_errorbar(
    aes(ymin = mean_count - sem_count,
        ymax = mean_count + sem_count),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  
  geom_jitter(
    data = relative_count_at_19_diet,
    aes(x = DIET_FORMULA,
        y = relative_total_count_19,
        color = BPA_EXPOSURE),
    position = position_jitterdodge(jitter.width = 0.15,
                                    dodge.width = 0.7),
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ SEX) +
  
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  
  labs(
    x = "Diet",
    y = "Counts over 24h",
    fill = "BPA Exposure",
    color = "BPA Exposure"
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

sf8b

## locomotion analysis separated by light period ----

ical_phase <- ical_long_all %>%
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  mutate(
    phase = case_when(
      hour >= 20 | hour < 6  ~ "dark",   # 20:00–05:59
      hour >= 6  & hour < 20 ~ "light"   # 06:00–19:59
    )
  )

ical_phase <- ical_phase %>%
  mutate(
    phase_hour_order = case_when(
      phase == "dark"  & hour >= 20 ~ hour,
      phase == "dark"  & hour < 6   ~ hour + 24,
      phase == "light"              ~ hour
    )
  ) %>%
  arrange(ID, exp_day, phase, phase_hour_order)

ical_phase <- ical_phase %>%
  group_by(ID, exp_day, phase) %>%
  mutate(
    phase_cumsum = cumsum(count),
    relative_phase_cumsum = phase_cumsum - first(phase_cumsum)
  ) %>%
  ungroup()

ical_phase <- ical_phase %>%
  mutate(
    lights_phase = if_else(phase == "dark", "lights_off", "lights_on")
  )

phase_summary <- ical_phase %>%
  group_by(ID, SEX, BPA_EXPOSURE, lights_phase) %>%
  summarise(
    total_relative_counts = max(relative_phase_cumsum, na.rm = TRUE),
    .groups = "drop"
  )

##### STATS for females in lights off----
dark_fem <- phase_summary %>%
  filter(
    SEX == "F",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = dark_fem
)

##### STATS for females in lights on ----
light_fem <- phase_summary %>%
  filter(
    SEX == "F",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = light_fem
)

##### STATS for males in lights off----
dark_male <- phase_summary %>%
  filter(
    SEX == "M",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = dark_male
)

##### STATS for males in lights on----
light_male <- phase_summary %>%
  filter(
    SEX == "M",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = light_male
)

phase_plot_df <- phase_summary %>%
  group_by(SEX, lights_phase, BPA_EXPOSURE) %>%
  summarise(
    mean_counts = mean(total_relative_counts, na.rm = TRUE),
    sem_counts  = sd(total_relative_counts, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

phase_plot_df <- phase_plot_df %>%
  mutate(
    lights_phase = factor(
      lights_phase,
      levels = c("lights_off", "lights_on"),
      labels = c("Lights OFF", "Lights ON")
    )
  )

phase_summary <- phase_summary %>%
  mutate(
    lights_phase = factor(
      lights_phase,
      levels = c("lights_off", "lights_on"),
      labels = c("Lights OFF", "Lights ON")
    )
  )

## supplementary figure 8C (SF8C) locomotion separated by light and not for diet ----
sf8c <- ggplot(
  phase_plot_df,
  aes(
    x = lights_phase,
    y = mean_counts,
    fill = BPA_EXPOSURE
  )
) +
  # Bars (means)
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.7
  ) +
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = mean_counts - sem_counts,
      ymax = mean_counts + sem_counts
    ),
    width = 0.2,
    position = position_dodge(width = 0.75)
  ) +
  # Individual animal points
  geom_jitter(
    data = phase_summary,
    aes(
      x = lights_phase,
      y = total_relative_counts,
      color = BPA_EXPOSURE
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 2,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ SEX) +
  labs(
    x = "Light cycle",
    y = "Counts",
    fill = "BPA exposure",
    color = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  # Set BPA_EXPOSURE colors manually
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray"))

sf8c

# supplementary figure 8 (SF8A + SF8B + SF8V)----

sf8a <- sf8a + labs(tag = "A")
sf8b <- sf8b + labs(tag = "B")
sf8c <- sf8c + labs(tag = "C")

combined_plot_locomotion <- sf8a / sf8b / sf8c
combined_plot_locomotion


#CORRELATION (ALTERNATIVE ANALYSIS ) ----

echoMRI_data_cor <- echoMRI_data_comparisons %>%
  filter(n_measurement %in% c("0 wks", "22 wks")) %>%
  select(ID, SEX, BPA_EXPOSURE, n_measurement, Lean,COHORT) %>%
  pivot_wider(
    names_from  = n_measurement,
    values_from = Lean
  ) %>%
  mutate(
    delta_lean = `22 wks` - `0 wks`
  ) %>% 
  drop_na()


delta_loco <- relative_count_at_19  %>%
  select(ID, relative_total_count_19, SEX, BPA_EXPOSURE, cohort) 

delta_loco %>% count(SEX, BPA_EXPOSURE)

cor_df <- delta_loco %>%
  left_join(
    echoMRI_data_cor %>%
      select(ID, delta_lean),
    by = "ID"
  ) %>%
  drop_na(delta_lean)

cor_df %>%
  summarise(
    n_ID = n(),
    missing_lean = sum(is.na(delta_lean))
  )

#### females exposed to BPA----

cor_fem_bpa <- cor_df %>%
  filter(SEX == "F") %>% 
  filter(BPA_EXPOSURE =="YES")

cor.test(
  cor_fem_bpa$relative_total_count_19,
  cor_fem_bpa$delta_lean,
  method = "pearson"
)

#### females NON exposed to BPA----
cor_fem_no_bpa <- cor_df %>%
  filter(SEX == "F") %>% 
  filter(BPA_EXPOSURE =="NO")

cor.test(
  cor_fem_no_bpa$relative_total_count_19,
  cor_fem_no_bpa$delta_lean,
  method = "pearson"
)

alternative_a <-ggplot(cor_fem_bpa, aes(x = delta_lean, y = relative_total_count_19, color = BPA_EXPOSURE)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  labs(
    x = " Δ lean mass",
    y = "Counts over 24 hours",
    color = "BPA exposure"
  ) +
  theme_minimal()
alternative_a

alternative_b <- ggplot(cor_fem_no_bpa, aes(x = delta_lean, y = relative_total_count_19, color = BPA_EXPOSURE)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  labs(
    x = " Δ lean mass",
    y = "Counts over 24 hours",
    color = "BPA exposure"
  ) +
  theme_minimal()
alternative_b

# alternative figure 1 (alternative_a + alternative_b )----

alternative_a <- alternative_a + labs(tag = "A")
alternative_b <- alternative_b + labs(tag = "B")
combined_plot_cor <- alternative_a | alternative_b 
combined_plot_cor


##Total energy expenditure (TEE) ####

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  mutate(ID= as.numeric(ID))

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

ical_long_allheat <- bind_rows(ical_long15heat, ical_long16heat) %>% #this is key, here we combined
#filter(!ID %in% c(9404, 9403)) #these animals are responsible for a skew behavior of the normal curve in locomotion data
  filter(!ID ==9406)    #9406 has a  weird pattern in locomotion
  

ical_long_allheat %>% 
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(BPA_EXPOSURE,SEX) %>%
  summarise(n_ID = n_distinct(ID)) 

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
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(ID, exp_day) %>%
  mutate(
    relative_total_kcal_hr = kcal_hr_total - first(kcal_hr_total))

ical_long_allheat %>%
  filter(ID == min(ID)) %>%
  select(datetime, exp_day, hour, hour_ordered, kcal_hr_total) %>%
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

#### individual kcal_hr data----

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
    y = "relative kcal per hr"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 8)
  )

ical_long_allgroupedheat <- ical_long_allheat %>% 
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(hour_label, SEX, BPA_EXPOSURE,DIET_FORMULA) %>% 
  summarise(
    mean_kcal_hr = mean(relative_total_kcal_hr, na.rm = TRUE),
    sem_kcal_hr  = sd(relative_total_kcal_hr, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

dark_phase <- data.frame(
  xmin = 1,
  xmax = 11,
  ymin = -Inf,
  ymax = Inf
)

#### cumulative kcal_hr separated by diet----

ggplot(
  ical_long_allgroupedheat,
  aes(
    x = as.numeric(hour_label),
    y = mean_kcal_hr,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  # DARK PHASE BACKGROUND
  geom_rect(
    data = dark_phase,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.5
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
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_x_continuous(
    breaks = 1:24,
    labels = levels(ical_long_allgrouped$hour_label)
  ) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "Relative kcal per hr (mean ± SEM)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 10)
  )

# Summarize relative_total_kcal_hr at hour == 19 including SEX and BPA_EXPOSURE so this means kcal over 24h
relative_kcal_hr_at_19 <- ical_long_allheat %>%
  filter(hour == 19) %>%                # keep only rows at hour 19
  group_by(ID) %>%
  summarise(
    relative_total_kcal_hr_19 = first(relative_total_kcal_hr),  # value at 19:00
    cohort = first(cohort),
    SEX = first(SEX),
    BPA_EXPOSURE = first(BPA_EXPOSURE),
    .groups = "drop"
  ) %>% 
  filter(!ID ==9406)    #9406 has a  weird pattern in locomotion
  
relative_kcal_hr_at_19%>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 

####histogram check data distribution-----

ggplot(relative_kcal_hr_at_19, aes(x =  relative_total_kcal_hr_19 )) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  theme_classic() #so clearly we have left skewed curve

ggplot(relative_kcal_hr_at_19, aes(sample = relative_total_kcal_hr_19)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(SEX ~ BPA_EXPOSURE) +
  theme_classic()

shapiro.test(relative_kcal_hr_at_19$relative_total_kcal_hr_19) #our variable relative_total_kcal_hr_19 is not normally distributed.

#median
# Compute median and 25/75th percentiles for each group
summary_relative_kcal_hr_at_19 <- relative_kcal_hr_at_19 %>%
  group_by(BPA_EXPOSURE, SEX) %>%
  summarise(
    med = median(relative_total_kcal_hr_19),
    ymin = quantile(relative_total_kcal_hr_19, 0.25),
    ymax = quantile(relative_total_kcal_hr_19, 0.75)
  )

# Plot median
plot_summary_relative_kcal_hr_at_19 <- ggplot(relative_kcal_hr_at_19, aes(x = BPA_EXPOSURE, y = relative_total_kcal_hr_19)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_point(data = summary_relative_kcal_hr_at_19, aes(x = BPA_EXPOSURE, y = med), size = 3) +
  geom_errorbar(data = summary_relative_kcal_hr_at_19, aes(x = BPA_EXPOSURE, y = med, ymin = ymin, ymax = ymax), width = 0.2) +
  facet_wrap(~ SEX) +
  labs(
    x = "BPA exposure",
    y = "median TEE (kcal over 24g)"
  ) +
  theme_classic()
plot_summary_relative_kcal_hr_at_19 

#median collapsed by sex
# Compute median and 25/75th percentiles for each group
summary_relative_kcal_hr_at_19 <- relative_kcal_hr_at_19 %>%
  ungroup() %>% 
  group_by(SEX) %>%
  summarise(
    med = median(relative_total_kcal_hr_19),
    ymin = quantile(relative_total_kcal_hr_19, 0.25),
    ymax = quantile(relative_total_kcal_hr_19, 0.75)
  )

# Plot median
plot_summary_relative_kcal_hr_at_19 <- ggplot(relative_kcal_hr_at_19, aes(x = SEX, y = relative_total_kcal_hr_19)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_point(data = summary_relative_kcal_hr_at_19, aes(x = SEX, y = med), size = 3) +
  geom_errorbar(data = summary_relative_kcal_hr_at_19, aes(x = SEX, y = med, ymin = ymin, ymax = ymax), width = 0.2) +
 # facet_wrap(~ SEX) +
  labs(
    x = "BPA exposure",
    y = "median TEE (kcal over 24g)"
  ) +
  theme_classic()
plot_summary_relative_kcal_hr_at_19 #we expected males spent more calories than females because they tend to have more lean mass


# STATS kcal hr Run unpaired t-test within each sex ----

# Check variance within each sex
variance_checkheat <- relative_kcal_hr_at_19 %>%
  group_by(SEX) %>%
  summarise(
    var_test = list(
      var.test(
        relative_total_kcal_hr_19 ~ BPA_EXPOSURE, 
        data = cur_data()
      )
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    f_statistic = var_test$statistic,
    df1 = var_test$parameter[1],
    df2 = var_test$parameter[2],
    p_value = var_test$p.value
  ) %>%
  select(SEX, f_statistic, df1, df2, p_value)

variance_checkheat #so because the assumption of equal variance is violated, so Welch’s t-test is justified.

t_test_resultsheat <- relative_kcal_hr_at_19 %>%
  group_by(SEX) %>%
  summarise(
    t_test = list(t.test(
      relative_total_kcal_hr_19 ~ BPA_EXPOSURE,  # formula: group comparison
      data = cur_data(),
      var.equal = FALSE                        # Welch's t-test (does not assume equal variance)
    )),
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),
    df = map_dbl(t_test, ~ .x$parameter),
    p_value = map_dbl(t_test, ~ .x$p.value),
    mean_NO = map_dbl(t_test, ~ .x$estimate[1]),   # mean for BPA_EXPOSURE == NO
    mean_YES = map_dbl(t_test, ~ .x$estimate[2])   # mean for BPA_EXPOSURE == YES
  ) %>%
  select(SEX, mean_NO, mean_YES, t_statistic, df, p_value)

t_test_resultsheat

# Compute mean and SEM per SEX × BPA_EXPOSURE
summary_19heat <- relative_kcal_hr_at_19 %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_kcal_hr = mean(relative_total_kcal_hr_19),
    sem_kcal_hr = sd(relative_total_kcal_hr_19) / sqrt(n()),
    .groups = "drop"
  )

# T-test results
t_test_summaryheat <- relative_kcal_hr_at_19 %>%
  group_by(SEX) %>%
  do(tidy(t.test(relative_total_kcal_hr_19 ~ BPA_EXPOSURE, data = .))) %>%
  ungroup() %>%
  select(SEX, estimate1, estimate2, statistic, p.value)

t_test_summaryheat

## supplementary figure 9A (SF9A) kcal over 24h ----
sf9a <- ggplot(summary_19heat, aes(x = BPA_EXPOSURE, y = mean_kcal_hr, fill = BPA_EXPOSURE)) +
  geom_col(alpha = 0.7, width = 0.6) +  # bars
  geom_errorbar(aes(ymin = mean_kcal_hr - sem_kcal_hr, ymax = mean_kcal_hr + sem_kcal_hr),
                width = 0.2) +           # SEM
  geom_jitter(
    data = relative_kcal_hr_at_19,
    aes(x = BPA_EXPOSURE, y = relative_total_kcal_hr_19, color = BPA_EXPOSURE),
    width = 0.15,
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ SEX) +
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  labs(
    x = "BPA Exposure",
    y = "kcal over 24h",
    fill = "BPA Exposure",
    color = "BPA Exposure"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )+
  coord_cartesian(ylim = c(0, 20))  # <-- set y-axis limits

sf9a


# Summarize relative_total_kcal_hr at hour == 19 including SEX and BPA_EXPOSURE and DIET so this means kcal spend over 24h
relative_kcal_hr_at_19_diet <- ical_long_allheat %>%
  filter(hour == 19) %>%
  group_by(ID) %>%
  summarise(
    relative_total_kcal_hr_19 = first(relative_total_kcal_hr),
    cohort = first(cohort),
    SEX = first(SEX),
    BPA_EXPOSURE = first(BPA_EXPOSURE),
    DIET_FORMULA = first(DIET_FORMULA),
    .groups = "drop"
  ) %>%
  mutate(
    DIET_FORMULA = dplyr::recode(DIET_FORMULA,
                                 "D12450Hi" = "HCD",
                                 "D12451i"  = "HFD")
  )


# females
model_fheat <- lm(relative_total_kcal_hr_19 ~ BPA_EXPOSURE * DIET_FORMULA,
              data = relative_kcal_hr_at_19_diet %>% filter(SEX == "F"))

anova(model_fheat)
summary(model_fheat)

# males
model_mheat <- lm(relative_total_kcal_hr_19 ~ BPA_EXPOSURE * DIET_FORMULA,
              data = relative_kcal_hr_at_19_diet %>% filter(SEX == "M"))

anova(model_mheat)
summary(model_mheat)


# females
emmeans(model_fheat, pairwise ~ BPA_EXPOSURE | DIET_FORMULA)

# males
emmeans(model_mheat, pairwise ~ BPA_EXPOSURE | DIET_FORMULA)

## supplementary figure 9B (SF9B) kcal over 24h separated by diet ----
# Mean and SEM per SEX × BPA × DIET
summary_19_dietheat <- relative_kcal_hr_at_19_diet %>%
  group_by(SEX, DIET_FORMULA, BPA_EXPOSURE) %>%
  summarise(
    mean_kcal_hr = mean(relative_total_kcal_hr_19),
    sem_kcal_hr = sd(relative_total_kcal_hr_19) / sqrt(n()),
    .groups = "drop"
  )

sf9b <- ggplot(summary_19_dietheat,
               aes(x = DIET_FORMULA,
                   y = mean_kcal_hr,
                   fill = BPA_EXPOSURE)) +
  
  geom_col(position = position_dodge(width = 0.7),
           alpha = 0.7,
           width = 0.6) +
  
  geom_errorbar(
    aes(ymin = mean_kcal_hr - sem_kcal_hr,
        ymax = mean_kcal_hr + sem_kcal_hr),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  
  geom_jitter(
    data = relative_kcal_hr_at_19_diet,
    aes(x = DIET_FORMULA,
        y = relative_total_kcal_hr_19,
        color = BPA_EXPOSURE),
    position = position_jitterdodge(jitter.width = 0.15,
                                    dodge.width = 0.7),
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ SEX) +
  
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  
  labs(
    x = "Diet",
    y = "kcal over 24h",
    fill = "BPA Exposure",
    color = "BPA Exposure"
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

sf9b

## kcal analysis separated by light period ----

ical_phaseheat <- ical_long_allheat %>%
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  mutate(
    phase = case_when(
      hour >= 20 | hour < 6  ~ "dark",   # 20:00–05:59
      hour >= 6  & hour < 20 ~ "light"   # 06:00–19:59
    )
  )

ical_phaseheat <- ical_phaseheat %>%
  mutate(
    phase_hour_order = case_when(
      phase == "dark"  & hour >= 20 ~ hour,
      phase == "dark"  & hour < 6   ~ hour + 24,
      phase == "light"              ~ hour
    )
  ) %>%
  arrange(ID, exp_day, phase, phase_hour_order)

ical_phaseheat <- ical_phaseheat %>%
  group_by(ID, exp_day, phase) %>%
  mutate(
    phase_cumsum = cumsum(kcal_hr),
    relative_phase_cumsum = phase_cumsum - first(phase_cumsum)
  ) %>%
  ungroup()

ical_phaseheat <- ical_phaseheat %>%
  mutate(
    lights_phase = if_else(phase == "dark", "lights_off", "lights_on")
  )

phase_summaryheat <- ical_phaseheat %>%
  group_by(ID, SEX, BPA_EXPOSURE, lights_phase) %>%
  summarise(
    total_relative_kcal_hr = max(relative_phase_cumsum, na.rm = TRUE),
    .groups = "drop"
  )

##### STATS for females in lights off----
dark_femheat <- phase_summaryheat %>%
  filter(
    SEX == "F",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_kcal_hr ~ BPA_EXPOSURE,
  data = dark_femheat
)

##### STATS for females in lights on ----
light_femheat <- phase_summaryheat %>%
  filter(
    SEX == "F",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_kcal_hr ~ BPA_EXPOSURE,
  data = light_femheat
)

##### STATS for males in lights off----
dark_maleheat <- phase_summaryheat %>%
  filter(
    SEX == "M",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_kcal_hr ~ BPA_EXPOSURE,
  data = dark_maleheat
)

##### STATS for males in lights on----
light_maleheat <- phase_summaryheat %>%
  filter(
    SEX == "M",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_kcal_hr ~ BPA_EXPOSURE,
  data = light_maleheat
)

phase_plot_dfheat <- phase_summaryheat %>%
  group_by(SEX, lights_phase, BPA_EXPOSURE) %>%
  summarise(
    mean_kcal_hr = mean(total_relative_kcal_hr, na.rm = TRUE),
    sem_kcal_hr  = sd(total_relative_kcal_hr, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

phase_plot_dfheat <- phase_plot_dfheat %>%
  mutate(
    lights_phase = factor(
      lights_phase,
      levels = c("lights_off", "lights_on"),
      labels = c("Lights OFF", "Lights ON")
    )
  )

phase_summaryheat <- phase_summaryheat %>%
  mutate(
    lights_phase = factor(
      lights_phase,
      levels = c("lights_off", "lights_on"),
      labels = c("Lights OFF", "Lights ON")
    )
  )

## supplementary figure 9C (SF9C) kcal separated by light and not for diet ----
sf9c <- ggplot(
  phase_plot_dfheat,
  aes(
    x = lights_phase,
    y = mean_kcal_hr,
    fill = BPA_EXPOSURE
  )
) +
  # Bars (means)
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.7
  ) +
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = mean_kcal_hr - sem_kcal_hr,
      ymax = mean_kcal_hr + sem_kcal_hr
    ),
    width = 0.2,
    position = position_dodge(width = 0.75)
  ) +
  # Individual animal points
  geom_jitter(
    data = phase_summaryheat,
    aes(
      x = lights_phase,
      y = total_relative_kcal_hr,
      color = BPA_EXPOSURE
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 2,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ SEX) +
  labs(
    x = "Light cycle",
    y = "Kcal",
    fill = "BPA exposure",
    color = "BPA exposure"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  # Set BPA_EXPOSURE colors manually
  scale_fill_manual(values = c("NO" = "black", "YES" = "gray")) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray"))

sf9c

# supplementary figure 9 (SF9A + SF9B + SF9C)----

sf9a <- sf9a + labs(tag = "A")
sf9b <- sf9b + labs(tag = "B")
sf9c <- sf9c + labs(tag = "C")

combined_plot_kcal_hr <- sf9a / sf9b / sf9c
combined_plot_kcal_hr


##Respiratory Exchange Ratio analysis####
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
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(BPA_EXPOSURE,SEX) %>%
  summarise(n_ID = n_distinct(ID)) 

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
  group_by(ID, exp_day) %>%
  mutate(RER_mean = mean(RER)) %>%
  ungroup() %>% 
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(ID, exp_day) 

ical_long_allRER %>%
  filter(ID == min(ID)) %>%
  select(datetime, exp_day, hour, hour_ordered, RER_mean) %>%
  arrange(datetime)

ical_long_allRER <- ical_long_allRER%>%
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

#### individual RER data----

ggplot( 
  ical_long_allRER,
  aes(
    x = hour_label,
    y = RER_mean,
    group = interaction(ID, exp_day)
  )
) +
  geom_line(size = 1) +
  facet_wrap(~ ID) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "mean RER"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 8)
  )

ical_long_allgroupedRER <- ical_long_allRER %>% 
  filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
  group_by(hour_label, SEX, BPA_EXPOSURE,DIET_FORMULA) %>% 
  summarise(
    mean_RER = mean(RER, na.rm = TRUE),
    sem_RER  = sd(RER, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n_ID = n_distinct(ID),
    .groups = "drop"
  )

dark_phase <- data.frame(
  xmin = 1,
  xmax = 11,
  ymin = -Inf,
  ymax = Inf
)

#### cumulative RER separated by diet----

ggplot(
  ical_long_allgroupedRER,
  aes(
    x = as.numeric(hour_label),
    y = mean_RER,
    group = BPA_EXPOSURE,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  # DARK PHASE BACKGROUND
  geom_rect(
    data = dark_phase,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.5
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
  facet_wrap(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  scale_x_continuous(
    breaks = 1:24,
    labels = levels(ical_long_allgroupedRER$hour_label)
  ) +
  labs(
    x = "Time (20:00 → 19:00)",
    y = "RER (mean ± SEM)",
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

#Novel object recognition test ----


####Data import----

data <- read_csv("~/Documents/GitHub/data/data/Data.csv")

# exploration of the data about what is happening with the behavior of mice in each stage

# hypothesis N1: animals prefers periphery than center ----
# we can expect the animals will not have a very clear pattern in terms of spend more time in the center or in the periphery
#however, mice trend to avoid open spaces so naturally we should expect that they will spend more time in the periphery area

data_raw <- read_csv("~/Documents/GitHub/data/data/Data.csv") %>% 
  select(Animal, Stage, `Segment of test`,
         `Center : time (s)`, `Periphery : time (s)`)
data_summary <- data_raw %>% 
  group_by(Stage, `Segment of test`) %>% 
  summarise(
    n = n(),
    mean_time_center = mean(`Center : time (s)`, na.rm = TRUE),
    sem_time_center  = sd(`Center : time (s)`, na.rm = TRUE)/sqrt(n),
    mean_time_periphery = mean(`Periphery : time (s)`, na.rm = TRUE),
    sem_time_periphery  = sd(`Periphery : time (s)`, na.rm = TRUE)/sqrt(n),
    .groups = "drop"
  )
data_raw_long <- data_raw %>% 
  pivot_longer(
    cols = c(`Center : time (s)`, `Periphery : time (s)`),
    names_to = "Location",
    values_to = "time"
  )

data_summary_long <- data_summary %>%
  pivot_longer(
    cols = -c(Stage, `Segment of test`, n),
    names_to = c(".value", "Location"),
    names_pattern = "(mean_time|sem_time)_(.*)")
    
    data_raw_long$Stage <- factor(data_raw_long$Stage,
                                  levels = c("Habituation",
                                             "Familiarization",
                                             "Recognition"))
    
    data_summary_long$Stage <- factor(data_summary_long$Stage,
                                      levels = c("Habituation",
                                                 "Familiarization",
                                                 "Recognition"))
ggplot() +
  
  geom_col(data = data_summary_long,
           aes(x = `Segment of test`,
               y = mean_time,
               fill = Location),
           position = position_dodge(width = 0.8),
           alpha = 0.6) +
  
  geom_errorbar(data = data_summary_long,
                aes(x = `Segment of test`,
                    ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time,
                    group = Location),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  
  geom_jitter(data = data_raw_long,
              aes(x = `Segment of test`,
                  y = time,
                  color = Location),
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.8
              ),
              alpha = 0.8,
              size = 2,
              show.legend = FALSE) +
  facet_wrap(~Stage) +
  labs(y = "Time (seconds)",
       x = "Segment of the test",
       fill = "Zone") +
  theme_classic()

#conclusion> so great animals spent more time in the periphery than in the center of the arena
#conclusion 2 > in the recognition stage animals spend more time interacting with the objects that are in the center

# hypothesis N2: In the familiarization stage animals will spend the same time with both objects ----
# we can expect the animals will spend the same amount of time with the two objects IF THOSE ARE THE SAME OBJECTS

data_fam <- read_csv("~/Documents/GitHub/data/data/Data.csv") %>% 
  select(Animal, Stage, `Segment of test`,
         `Familiar Object 1 : time investigating (s)`,
         `Familiar Object 2 : time investigating (s)`) %>% 
  filter(Stage =='Familiarization') %>% 
  drop_na()

data_summary_fam <- data_fam %>% 
  group_by(`Segment of test`) %>% 
  summarise(
    n = n(),
    mean_time_fam1 = mean( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE),
    sem_time_fam1  = sd( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
    mean_time_fam2 = mean( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE),
    sem_time_fam2  = sd( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
    .groups = "drop"
  )
data_raw_long_fam <- data_fam %>% 
  pivot_longer(
    cols = c(`Familiar Object 1 : time investigating (s)`,  `Familiar Object 2 : time investigating (s)`),
    names_to = "Object",
    values_to = "time"
  )

data_raw_long_fam <- data_raw_long_fam %>%
  mutate(Object = as.character(Object),
         Object = case_when(
           Object == "Familiar Object 1 : time investigating (s)" ~ "1",
           Object == "Familiar Object 2 : time investigating (s)" ~ "2",
           TRUE ~ Object  # leave any unexpected values as-is
         ),
         Object = factor(Object, levels = c("1", "2")))

data_summary_long_fam <- data_summary_fam %>%
  pivot_longer(
    cols = -`Segment of test`,
    names_to = c(".value", "Object"),
    names_pattern = "(mean_time|sem_time)_fam(\\d)"
  ) %>%
  mutate(Object = case_when(
    Object == "1" ~ "1",
    Object == "2" ~ "2"
  ),
  Object = factor(Object, levels = c("1", "2"))) %>% 
  drop_na()


ggplot(data_summary_long_fam, aes(x = Object, y = mean_time, fill = Object)) +
  geom_col(alpha = 0.6) +
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.2) +
  geom_jitter(data = data_raw_long_fam,
              aes(x = Object, y = time, color = Object),
              width = 0.15, size = 2, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~`Segment of test`) +
  labs(x = "Object", y = "Time Investigating (s)", fill = "Object") +
  theme_classic()

##### STATS----
# Make sure we are using only Familiarization stage
data_fam_wide <- data_raw_long_fam %>%
  select(Animal, `Segment of test`, Object, time) %>%
  pivot_wider(
    names_from = Object,
    values_from = time,
    names_prefix = "Object_"
  )
data_fam_wide

t_test_results <- data_fam_wide %>%
  group_by(`Segment of test`) %>%
  summarise(
    t_test = list(t.test(Object_1, Object_2, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(t_test = map(t_test, broom::tidy)) %>%
  unnest(t_test)
t_test_results

#conclusion 1> For both segments, p > 0.05, so mice did not spend significantly different time investigating Object 1 vs Object 2.
#conclusion 2> The small positive estimate just means Object 1 time was slightly higher on average, but not significant.
#conclusion 3> animals explore the same amount of time the same object which is what we expected so they are ready for the recognition stage


#projections> to calculate discrimination indexes (DI) defined as (Tn-Tf)/(Tn+Tf) and/or to compare the novel object preference (%) = ((Tn)/(Tn+Tf))*100
# Tn: time spent with the novel object
# Tf: time spent with the familiar object

#hypothesis N3: In the familiarization stage animals have a discrimination index (Di) close to zero ----

dat_fam_di <- data_fam %>% 
  ungroup() %>% 
  group_by(Animal,`Segment of test`) %>% 
  mutate(Di =(`Familiar Object 1 : time investigating (s)` - `Familiar Object 2 : time investigating (s)`)/(`Familiar Object 1 : time investigating (s)` + `Familiar Object 2 : time investigating (s)`))

data_summary_fam_di <- dat_fam_di %>% 
  ungroup() %>% 
  group_by(`Segment of test`) %>% 
  summarise(
    n = n(),
    mean_time_Di = mean(Di),
    sem_time_Di  = sd(Di)/sqrt(n),
    .groups = "drop"
  ) %>% 
  ungroup()

data_raw_long_fam_di <- dat_fam_di %>% 
  pivot_longer(
    cols = c(Di),
    names_to = "Object",
    values_to = "Di_value"
  )

ggplot(data_summary_fam_di, aes(x = `Segment of test`, y = mean_time_Di)) +
  geom_col(alpha = 0.6, fill = "grey70") +
  geom_errorbar(aes(ymin = mean_time_Di - sem_time_Di,
                    ymax = mean_time_Di + sem_time_Di),
                width = 0.2) +
  geom_jitter(data = dat_fam_di,
              aes(x = `Segment of test`, y = Di),
              width = 0.1, size = 2, alpha = 0.8) +
  labs(x = "Segment of test",
       y = "Discrimination index (Di)") +
  theme_classic()

#####STATS----

data_fam_di_wide <- dat_fam_di %>%
  ungroup() %>%
  select(Animal, `Segment of test`, Di) %>%
  distinct() %>%   # important in case mutate created duplicates
  pivot_wider(
    names_from = `Segment of test`,
    values_from = Di
  )
t_test_di <- t.test(
  data_fam_di_wide$`0 - 300 secs.`,
  data_fam_di_wide$`300 - 600 secs.`,
  paired = TRUE
)

t_test_di
#The discrimination index remains stable across the two time segments during familiarization.

# Since this is the familiarization phase, your real hypothesis might actually be: Di ≈ 0 (no preference)

dat_fam_di %>%
  group_by(`Segment of test`) %>%
  summarise(
    t_test = list(t.test(Di, mu = 0)),
    .groups = "drop"
  ) %>%
  mutate(t_test = purrr::map(t_test, broom::tidy)) %>%
  tidyr::unnest(t_test)

#During familiarization, mice did not show a preference for either object in either time segment.

####DATA IMPORT COHORT 16 4 MONTH OLD----


METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  rename(Animal = ID)

data2 <- read_csv("~/Documents/GitHub/data/data/Data-2.csv") %>% 
  mutate(Animal = gsub("\\.", "-", as.character(Animal))) %>% 
left_join(METABPA, by= "Animal") 

# hypothesis N1: animals prefers periphery than center ----
# we can expect the animals will trend to avoid open spaces so we expect that they will spend more time in the periphery area


data2_summary <- data2 %>% 
  group_by(Stage, `Segment of test`) %>% 
  summarise(
    n = n(),
    mean_time_center = mean(`Center : time (s)`, na.rm = TRUE),
    sem_time_center  = sd(`Center : time (s)`, na.rm = TRUE)/sqrt(n),
    mean_time_periphery = mean(`Periphery : time (s)`, na.rm = TRUE),
    sem_time_periphery  = sd(`Periphery : time (s)`, na.rm = TRUE)/sqrt(n),
    .groups = "drop"
  )

data2_raw_long <- data2 %>% 
  pivot_longer(
    cols = c(`Center : time (s)`, `Periphery : time (s)`),
    names_to = "Location",
    values_to = "time"
  )

data2_summary_long <- data2_summary %>%
  pivot_longer(
    cols = -c(Stage, `Segment of test`, n),
    names_to = c(".value", "Location"),
    names_pattern = "(mean_time|sem_time)_(.*)")

data2_raw_long$Stage <- factor(data2_raw_long$Stage,
                              levels = c("Habituation",
                                         "Familiarization",
                                         "Recognition"))

data2_summary_long$Stage <- factor(data2_summary_long$Stage,
                                  levels = c("Habituation",
                                             "Familiarization",
                                             "Recognition"))
ggplot() +
  geom_col(data = data2_summary_long,
           aes(x = `Segment of test`,
               y = mean_time,
               fill = Location),
           position = position_dodge(width = 0.8),
           alpha = 0.6) +
  
  geom_errorbar(data = data2_summary_long,
                aes(x = `Segment of test`,
                    ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time,
                    group = Location),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  
  geom_jitter(data = data2_raw_long,
              aes(x = `Segment of test`,
                  y = time,
                  color = Location),
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.8
              ),
              alpha = 0.8,
              size = 2,
              show.legend = FALSE) +
  facet_grid(~Stage*BPA_EXPOSURE) +
  labs(y = "Time (seconds)",
       x = "Segment of the test",
       fill = "Zone") +
  theme_classic()

#conclusion> so great animals spent more time in the periphery than in the center of the arena

# hypothesis N2: In the familiarization stage animals will spend the same time with both objects ----
# we can expect the animals will spend the same amount of time with the two objects

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
  rename(Animal = ID)

data2 <- read_csv("~/Documents/GitHub/data/data/Data-2.csv") %>% 
  mutate(Animal = gsub("\\.", "-", as.character(Animal))) %>% 
  left_join(METABPA, by= "Animal") 

data2_fam <- data2 %>% 
  select(Animal, Stage, `Segment of test`,
         `Familiar Object 1 : time investigating (s)`,
         `Familiar Object 2 : time investigating (s)`,
         SEX, BPA_EXPOSURE) %>% 
  filter(Stage =='Familiarization') %>% 
  drop_na()

data2_summary_fam <- data2_fam %>% 
  group_by(`Segment of test`, BPA_EXPOSURE) %>% 
  summarise(
    n = n(),
    mean_time_fam1 = mean( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE),
    sem_time_fam1  = sd( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
    mean_time_fam2 = mean( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE),
    sem_time_fam2  = sd( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
    .groups = "drop"
  )
data2_raw_long_fam <- data2_fam %>% 
  pivot_longer(
    cols = c(`Familiar Object 1 : time investigating (s)`,  `Familiar Object 2 : time investigating (s)`),
    names_to = "Object",
    values_to = "time"
  )

data2_raw_long_fam <- data2_raw_long_fam %>%
  mutate(Object = as.character(Object),
         Object = case_when(
           Object == "Familiar Object 1 : time investigating (s)" ~ "1",
           Object == "Familiar Object 2 : time investigating (s)" ~ "2",
           TRUE ~ Object  # leave any unexpected values as-is
         ),
         Object = factor(Object, levels = c("1", "2")))

data2_summary_long_fam <- data2_summary_fam %>%
  pivot_longer(
    cols = -c(`Segment of test`, BPA_EXPOSURE),  
    names_to = c(".value", "Object"),
    names_pattern = "(mean_time|sem_time)_fam(\\d)"
  ) %>%
  mutate(Object = factor(Object, levels = c("1", "2"))) %>%
  drop_na()

ggplot(data2_summary_long_fam, aes(x = Object, y = mean_time, fill = Object)) +
  geom_col(alpha = 0.6) +
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.2) +
  geom_jitter(data = data2_raw_long_fam,
              aes(x = Object, y = time, color = Object),
              width = 0.15, size = 2, alpha = 0.8, show.legend = FALSE) +
  facet_grid(~`Segment of test`*BPA_EXPOSURE) +
  labs(x = "Object", y = "Time Investigating (s)", fill = "Object") +
  theme_classic()

##### STATS----

data2_fam_wide <- data2_raw_long_fam %>%
  select(Animal, BPA_EXPOSURE, `Segment of test`, Object, time) %>%
  pivot_wider(
    names_from = Object,
    values_from = time,
    names_prefix = "Object_"
  )
t2_test_results <- data2_fam_wide %>%
  group_by(`Segment of test`, BPA_EXPOSURE) %>%
  summarise(
    t_test = list(t.test(Object_1, Object_2, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(t_test = map(t_test, tidy)) %>%
  unnest(t_test)

t2_test_results


#conclusion 1> This is not a reliable statistical test with n = 2. The test has almost no power.
#conclusion 2> Multiple paired t-tests are not ideal, p-values are meaningless with n = 2

#hypothesis N3: In the familiarization stage animals have a discrimination index (Di) close to zero ----

dat2_fam_di <- data2_fam %>% 
  ungroup() %>% 
  group_by(Animal,`Segment of test`,BPA_EXPOSURE) %>% 
  mutate(Di =(`Familiar Object 1 : time investigating (s)` - `Familiar Object 2 : time investigating (s)`)/(`Familiar Object 1 : time investigating (s)` + `Familiar Object 2 : time investigating (s)`))

data2_summary_fam_di <- dat2_fam_di %>% 
  ungroup() %>% 
  group_by(`Segment of test`,BPA_EXPOSURE) %>% 
  summarise(
    n = n(),
    mean_time_Di = mean(Di),
    sem_time_Di  = sd(Di)/sqrt(n),
    .groups = "drop"
  ) %>% 
  ungroup()

data2_raw_long_fam_di <- dat2_fam_di %>% 
  pivot_longer(
    cols = c(Di),
    names_to = "Object",
    values_to = "Di_value"
  )

ggplot(data2_summary_fam_di, aes(x = `Segment of test`, y = mean_time_Di)) +
  geom_col(alpha = 0.6, fill = "grey70") +
  geom_errorbar(aes(ymin = mean_time_Di - sem_time_Di,
                    ymax = mean_time_Di + sem_time_Di),
                width = 0.2) +
  geom_jitter(data = dat2_fam_di,
              aes(x = `Segment of test`, y = Di),
              width = 0.1, size = 2, alpha = 0.8) +
  labs(x = "Segment of test",
       y = "Discrimination index (Di)") +
  facet_wrap(~BPA_EXPOSURE)
  theme_classic()
  
  
##ALL COHORTS ANALYSIS----

#Read me: 
#Data.csv- raw data for cohort 15 (7/15/25)
#Data-2.csv- Raw data for cohort 16 (11/10/2025) #Zander comment (12/21/25): The Data-2 was cut short. Due to technical issues we had to start late and we had to rush because another lab group was telling us they needed the room. We for sure got the 0-5 minute window for all mice though.
#Data-3.csv- Raw data for cohort 15 and cohort 16 (10/13/2025) 
#Data-4.csv- raw data for cohort 16 (1/15/2026) #habituation was not correctly recordered. Zander did not recorded 300 - 600 secs section of the test in the habituation phase. Also he ran the same day habituation, familiarization and recognition phases

  
####Data import----
  
METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv") %>% 
filter(grepl("-", ID) | ID %in% c("9408", "9409"))   #I want to keep just IDs that were measured for NORT for cohort 1, 2 and 3

  
data<- read_csv("~/Documents/GitHub/data/data/Data.csv") %>% #familiarization data from 9380-3 and 9391-1 was not recorded property 
  rename(ID = Animal) %>% 
  mutate(ID= as.character(ID),ID = case_when(
    ID == "1"  ~ "9389-1",
    ID == "2"  ~ "9389-2",
    ID == "3"  ~ "9389-3",
    ID == "4"  ~ "9389-4",
    ID == "5"  ~ "9380-1",
    ID == "6"  ~ "9380-2",
    ID == "7"  ~ "9380-3",
    ID == "8"  ~ "9391-1",
    ID == "9"  ~ "9391-2",
    ID == "10" ~ "9392-1",
    ID == "11" ~ "9392-2",
    ID == "12" ~ "9392-3",
    ID == "13" ~ "9385-1",
    ID == "14" ~ "9385-2",
    ID == "15" ~ "9385-3",
    ID == "16" ~ "9385-4",
    TRUE       ~ ID))%>% 
  left_join(METABPA, by= "ID")
  
data2 <- read_csv("~/Documents/GitHub/data/data/Data-2.csv") %>% 
  rename(ID = Animal) %>% 
    mutate(ID = gsub("\\.", "-", as.character(ID))) %>% 
    left_join(METABPA, by= "ID")

data3 <- read_csv("~/Documents/GitHub/data/data/Data-3.csv") %>% 
  rename(ID = Animal) %>% 
  mutate(ID = gsub("\\.", "-", as.character(ID))) %>% 
  left_join(METABPA, by= "ID")

data4 <- read_csv("~/Documents/GitHub/data/data/Data-4.csv") %>% #habituation data was not recorded property 
  rename(ID = Animal) %>% 
  mutate(ID = gsub("\\.", "-", as.character(ID))) %>% 
  left_join(METABPA, by= "ID")

all_data <- bind_rows(data, data2, data3, data4) %>% 
  mutate(
    Date = mdy(Date),
    DOB  = mdy(DOB),
    age_weeks = as.numeric(difftime(Date, DOB, units = "weeks")))%>%
  mutate(age_week_round = round(age_weeks),
         age_week_round = case_when(
           age_week_round %in% c(13, 14) ~ 11,
           age_week_round ==17 ~ 14,# this is 17 wks old, 14 wks with chow
           age_week_round ==27 ~ 24, # this is 27 wks old, 24 wks with cho
           TRUE ~ age_week_round               # keep other ages as they are
         )) 

all_data   %>% 
  group_by(SEX,BPA_EXPOSURE,COHORT) %>%
  summarise(n_ID = n_distinct(ID))

# hypothesis N1: animals prefers periphery than center ----
# we can expect the animals will trend to avoid open spaces so they will spend more time in the periphery area
#this behavior will be independent of the age, BPA exposure, sex and Stage
  
data_raw <- all_data %>% 
 select(ID, Stage, `Segment of test`,
           `Center : time (s)`, `Periphery : time (s)`, BPA_EXPOSURE, age_week_round,SEX,COHORT)

data_raw   %>% 
  group_by(SEX,BPA_EXPOSURE,COHORT) %>%
  summarise(n_ID = n_distinct(ID))

  data_summary <- data_raw %>% 
    group_by(Stage, `Segment of test`,age_week_round) %>% #there is not enough females to split the data by sex
    summarise(
      n = n(),
      mean_time_center = mean(`Center : time (s)`, na.rm = TRUE),
      sem_time_center  = sd(`Center : time (s)`, na.rm = TRUE)/sqrt(n),
      mean_time_periphery = mean(`Periphery : time (s)`, na.rm = TRUE),
      sem_time_periphery  = sd(`Periphery : time (s)`, na.rm = TRUE)/sqrt(n),
      .groups = "drop"
    )

  data_raw_long <- data_raw %>% 
    pivot_longer(
      cols = c(`Center : time (s)`, `Periphery : time (s)`),
      names_to = "Location",
      values_to = "time"
    )
  
  data_summary_long <- data_summary %>%
    pivot_longer(
      cols = -c(Stage, `Segment of test`, n, age_week_round),
      names_to = c(".value", "Location"),
      names_pattern = "(mean_time|sem_time)_(.*)")
  
  data_raw_long$Stage <- factor(data_raw_long$Stage,
                                levels = c("Habituation",
                                           "Familiarization",
                                           "Recognition"))
  
  data_summary_long$Stage <- factor(data_summary_long$Stage,
                                    levels = c("Habituation",
                                               "Familiarization",
                                               "Recognition"))
  ggplot() +
    
    geom_col(data = data_summary_long,
             aes(x = `Segment of test`,
                 y = mean_time,
                 fill = Location),
             position = position_dodge(width = 0.8),
             alpha = 0.6) +
    
    geom_errorbar(data = data_summary_long,
                  aes(x = `Segment of test`,
                      ymin = mean_time - sem_time,
                      ymax = mean_time + sem_time,
                      group = Location),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    
    geom_jitter(data = data_raw_long,
                aes(x = `Segment of test`,
                    y = time,
                    color = Location),
                position = position_jitterdodge(
                  jitter.width = 0.15,
                  dodge.width = 0.8
                ),
                alpha = 0.8,
                size = 2,
                show.legend = FALSE) +
    facet_wrap(~Stage*age_week_round) +
    labs(y = "Time (seconds)",
         x = "Segment of the test",
         fill = "Zone") +
    theme_classic()
  
#conclusion> so great animals spent more time in the periphery than in the center of the arena in all Stages
#conclusion 2 > In the recognition stage for the 27-week-old animals, the 5–10 minute period was not fully recorded. Zander mentioned that this was because another group was scheduled to use the room
  
# hypothesis N2: In the familiarization stage animals will spend the same time with both objects ----
# we can expect the animals will spend the same amount of time with the two identical objects independent of age, Segment of the test and BPA exposure 
  
data_fam <- all_data %>% 
    select(ID, Stage, `Segment of test`,
           `Familiar Object 1 : time investigating (s)`,
           `Familiar Object 2 : time investigating (s)`,
           BPA_EXPOSURE, SEX, COHORT, age_week_round) %>% 
    filter(Stage =='Familiarization') %>% 
    drop_na()
  
  data_fam   %>% 
    group_by(SEX,BPA_EXPOSURE,COHORT) %>%
    summarise(n_ID = n_distinct(ID))  
  
  
  data_summary_fam <- data_fam %>% 
    group_by(`Segment of test`, age_week_round) %>% 
    summarise(
      n = n(),
      mean_time_fam1 = mean( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE),
      sem_time_fam1  = sd( `Familiar Object 1 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
      mean_time_fam2 = mean( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE),
      sem_time_fam2  = sd( `Familiar Object 2 : time investigating (s)`, na.rm = TRUE)/sqrt(n),
      .groups = "drop"
    )
  
  data_raw_long_fam <- data_fam %>% 
    pivot_longer(
      cols = c(`Familiar Object 1 : time investigating (s)`,  `Familiar Object 2 : time investigating (s)`),
      names_to = "Object",
      values_to = "time"
    )
  
  data_raw_long_fam <- data_raw_long_fam %>%
    mutate(Object = as.character(Object),
           Object = case_when(
             Object == "Familiar Object 1 : time investigating (s)" ~ "1",
             Object == "Familiar Object 2 : time investigating (s)" ~ "2",
             TRUE ~ Object  # leave any unexpected values as-is
           ),
           Object = factor(Object, levels = c("1", "2")))
  
  data_summary_long_fam <- data_summary_fam %>%
    pivot_longer(
      cols = -c(`Segment of test`, age_week_round),
      names_to = c(".value", "Object"),
      names_pattern = "(mean_time|sem_time)_fam(\\d)"
    ) %>%
    mutate(Object = case_when(
      Object == "1" ~ "1",
      Object == "2" ~ "2"
    ),
    Object = factor(Object, levels = c("1", "2"))) %>% 
    drop_na()
  
  
  ggplot(data_summary_long_fam, aes(x = Object, y = mean_time, fill = Object)) +
    geom_col(alpha = 0.6) +
    geom_errorbar(aes(ymin = mean_time - sem_time,
                      ymax = mean_time + sem_time),
                  width = 0.2) +
    geom_jitter(data = data_raw_long_fam,
                aes(x = Object, y = time, color = Object),
                width = 0.15, size = 2, alpha = 0.8, show.legend = FALSE) +
    facet_wrap(~`Segment of test`*age_week_round) +
    labs(x = "Object", y = "Time Investigating (s)", fill = "Object") +
    theme_classic()
  
##### STATS----
# Make sure we are using only Familiarization stage
  data_fam_wide <- data_raw_long_fam %>%
    select(ID, `Segment of test`, Object, time,age_week_round) %>%
    pivot_wider(
      names_from = Object,
      values_from = time,
      names_prefix = "Object_"
    )
  data_fam_wide
  
  t_test_results <- data_fam_wide %>%
    group_by(`Segment of test`,age_week_round) %>%
    summarise(
      t_test = list(t.test(Object_1, Object_2, paired = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(t_test = map(t_test, broom::tidy)) %>%
    unnest(t_test)
  t_test_results
  
#conclusion: Across age groups, BPA exposure conditions, and test segments, 
#animals did not show a significant preference for either familiar object 
#during the familiarization phase.
  
#projections> to calculate discrimination indexes (DI) defined as (Tn-Tf)/(Tn+Tf) and/or to compare the novel object preference (%) = ((Tn)/(Tn+Tf))*100
# Tn: time spent with the novel object
# Tf: time spent with the familiar object
  
#hypothesis N3: In the familiarization stage animals have a discrimination index (Di) close to zero independent of age, Segment of the test and BPA exposure ----
  
dat_fam_di <- data_fam %>% 
    ungroup() %>% 
    group_by(ID,`Segment of test`,BPA_EXPOSURE) %>% 
    mutate(Di =(`Familiar Object 1 : time investigating (s)` - `Familiar Object 2 : time investigating (s)`)/(`Familiar Object 1 : time investigating (s)` + `Familiar Object 2 : time investigating (s)`))
  
  data_summary_fam_di <- dat_fam_di %>% 
    ungroup() %>% 
    group_by(`Segment of test`,age_week_round,BPA_EXPOSURE) %>% 
    summarise(
      n = n(),
      mean_time_Di = mean(Di),
      sem_time_Di  = sd(Di)/sqrt(n),
      .groups = "drop"
    ) %>% 
    ungroup()
  
  data_raw_long_fam_di <- dat_fam_di %>% 
    pivot_longer(
      cols = c(Di),
      names_to = "Object",
      values_to = "Di_value"
    )
  
  ggplot(data_summary_fam_di, aes(x = `Segment of test`, y = mean_time_Di)) +
    geom_col(alpha = 0.6, fill = "grey70") +
    geom_errorbar(aes(ymin = mean_time_Di - sem_time_Di,
                      ymax = mean_time_Di + sem_time_Di),
                  width = 0.2) +
    geom_jitter(data = dat_fam_di,
                aes(x = `Segment of test`, y = Di),
                width = 0.1, size = 2, alpha = 0.8) +
    labs(x = "Segment of test",
         y = "Discrimination index (Di)") +
 facet_wrap(~age_week_round*BPA_EXPOSURE, ncol = 2, nrow = 3) +
    theme_classic()
  
  #####STATS----
  
  t_test_di <- dat_fam_di %>%
    group_by(`Segment of test`, age_week_round, BPA_EXPOSURE) %>%
    summarise(
      t_test = list(t.test(Di, mu = 0)),
      .groups = "drop"
    ) %>%
    mutate(t_test = map(t_test, tidy)) %>%
    unnest(t_test)
  
  t_test_di
  
#hypothesis N4: In the recognition stage animals exposed to BPA will decrease Novel object preference % (NOP %) ----
#NOP (%) = ((TN)/(TN+TF))*100 where TN = time spent interacting with novel object
#TF = Time spent interacting with familiar object. 
#Reductions in preference is indicative of impaired learning and memory.

#I took the 0-300 sec because recognition phase for data-2.csv was cut short
  
dat_recog <- all_data%>% 
    ungroup() %>% 
    group_by(`Segment of test`) %>% 
    filter(Stage =='Recognition') %>% 
    mutate(NOP = ((`Novel Object : time investigating (s)`) /
                    (`Familiar Object 1 : time investigating (s)` + 
                       `Novel Object : time investigating (s)`)) * 100) %>% 
    select(ID, `Segment of test`, BPA_EXPOSURE, SEX, COHORT, 
           age_week_round, NOP) %>% 
    filter(SEX=="F")
  
data_summary_recog <- dat_recog %>% 
    ungroup() %>% 
    group_by(`Segment of test`,age_week_round,BPA_EXPOSURE) %>% 
    summarise(
      n = n(),
      mean_time_NOP = mean(NOP),
      sem_time_NOP  = sd(NOP)/sqrt(n),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
  filter(`Segment of test` == "0 - 300 secs.") #because this person recorded just the first part of the test and not the entire test =/

data_raw_long_recog <- dat_recog %>% 
  pivot_longer(
    cols = c(NOP),
    names_to = "Object",
    values_to = "NOP_value"
  ) %>% 
  filter(`Segment of test` == "0 - 300 secs.") #because this person recorded just the first part of the test and not the entire test =/
  

ggplot(data_summary_recog, aes(x = BPA_EXPOSURE, y = mean_time_NOP)) +
  geom_col(alpha = 0.6, fill = "grey70") +
  geom_errorbar(aes(ymin = mean_time_NOP - sem_time_NOP,
                    ymax = mean_time_NOP + sem_time_NOP),
                width = 0.2) +
  geom_jitter(data = dat_recog,
              aes(x = BPA_EXPOSURE, y = NOP),
              width = 0.1, size = 2, alpha = 0.8) +
  labs(x = "BPA perinatal exposure",
       y = "% of novel object preference  (NOP%)") +
  facet_wrap(~age_week_round) +
  theme_classic()

t_test_NOP <- dat_recog %>%
  group_by(age_week_round) %>%
  summarise(t_test = list(t.test(NOP ~ BPA_EXPOSURE)),
            .groups = "drop") %>%
  mutate(t_test = map(t_test, broom::tidy)) %>%
  unnest(t_test)  
t_test_NOP 

##match bw baseline----


