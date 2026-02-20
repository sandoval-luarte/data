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

##-- stats of BPA on dams----
# List of parameters
parameters <- c("gestation_lenght_days", "Liter_size",
                "female_pups_number", "male_pups_number", "F_M_sexratio")

# Function to run appropriate test comparing YES vs NO
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

## BW diet collapsed by diet----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

BW_data_collapsed <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(15, 16,18)) %>% 
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
    -DIET_FORMULA.y
  ) %>% 
  rename(
    SEX = SEX.x,
    DIET_FORMULA = DIET_FORMULA.x
  )

BW_data_collapsed  %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 


BW_summary_collapsed <- BW_data_collapsed %>%
  group_by(day_rel,BPA_EXPOSURE,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## STATS females in OD (HCD+HFD) ----
bw_fem_collapsed <- BW_data_collapsed %>%
  filter(
    SEX == "F",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_fem_collapsed)   # should be > 0
table(bw_fem_collapsed$BPA_EXPOSURE)

fem_bw_collapsed <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_fem_collapsed
)

anova(fem_bw_collapsed)

emm_day_fem_collapsed <- emmeans(
  fem_bw_collapsed,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_fem_collapsed$day_rel)))
)

day_contrasts_fem_collapsed <- contrast(
  emm_day_fem_collapsed,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_fem_collapsed <- as.data.frame(day_contrasts_fem_collapsed) 
day_stats_fem_collapsed 
#so for females in OD diets (HCD + HFD) after day 77 BPA females are heavier than non BPA females (week 11)

## STATS males in OD (HCD+HFD)----
bw_m_collapsed <- BW_data_collapsed %>%
  filter(
    SEX == "M",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_m_collapsed)   # should be > 0
table(bw_m_collapsed$BPA_EXPOSURE)

m_bw_collapsed <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_m_collapsed
)

anova(m_bw_collapsed)

emm_day_m_collapsed <- emmeans(
  m_bw_collapsed,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_m_collapsed$day_rel)))
)

day_contrasts_m_collapsed <- contrast(
  emm_day_m_collapsed,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_m_collapsed <- as.data.frame(day_contrasts_m_collapsed) 
day_stats_m_collapsed
#so for males in OD diets (HCD + HFD) there is none day in which BPA exposed males are heavier than non BPA males

## Figure 1A BW collapsed by diet----

BW_data_collapsed <- BW_data_collapsed %>%
  mutate(
    week_rel = day_rel / 7
  )
BW_data_collapsed <- BW_data_collapsed %>%
  mutate(
    week_rel = floor(day_rel / 7)
  )
BW_summary_collapsed <- BW_data_collapsed %>%
  group_by(week_rel, BPA_EXPOSURE, SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
shade_df <- tibble(
  SEX = c("F"),
  xmin = c(11),
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)

plot_bw_sex_collapsed <- ggplot(
  BW_summary_collapsed,
  aes(
    x = week_rel,
    y = mean_BW,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.35
  ) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW
    ),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(
    ~ SEX) +
  labs(
    x = "Weeks",
    y = "Body weight (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

plot_bw_sex_collapsed

## BW separated by diet ----
METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(15, 16)) %>% 
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
    -DIET_FORMULA.y
  ) %>% 
  rename(
    SEX = SEX.x,
    DIET_FORMULA = DIET_FORMULA.x
  )

BW_data  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 


BW_summary <- BW_data %>%
  group_by(day_rel,BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## STATS females in HCD----
bw_fem_hcd <- BW_data %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_fem_hcd)   # should be > 0
table(bw_fem_hcd$BPA_EXPOSURE)

fem_bw_hcd <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_fem_hcd
)

anova(fem_bw_hcd)

emm_day_fem_hcd <- emmeans(
  fem_bw_hcd,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_fem_hcd$day_rel)))
)

day_contrasts_fem_hcd <- contrast(
  emm_day_fem_hcd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_fem_hcd <- as.data.frame(day_contrasts_fem_hcd) 
day_stats_fem_hcd
#so for females in HCD after day 119 (week 17) BPA females are heavier than non BPA females

## STATS females in HFD----

bw_fem_hfd <- BW_data %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_fem_hfd)   # should be > 0
table(bw_fem_hfd$BPA_EXPOSURE)

fem_bw_hfd <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_fem_hfd
)

anova(fem_bw_hfd)

emm_day_fem_hfd <- emmeans(
  fem_bw_hfd,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_fem_hfd$day_rel)))
)

day_contrasts_fem_hfd <- contrast(
  emm_day_fem_hfd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_fem_hfd <- as.data.frame(day_contrasts_fem_hfd) 
day_stats_fem_hfd 
#so for females in HFD after day 35 (week 5) BPA females are heavier than non BPA females

## STATS males in HCD----

bw_male_hcd <- BW_data %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_male_hcd)   # should be > 0
table(bw_male_hcd$BPA_EXPOSURE)

male_bw_hcd <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_male_hcd
)

anova(male_bw_hcd)

emm_day_male_hcd <- emmeans(
  male_bw_hcd,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_male_hcd$day_rel)))
)

day_contrasts_male_hcd <- contrast(
  emm_day_male_hcd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_male_hcd <- as.data.frame(day_contrasts_male_hcd) 
day_stats_male_hcd
#so for there is no days in which BPA males are heavier than non-BPA males

## STATS males in HFD----

bw_male_hfd <- BW_data %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_male_hfd)   # should be > 0
table(bw_male_hfd$BPA_EXPOSURE)

male_bw_hfd <- lmer(
  BW ~ day_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_male_hfd
)

anova(male_bw_hfd)

emm_day_male_hfd <- emmeans(
  male_bw_hfd,
  ~ BPA_EXPOSURE | day_rel,
  at = list(day_rel = sort(unique(bw_male_hfd$day_rel)))
)

day_contrasts_male_hfd <- contrast(
  emm_day_male_hfd,
  method = "pairwise",
  adjust = "fdr"   # multiple testing correction
)

day_stats_male_hfd <- as.data.frame(day_contrasts_male_hfd) 
day_stats_male_hfd 
#so for there is no days in which BPA males are heavier than non-BPA males

## Figure 1B BW separated by diet----

BW_data <- BW_data %>%
  mutate(
    week_rel = day_rel / 7
  )
BW_data <- BW_data %>%
  mutate(
    week_rel = floor(day_rel / 7)
  )
BW_summary <- BW_data %>%
  group_by(week_rel, BPA_EXPOSURE, SEX, DIET_FORMULA) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
shade_df <- tibble(
  SEX = c("F", "F"),
  DIET_FORMULA = c("D12450Hi", "D12451i"),
  xmin = c(17, 5),
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
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
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.35
  ) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW
    ),
    alpha = 0.25,
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
    y = "Body weight (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

plot_bw_sex

# Figure 1 (Figure 1A+ Figure 1B) ----
plot_bw_sex_collapsed <- plot_bw_sex_collapsed + labs(tag = "A")
plot_bw_sex <- plot_bw_sex + labs(tag = "B")

combined_plot <- plot_bw_sex_collapsed | plot_bw_sex
combined_plot

## BW at day 0 with HCD or HFD----

BW_day0_raw <- BW_data %>% 
  filter(bw_rel == 0)

BW_day0_sum <- BW_data %>% 
  filter(bw_rel == 0) %>% 
  group_by(BPA_EXPOSURE, SEX) %>% 
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

#### supplementary figure 1----

ggplot(BW_day0_sum,
       aes(x = BPA_EXPOSURE, y = mean_BW, fill = BPA_EXPOSURE)) +
  
  # Bars (means)
  geom_col(width = 0.6, color = "black", alpha = 0.7) +
  
  # SEM error bars
  geom_errorbar(
    aes(ymin = mean_BW - sem_BW,
        ymax = mean_BW + sem_BW),
    width = 0.2,
    linewidth = 0.8
  ) +
  
  # Individual points
  geom_jitter(
    data = BW_day0_raw,
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
    y = "Body weight (g)",
    x = "BPA exposure"
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

### STATS----
bw_day0 <- BW_data %>% 
  filter(bw_rel == 0)

lm_bw_day0 <- lm(
  BW ~ SEX * BPA_EXPOSURE * DIET_FORMULA,
  data = bw_day0
)
anova(lm_bw_day0)

bw_day0 <- BW_data %>% 
  filter(bw_rel == 0) %>% 
  distinct(ID, .keep_all = TRUE)

# Females
t_female <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_day0 %>% filter(SEX == "F")
)

# Males
t_male <- t.test(
  BW ~ BPA_EXPOSURE,
  data = bw_day0 %>% filter(SEX == "M")
)

t_female
t_male

## BW gain over time data----

BW_gain <- BW_data %>% 
 # filter(SEX =="F") %>% 
  select(ID,day_rel,bw_rel,BPA_EXPOSURE,DIET_FORMULA,SEX) %>% 
  mutate(week_rel = floor(day_rel / 7))

BW_gain %>% 
  group_by(BPA_EXPOSURE,DIET_FORMULA,SEX) %>%
  summarise(n_ID = n_distinct(ID)) 

BW_gainsummary <- BW_gain %>%
  group_by(week_rel, BPA_EXPOSURE, DIET_FORMULA,SEX) %>%
  summarise(
    mean_BWgain = mean(bw_rel, na.rm = TRUE),
    sem_BWgain  = sd(bw_rel, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## STATS----
#### females in HCD----
bw_gain_fem_hcd <- BW_gain %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_gain_fem_hcd)
table(bw_gain_fem_hcd$BPA_EXPOSURE)

fem_gain_hcd <- lmer(
  bw_rel ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_gain_fem_hcd
)

anova(fem_gain_hcd)

emm_week_fem_hcd <- emmeans(
  fem_gain_hcd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_gain_fem_hcd$week_rel)))
)

week_contrasts_fem_hcd <- contrast(
  emm_week_fem_hcd,
  method = "pairwise",
  adjust = "fdr"
)

week_stats_fem_hcd <- as.data.frame(week_contrasts_fem_hcd)
week_stats_fem_hcd 
#There were no weeks in which BPA-exposed females fed HCD exhibited a greater 
#percentage body-weight gain than non-BPA females

#### females in HFD----
bw_gain_fem_hfd <- BW_gain %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_gain_fem_hfd)
table(bw_gain_fem_hfd$BPA_EXPOSURE)

fem_gain_hfd <- lmer(
  bw_rel ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_gain_fem_hfd
)

anova(fem_gain_hfd)

emm_week_fem_hfd <- emmeans(
  fem_gain_hfd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_gain_fem_hfd$week_rel)))
)

week_contrasts_fem_hfd <- contrast(
  emm_week_fem_hfd,
  method = "pairwise",
  adjust = "fdr"
)

week_stats_fem_hfd <- as.data.frame(week_contrasts_fem_hfd)
week_stats_fem_hfd 
#There were no weeks in which BPA-exposed females fed HFD exhibited a greater 
#percentage body-weight gain than non-BPA females

#### males in HCD----
bw_gain_m_hcd <- BW_gain %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_gain_m_hcd)
table(bw_gain_m_hcd$BPA_EXPOSURE)

m_gain_hcd <- lmer(
  bw_rel ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_gain_m_hcd
)

anova(m_gain_hcd)

emm_week_m_hcd <- emmeans(
  m_gain_hcd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_gain_m_hcd$week_rel)))
)

week_contrasts_m_hcd <- contrast(
  emm_week_m_hcd,
  method = "pairwise",
  adjust = "fdr"
)

week_stats_m_hcd <- as.data.frame(week_contrasts_m_hcd)
week_stats_m_hcd 
#There were no weeks in which BPA-exposed males fed HCD exhibited a greater 
#percentage body-weight gain than non-BPA males

#### males in HFD----
bw_gain_m_hfd <- BW_gain %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

nrow(bw_gain_m_hfd)
table(bw_gain_m_hfd$BPA_EXPOSURE)

m_gain_hfd <- lmer(
  bw_rel ~ week_rel * BPA_EXPOSURE + (1 | ID),
  data = bw_gain_m_hfd
)

anova(m_gain_hfd)

emm_week_m_hfd <- emmeans(
  m_gain_hfd,
  ~ BPA_EXPOSURE | week_rel,
  at = list(week_rel = sort(unique(bw_gain_m_hfd$week_rel)))
)

week_contrasts_m_hfd <- contrast(
  emm_week_m_hfd,
  method = "pairwise",
  adjust = "fdr"
)

week_stats_m_hfd <- as.data.frame(week_contrasts_m_hfd)
week_stats_m_hfd 
#There were no weeks in which BPA-exposed males fed HFD exhibited a greater 
#percentage body-weight gain than non-BPA males

### supplementary figure 2 ----

plot_bw_gain <- ggplot(
  BW_gainsummary,
  aes(
    x = week_rel,
    y = mean_BWgain,
    color = BPA_EXPOSURE,
    fill  = BPA_EXPOSURE
  )
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_BWgain - sem_BWgain,
      ymax = mean_BWgain + sem_BWgain
    ),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(
    ~ SEX*DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  ) +
  labs(
    x = "Weeks",
    y = "Body weight gain (%)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

plot_bw_gain

## Speed = rate of change of BW over time----

BW_speed <- BW_data %>%
  group_by(ID, SEX, BPA_EXPOSURE,DIET_FORMULA) %>%
  do(tidy(lm(BW ~ day_rel, data = .))) %>%
  filter(term == "day_rel") %>%
  rename(
    speed_g_per_day = estimate,
    se = std.error
  )

BW_speed %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

# plot 1 BPA within sex 
p1 <- ggplot(BW_speed,
             aes(x = BPA_EXPOSURE,
                 y = speed_g_per_day,
                 fill = BPA_EXPOSURE)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_grid(
    ~ SEX * DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  )+
  labs(
    y = "BW gain speed (g/day)",
    x = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

p1


### STATS----

lmer_speed <- lmer(
  BW ~ day_rel * BPA_EXPOSURE * SEX + (day_rel | ID),
  data = BW_data
)
summary(lmer_speed)


#Male mice exhibited a significantly faster rate of body
#weight gain than females (≈0.08 g/day greater; t = 4.77).
#BPA exposure did not significantly alter growth rate in 
#either sex, and no significant BPA × sex interaction on
#growth speed was observed.

#Although the estimated growth rate was slightly higher in BPA-exposed females,
#this difference did not reach statistical significance (t = 1.36).

### Could BPA affect growth speed only later in life?----

lmer(
  BW ~ day_rel * BPA_EXPOSURE * SEX + (day_rel | ID),
  data = BW_data %>% filter(day_rel >= 35) #the day in which we saw differences in females HFD
)

#Even when restricting the analysis to the post-divergence period (day ≥ 35)
#BPA exposure did not significantly alter the rate of body-weight gain in either sex.

# Plot 2: BPA effects collapsed by sex 
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
  facet_wrap(
    ~  DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD"
      )
    )
  )+
  theme_classic(base_size = 14)

p2

###supplementary figure 3----

p1 <- p1 + labs(tag = "A")
p2 <- p2 + labs(tag = "B")

combined_plot <- p1 | p2
combined_plot

#so Sex was the primary determinant of growth rate
#with females growing more slowly than males, regardless of BPA exposure.
#However, BPA exposure was associated with higher BW in females, specially after HFD exposure
#without a corresponding increase in growth rate.

#CONCLUSIONS (within each sex BPA exposed vs non BPA exposed)
# Same starting BW
# Same % BW gain
# Same growth rate
# Higher final BW in BPA females in HCD and HFD but stronger in HFD

# BODY COMPOSITION ANALYSIS----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16,18)) %>% 
  filter(!ID %in% c("9406", "9354")) %>%  #9406 has a  weird pattern in locomotion and 9354 lack body comp measurements for 18 weeks
group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index,COHORT,DIET_FORMULA) %>% 
  left_join(METABPA, by= "ID") %>% 
  ungroup() %>% 
  select(
    -DIET_FORMULA.y) %>% 
  rename(
    DIET_FORMULA = DIET_FORMULA.x)

## Adiposity index collapsed by diet ----

echoMRI_data %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have the same amount of data than for BW (i.e 18 F and 20 M)

echoMRI_data_comparisons_collapsed <- echoMRI_data %>% 
  mutate(
    n_measurement= case_when(
      COHORT == 15 & Date == "2025-04-29" ~ "0 wks",
      COHORT == 16 & Date == "2025-08-04" ~ "0 wks",
      COHORT == 18 & Date %in% c("2026-01-21", "2026-01-23") ~ "0 wks",
      COHORT == 15 & Date == "2025-05-28" ~ "4 wks",
      COHORT == 16 & Date == "2025-09-09" ~ "4 wks", #really this is 5 wks
      COHORT == 15 & Date == "2025-07-07" ~ "10 wks",
      COHORT == 16 & Date == "2025-10-07" ~ "10 wks", #really 9 wks
      COHORT == 15 & Date == "2025-08-01" ~ "13 wks",
      COHORT == 16 & Date == "2025-11-17" ~ "13 wks", #really this is 15 wks
      COHORT == 15 & Date == "2025-09-09" ~ "19 wks",
      COHORT == 16 & Date == "2025-12-18" ~ "19 wks",
      COHORT == 15 & Date == "2025-10-07" ~ "23 wks",
      COHORT == 16 & Date == "2026-01-09" ~ "23 wks" #really this is 22 wks
    )) %>% 
  mutate(
    n_measurement = factor(
      n_measurement,
      levels = c("0 wks", "4 wks", "10 wks", "13 wks", "19 wks", "23 wks")
    )
  )


echoMRI_data_comparisons_collapsed  %>% 
  group_by(SEX,n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have the same amount of data than for BW (i.e 17 F and 20 M) we removed ID 9406 and 9354 

bodycomp_summary_collapsed <- echoMRI_data_comparisons_collapsed %>%
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
  mutate(n_measurement = factor(n_measurement, levels = c("0 wks",
                                                          "3 wks",
                                                          #"5 wks",
                                                          "9 wks",
                                                          "12 wks",
                                                          #"15 wks",
                                                          "18 wks",
                                                          #   "19 wks",
                                                          "22 wks"))) %>% 
  ungroup()

### supplementary figure 4A (SF4A) adiposity index collapsed by diet ----
sf4a <- ggplot(bodycomp_summary_collapsed,
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
  facet_wrap( ~ SEX)+
  labs(
    x = "weeks (wks)",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf4a

####STATS for females----

ai_fem <- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "F",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_fem$BPA_EXPOSURE)

fem_ai<- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_fem
)

anova(fem_ai)

emm_fem_ai<- emmeans(
  fem_ai,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_ai <- contrast(
  emm_fem_ai,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_fem<- as.data.frame(contr_fem_ai)
ai_stats_fem

#there is no days in which females exposed to BPA have higher AI than control ones

####STATS for males----

ai_m <- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "M",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_m$BPA_EXPOSURE)

m_ai<- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_m
)

anova(m_ai)

emm_m_ai<- emmeans(
  m_ai,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_ai <- contrast(
  emm_m_ai,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_m<- as.data.frame(contr_m_ai)
ai_stats_m

#there is no days in which males exposed to BPA have higher AI than control ones

## Adiposity index separated by diet ---- 

echoMRI_data %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have the same amount of data than for BW (i.e 18 F and 20 M)

echoMRI_data_comparisons <- echoMRI_data %>% 
  mutate(
    n_measurement= case_when(
      COHORT == 15 & Date == "2025-04-29" ~ "0 wks",
      COHORT == 16 & Date == "2025-08-04" ~ "0 wks",
      COHORT == 15 & Date == "2025-05-28" ~ "3 wks",
      COHORT == 16 & Date == "2025-09-09" ~ "3 wks", #really this is 5 wks
      COHORT == 15 & Date == "2025-07-07" ~ "9 wks",
      COHORT == 16 & Date == "2025-10-07" ~ "9 wks",
      COHORT == 15 & Date == "2025-08-01" ~ "12 wks",
      COHORT == 16 & Date == "2025-11-17" ~ "12 wks", #really this is 15 wks
      COHORT == 15 & Date == "2025-09-09" ~ "18 wks",
      COHORT == 16 & Date == "2025-12-18" ~ "18 wks",#really this is 19 wks
      COHORT == 15 & Date == "2025-10-07" ~ "22 wks",
      COHORT == 16 & Date == "2026-01-09" ~ "22 wks"
     )) %>% 
  mutate(
    n_measurement = factor(
      n_measurement,
      levels = c("0 wks", "3 wks", "9 wks", "12 wks", "18 wks", "22 wks")
    )
  )


echoMRI_data_comparisons  %>% 
  group_by(SEX,n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have the same amount of data than for BW (i.e 17 F and 20 M) we removed ID 9406 and 9354 

bodycomp_summary <- echoMRI_data_comparisons %>%
  group_by(n_measurement, BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
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
  mutate(n_measurement = factor(n_measurement, levels = c("0 wks",
                                                          "3 wks",
                                                         #"5 wks",
                                                          "9 wks",
                                                          "12 wks",
                                                         #"15 wks",
                                                          "18 wks",
                                                       #   "19 wks",
                                                          "22 wks"))) %>% 
  ungroup()


## supplementary figure 4B (SF4B) adiposity index separated by diet ----
sf4b <- ggplot(bodycomp_summary,
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
  facet_wrap( ~ SEX*DIET_FORMULA,
    labeller = labeller(
      DIET_FORMULA = c(
        "D12450Hi" = "HCD",
        "D12451i"  = "HFD" ))) +
  labs(
    x = "weeks (wks)",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf4b

###STATS for females HCD----

ai_fem_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_fem_hcd$BPA_EXPOSURE)

fem_ai_hcd <- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_fem_hcd
)

anova(fem_ai_hcd)

emm_fem_ai_hcd <- emmeans(
  fem_ai_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_ai_hcd <- contrast(
  emm_fem_ai_hcd,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_fem_hcd <- as.data.frame(contr_fem_ai_hcd)
ai_stats_fem_hcd

#there is no days in which females in HCD exposed to BPA have higher AI than control ones

###STATS for females HFD----

ai_fem_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_fem_hfd$BPA_EXPOSURE)

fem_ai_hfd <- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_fem_hfd
)

anova(fem_ai_hfd)

emm_fem_ai_hfd <- emmeans(
  fem_ai_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_ai_hfd <- contrast(
  emm_fem_ai_hfd,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_fem_hfd <- as.data.frame(contr_fem_ai_hfd)
ai_stats_fem_hfd

#there is no days in which females in HFD exposed to BPA have higher AI than control ones

###STATS for males HCD----

ai_m_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_m_hcd$BPA_EXPOSURE)

m_ai_hcd <- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_m_hcd
)

anova(m_ai_hcd)

emm_m_ai_hcd <- emmeans(
  m_ai_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_ai_hcd <- contrast(
  emm_m_ai_hcd,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_m_hcd <- as.data.frame(contr_m_ai_hcd)
ai_stats_m_hcd

#there is no days in which males in HCD exposed to BPA have higher AI than control ones

###STATS for males HFD----

ai_m_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(ai_m_hfd$BPA_EXPOSURE)

m_ai_hfd <- lmer(
  adiposity_index ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = ai_m_hfd
)

anova(m_ai_hfd)

emm_m_ai_hfd <- emmeans(
  m_ai_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_ai_hfd <- contrast(
  emm_m_ai_hfd,
  method = "pairwise",
  adjust = "fdr"
)

ai_stats_m_hfd <- as.data.frame(contr_m_ai_hfd)
ai_stats_m_hfd

#there is no days in which males in HFD exposed to BPA have higher AI than control ones


# supplementary figure 4 (SF4A+ SF4B)----
sf4a <- sf4a + labs(tag = "A")
sf4b <- sf4b + labs(tag = "B")

combined_plot_ai <- sf4a  | sf4b 
combined_plot_ai

## supplementary figure 5A (SF5A) fat mass collapsed by diet ----
sf5a <- ggplot(bodycomp_summary_collapsed,
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
  facet_wrap( ~ SEX) +
  labs(
    x = "weeks (wks)",
    y = "Fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf5a

###STATS for females----

fat_fem<- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "F",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_fem$BPA_EXPOSURE)

fem_fat<- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_fem
)

anova(fem_fat)

emm_fem_fat<- emmeans(
  fem_fat,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_fat<- contrast(
  emm_fem_fat,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_fem<- as.data.frame(contr_fem_fat)
fat_stats_fem

#there is no days in which females exposed to BPA have higher fat mass than control ones

###STATS for males----
fat_m <- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "M",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_m$BPA_EXPOSURE)

m_fat <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_m
)

anova(m_fat)

emm_m_fat <- emmeans(
  m_fat,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_fat<- contrast(
  emm_m_fat,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_m <- as.data.frame(contr_m_fat)
fat_stats_m

#there is no days in which males exposed to BPA have higher fat mass than control ones

## supplementary figure 5B (SF5B) fat mass splitted by diet ----
  
sf5b <- ggplot(bodycomp_summary,
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
  facet_wrap( ~ SEX*DIET_FORMULA,
              labeller = labeller(
                DIET_FORMULA = c(
                  "D12450Hi" = "HCD",
                  "D12451i"  = "HFD" ))) +
  labs(
    x = "weeks (wks)",
    y = "Fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf5b 

###STATS for females HCD----

fat_fem_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_fem_hcd$BPA_EXPOSURE)

fem_fat_hcd <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_fem_hcd
)

anova(fem_fat_hcd)

emm_fem_fat_hcd <- emmeans(
  fem_fat_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_fat_hcd <- contrast(
  emm_fem_fat_hcd,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_fem_hcd <- as.data.frame(contr_fem_fat_hcd)
fat_stats_fem_hcd

#there is no days in which females in HCD exposed to BPA have higher fat mass than control ones

###STATS for females HFD----

fat_fem_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_fem_hfd$BPA_EXPOSURE)

fem_fat_hfd <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_fem_hfd
)

anova(fem_fat_hfd)

emm_fem_fat_hfd <- emmeans(
  fem_fat_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_fat_hfd <- contrast(
  emm_fem_fat_hfd,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_fem_hfd <- as.data.frame(contr_fem_fat_hfd)
fat_stats_fem_hfd

#there is no days in which females in HFD exposed to BPA have higher fat mass than control ones

###STATS for males HCD----

fat_m_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_m_hcd$BPA_EXPOSURE)

m_fat_hcd <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_m_hcd
)

anova(m_fat_hcd)

emm_m_fat_hcd <- emmeans(
  m_fat_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_fat_hcd <- contrast(
  emm_m_fat_hcd,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_m_hcd <- as.data.frame(contr_m_fat_hcd)
fat_stats_m_hcd

#there is no days in which females in HCD exposed to BPA have higher fat mass than control ones

###STATS for males HFD----

fat_m_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(fat_m_hfd$BPA_EXPOSURE)

m_fat_hfd <- lmer(
  Fat ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = fat_m_hfd
)

anova(m_fat_hfd)

emm_m_fat_hfd <- emmeans(
  m_fat_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_fat_hfd <- contrast(
  emm_m_fat_hfd,
  method = "pairwise",
  adjust = "fdr"
)

fat_stats_m_hfd <- as.data.frame(contr_m_fat_hfd)
fat_stats_m_hfd

#there is no days in which males in HFD exposed to BPA have higher fat mass than control ones

# supplementary figure 5 (SF4A+ SF4B)----
sf5a <- sf5a + labs(tag = "A")
sf5b <- sf5b + labs(tag = "B")

combined_plot_fm <- sf5a  | sf5b 
combined_plot_fm

## Figure 2A (F2A) lean mass collapsed by diet ----
f2a <- ggplot(bodycomp_summary_collapsed,
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
  facet_wrap( ~ SEX) +
  labs(
    x = "Days of measurement",
    y = "Lean mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
f2a

####STATS for females----

lean_fem<- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "F",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_fem$BPA_EXPOSURE)

fem_lean<- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_fem
)

anova(fem_lean)

emm_fem_lean<- emmeans(
  fem_lean,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_lean <- contrast(
  emm_fem_lean,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_fem<- as.data.frame(contr_fem_lean)
lean_stats_fem

# so after week 12 female mice exposed to BPA have more lean mass than non exposed one independent of the obesogenic diet

####STATS for males----

lean_m<- echoMRI_data_comparisons_collapsed %>%
  filter(
    SEX == "M",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_m$BPA_EXPOSURE)

m_lean<- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_m
)

anova(m_lean)

emm_m_lean<- emmeans(
  m_lean,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_lean <- contrast(
  emm_m_lean,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_m<- as.data.frame(contr_m_lean)
lean_stats_m

#there is no weeks in which males exposed to BPA have more lean mass than non exposed ones

## Figure 2B (F2B) lean mass by diet ----
f2b <- ggplot(bodycomp_summary,
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
  facet_wrap( ~ SEX*DIET_FORMULA,
              labeller = labeller(
                DIET_FORMULA = c(
                  "D12450Hi" = "HCD",
                  "D12451i"  = "HFD" ))) +
  labs(
    x = "Days of measurement",
    y = "Lean mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
f2b

###STATS for females HCD----

lean_fem_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_fem_hcd$BPA_EXPOSURE)

fem_lean_hcd <- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_fem_hcd
)

anova(fem_lean_hcd)

emm_fem_lean_hcd <- emmeans(
  fem_lean_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_lean_hcd <- contrast(
  emm_fem_lean_hcd,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_fem_hcd <- as.data.frame(contr_fem_lean_hcd)
lean_stats_fem_hcd

#after 18 weeks of HCD onwards BPA exposed females have more lean mass than control ones


###STATS for females HFD----

lean_fem_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "F",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_fem_hfd$BPA_EXPOSURE)

fem_lean_hfd <- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_fem_hfd
)

anova(fem_lean_hfd)

emm_fem_lean_hfd <- emmeans(
  fem_lean_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_fem_lean_hfd <- contrast(
  emm_fem_lean_hfd,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_fem_hfd <- as.data.frame(contr_fem_lean_hfd)
lean_stats_fem_hfd

#after 12 weeks of HFD onwards BPA exposed females have more lean mass than control ones

###STATS for males HCD----

lean_m_hcd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12450Hi",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_m_hcd$BPA_EXPOSURE)

m_lean_hcd <- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_m_hcd
)

anova(m_lean_hcd)

emm_m_lean_hcd <- emmeans(
  m_lean_hcd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_lean_hcd <- contrast(
  emm_m_lean_hcd,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_m_hcd <- as.data.frame(contr_m_lean_hcd)
lean_stats_m_hcd

#there is no weeks in which BPA exposed males have more lean mass than non exposed one after being fed with HCD

###STATS for males HFD----

lean_m_hfd <- echoMRI_data_comparisons %>%
  filter(
    SEX == "M",
    DIET_FORMULA == "D12451i",
    BPA_EXPOSURE %in% c("YES", "NO")
  )

table(lean_m_hfd$BPA_EXPOSURE)

m_lean_hfd <- lmer(
  Lean ~ n_measurement * BPA_EXPOSURE + (1 | ID),
  data = lean_m_hfd
)

anova(m_lean_hfd)

emm_m_lean_hfd <- emmeans(
  m_lean_hfd,
  ~ BPA_EXPOSURE | n_measurement
)

contr_m_lean_hfd <- contrast(
  emm_m_lean_hfd,
  method = "pairwise",
  adjust = "fdr"
)

lean_stats_m_hfd <- as.data.frame(contr_m_lean_hfd)
lean_stats_m_hfd

#there is no weeks in which BPA exposed males have more lean mass than non exposed one after being fed with HFD

# Figure 2 (F2A+ F2B)----
f2a <- f2a + labs(tag = "A")
f2b <- f2b + labs(tag = "B")

combined_plot_lm <- f2a  | f2b 
combined_plot_lm


# FOOD INTAKE ANALYSIS----

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
  ungroup() %>% 
  filter(!ID ==9406)  #9406 has a  weird pattern in locomotion

FI_data   %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

####FI collapsed by diet----

FI_summary <- FI_data %>%
  group_by(day_rel, SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_FI = mean(FIcumulative, na.rm = TRUE),
    sem_FI  = sd(FIcumulative, na.rm = TRUE) / sqrt(n_distinct(ID)),
    n = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  filter(day_rel == 138)

FI_plotA <- ggplot(FI_summary, aes(x = BPA_EXPOSURE, y = mean_FI, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_FI - sem_FI, ymax = mean_FI + sem_FI),
    width = 0.2
  ) +
  facet_wrap(~ SEX) +
  labs(
    y = "Cumulative food intake (kcal)",
    x = "BPA exposure",
    fill = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

FI_plotA

### STATS -----

FI_day138 <- FI_data %>%
  filter(day_rel == 138) %>%
  group_by(ID, SEX, BPA_EXPOSURE) %>%
  summarise(
    FIcumulative = mean(FIcumulative, na.rm = TRUE),
    .groups = "drop"
  )

###### STATS for females----
leveneTest(
  FIcumulative ~ BPA_EXPOSURE,
  data = FI_day138 %>% filter(SEX == "F")
) #p> 0.05, no evidence that variances differ between BPA vs control

t_female <- t.test(
  FIcumulative ~ BPA_EXPOSURE,
  data = FI_day138 %>% filter(SEX == "F"),
  var.equal = TRUE
)
t_female

###### STATS for males----
leveneTest(
  FIcumulative ~ BPA_EXPOSURE,
  data = FI_day138 %>% filter(SEX == "M")
) #p> 0.05, no evidence that variances differ between BPA vs control

t_male <- t.test(
  FIcumulative ~ BPA_EXPOSURE,
  data = FI_day138 %>% filter(SEX == "M"),
  var.equal = TRUE
)
t_male

####FI separated by diet----

FI_data  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID))

FI_summary <- FI_data %>%
  group_by(day_rel, SEX, BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(
    mean_FI = mean(FIcumulative, na.rm = TRUE),
    sem_FI  = sd(FIcumulative, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>% 
  filter(day_rel == 138)

FI_plotB <- ggplot(FI_summary, aes(x = BPA_EXPOSURE, y = mean_FI, fill = BPA_EXPOSURE)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_FI - sem_FI, ymax = mean_FI + sem_FI),
    width = 0.2
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
    y = "Cumulative food intake (kcal)",
    x = "BPA exposure",
    fill = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

FI_plotB

####FI separated sex and collapsed by diet----
##### just to see if males ate more than females

FI_data  %>% 
  group_by(SEX) %>%
  summarise(n_ID = n_distinct(ID))

FI_summary <- FI_data %>%
  group_by(day_rel, SEX) %>%
  summarise(
    mean_FI = mean(FIcumulative, na.rm = TRUE),
    sem_FI  = sd(FIcumulative, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>% 
  filter(day_rel == 138)
FI_summary

# Filter data for day 138 and keep SEX column
FI_day138 <- FI_data %>%
  filter(day_rel == 138) %>%
  group_by(ID, SEX) %>%
  summarise(FIcumulative = mean(FIcumulative, na.rm = TRUE), .groups = "drop")

# Check group sizes
table(FI_day138$SEX)

# Levene's test for equality of variances (optional but recommended)

leveneTest(FIcumulative ~ SEX, data = FI_day138)

# Run unpaired t-test (assume equal variances if Levene p > 0.05)
t_test_sex <- t.test(FIcumulative ~ SEX,
                     data = FI_day138,
                     var.equal = TRUE)  # set FALSE if variances unequal

# View results
t_test_sex


FI_plotC <- ggplot(FI_summary, aes(x = SEX, y = mean_FI)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_FI - sem_FI, ymax = mean_FI + sem_FI),
    width = 0.2
  ) +
  labs(
    y = "Cumulative food intake (kcal)",
    x = "SEX"
  ) +
  theme_classic(base_size = 14)

FI_plotC

# supplementary figure 6 (FIplotA + FIplotB + FIplotC)----
sf6a <- FI_plotC + labs(tag = "A")
sf6b <- FI_plotA  + labs(tag = "B")
sf6c <- FI_plotB  + labs(tag = "C")

combined_plot_fi <- sf6a  | sf6b | sf6c 
combined_plot_fi


# OGTT----

#OGTT WAS DONE FOR COHORT 15 AND 16 WHEN ALL ANIMALS WERE 27 WEEK OLD

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

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

####OGTT collapsed by diet----

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

ogtta <- ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
 facet_wrap(~ SEX) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()
ogtta

##### STATS for females----

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
  
####OGTT collapsed by diet----

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

sf7a <- ogtta + labs(tag = "A")
sf7b <- ogttc + labs(tag = "B")
sf7c <- ogttb + labs(tag = "C")

combined_plot_ogtt <- sf7a  | sf7b | sf7c 
combined_plot_ogtt

# INDIRECT CALORIMETRY / COLUMBUS DATA ANALYSIS ----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

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

ical_long_all <- bind_rows(ical_long15, ical_long16) #this is key, here we combined

ical_long_all %>% 
filter(!ID ==9406) %>%   #9406 has a  weird pattern in locomotion
group_by(BPA_EXPOSURE,SEX) %>%
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
  )

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
  )

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
  mutate(
    # Positions to display p-values on plot
    y_position = max(relative_count_at_19$relative_total_count_19) * 1.05
  ) %>%
  select(SEX, estimate1, estimate2, statistic, p.value, y_position)

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
  # Add p-values from t-test
  geom_text(
    data = t_test_summaryheat,
    aes(x = 1.5,  # midpoint between NO (1) and YES (2)
        y = y_position,
        label = paste0("p = ", signif(p.value, 3))),
    inherit.aes = FALSE
  ) +
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

#Read me: 
#Data.csv- raw data for cohort 1 at 4 months of age 
#Data-3.csv- Raw data for cohort 1 at 6 months of age and cohort 2 at 2 months
#Data-2.csv- Raw data for cohort 2 at 4 months of age

#so in conclusion the only data we can mix is data.csv and data-2.csv

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
library(broom)  # for tidy output

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


#projections> to calculate discrimination indexes defined as (Tn-Tf)/(Tn+Tf) and/or to compare the novel object preference (%) = ((Tn)/(Tn+Tf))*100
# Tn: time spent with the novel object
# Tf: time spent with the familiar object