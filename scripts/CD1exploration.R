# This script aims to explore changes in body weight (BW), body composition and behavior in females and males CD1 offspring
#exposed to perinatal BPA 50 ug/kg and fed with HFD (D12451i, research diets) or HCD (D12451Hi, research diets)

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
library(ggpattern)
library(car)

# BW diet collapsed----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

BW_data_collapsed <- read_csv("../data/BW.csv") %>% 
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

BW_data_collapsed  %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 


BW_summary_collapsed <- BW_data_collapsed %>%
  group_by(day_rel,BPA_EXPOSURE,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## STATS----
## females in OD (HCD+HFD) ----
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
#so for females in OD diets (HCD + HFD) after day 77 BPA females are heavier than non BPA females

## males in OD (HCD+HFD)----
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
#so for males in OD diets (HCD + HFD) there is none day in which BPA exposed males are heavier than non BPA males

# paper figure 1 plot collapsed by diet----

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

# BW over time data by diet ----
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
##STATS----

## females in HCD----
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
#so for females in HCD after day 119 BPA females are heavier than non BPA females

## females in HFD----

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
#so for females in HFD after day 35 BPA females are heavier than non BPA females

## males in HCD----

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
#so for there is no days in which BPA males are heavier than non-BPA males

## males in HFD----

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
#so for there is no days in which BPA males are heavier than non-BPA males

# paper figure 1 plot----

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

# figure 1 paper final----
plot_bw_sex_collapsed <- plot_bw_sex_collapsed + labs(tag = "A")
plot_bw_sex <- plot_bw_sex + labs(tag = "B")

combined_plot <- plot_bw_sex_collapsed | plot_bw_sex
combined_plot

# BW at day 0 with HCD or HFD----

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

# BW gain over time data----

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
#There were no weeks in which BPA-exposed males fed HFD exhibited a greater 
#percentage body-weight gain than non-BPA males

# paper supp figure 1 plot----

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

# Speed = rate of change of body weight over time data----

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


## STATS----

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

## Could BPA affect growth speed only later in life?----

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

# paper supp figure 2 plot----

p1 <- p1 + labs(tag = "B")
p2 <- p2 + labs(tag = "A")

combined_plot <- p2 | p1
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

# Body comp over time ----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16)) %>% 
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

## AI diet collapsed ----

echoMRI_data %>% 
  group_by(SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great we have the same amount of data than for BW (i.e 18 F and 20 M)

echoMRI_data_comparisons_collapsed <- echoMRI_data %>% 
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

## supplementary plot 3 AI collapsed by diet ----
sp3 <- ggplot(bodycomp_summary_collapsed,
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
    x = "Days of measurement",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sp3

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

## AI separated by diet ---- 

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


## supplementary plot 4? AI splitted by diet ----
sp4 <- ggplot(bodycomp_summary,
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
    x = "Days of measurement",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sp4

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

## supplementary plot 5? FM collapsed by diet ----
sp5 <- ggplot(bodycomp_summary_collapsed,
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
    x = "Days of measurement",
    y = "Fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sp5

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


----
  
  p2 <- ggplot(bodycomp_summary_collapsed,
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
    x = "Days of measurement",
    y = "Fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2


## plot 3 LEAN MASS males and females ----
p3 <- ggplot(bodycomp_summary,
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
p3

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



# OGTT----

#OGTT WAS DONE FOR COHORT 15 AND 16 WHEN ALL ANIMALS WERE 27 WEEK OLD

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

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
  group_by(ID, BPA_EXPOSURE,SEX,DIET_FORMULA) %>% 
  summarise(
    AUC = trapz(time_min, glucose),
    .groups = "drop"
  )

ggplot(auc_df, aes(x = BPA_EXPOSURE, y = AUC)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
 facet_wrap(~ SEX*DIET_FORMULA) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()

## STATS----
###females in HCD----

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

###females in HFD----

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

###males in HCD----

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

###males in HFD----

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

### Does BPA exposure affect OGTT AUC differently in females vs males, regardless of diet?----

auc_sex_bpa <- auc_df %>%
  filter(BPA_EXPOSURE %in% c("YES", "NO"))

auc_sex_bpa %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  summarise(n = n(), .groups = "drop")

model <- lm(AUC ~ SEX * BPA_EXPOSURE, data = auc_sex_bpa)

summary(model)
anova(model)

#A significant main effect of sex was observed, with males exhibiting higher AUC values than females
#independent of BPA exposure


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

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

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
  filter(!ID %in% c(9354, 9414,9406)) %>%  # 9367 and 9368 have confused data from 8/1/25 echoMRI 
  # 9354 lack of complete schedule of measurements in echoMRI
  # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
  # 9406 weird pattern in locomotion
  group_by(BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
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
  filter(!ID %in% c(9354, 9367, 9368, 9414,9406)) %>%  # 9367 and 9368 have confused data from 8/1/25 echoMRI 
  # 9354 lack of complete schedule of measurements in echoMRI
  # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
  # 9406 weird pattern in locomotion
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
  filter(!ID %in% c(9367, 9406)) %>% #check if something weird happened with these animals during the data collection, check with ZR
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
  facet_wrap(~ SEX*DIET_FORMULA) +
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

# Summarize relative_total_count at hour == 19 including SEX and BPA_EXPOSURE
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

# Run unpaired t-test within each sex

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


ggplot(summary_19, aes(x = BPA_EXPOSURE, y = mean_count, fill = BPA_EXPOSURE)) +
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
    y = "Total Counts over 24h",
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

## stats ----

fem_data <- ical_long_all %>%
  filter(
    SEX == "F",
  !ID %in% c(9354, 9414,9406)) %>%
      # 9354 lack of complete schedule of measurements in echoMRI
      # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
      # 9406 weird pattern in locomotion
  mutate(
    hour_ordered = factor(hour_ordered))  # treat hour as categorical

lmm_fem <- lmer(
  relative_total_count ~ BPA_EXPOSURE * hour_ordered + (1 | ID),
  data = fem_data
)

summary(lmm_fem)
anova(lmm_fem)

emm <- emmeans(
  lmm_fem,
  ~ BPA_EXPOSURE | hour_ordered
)

## lights split ----

ical_phase <- ical_long_all %>%
  filter(!ID %in% c(9354, 9367, 9368, 9414,9406)) %>%  # 9367 and 9368 have confused data from 8/1/25 echoMRI 
  # 9354 lack of complete schedule of measurements in echoMRI
  # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
  # 9406 weird pattern in locomotion
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

dark_fem <- phase_summary %>%
  filter(
    SEX == "F",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = dark_fem
)

light_fem <- phase_summary %>%
  filter(
    SEX == "F",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = light_fem
)

dark_male <- phase_summary %>%
  filter(
    SEX == "M",
    lights_phase == "lights_off"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = dark_male
)

light_male <- phase_summary %>%
  filter(
    SEX == "M",
    lights_phase == "lights_on"
  )

t.test(
  total_relative_counts ~ BPA_EXPOSURE,
  data = light_male
)


deltaf <- phase_summary %>%
  filter(SEX == "F") %>%
  pivot_wider(
    names_from = lights_phase,
    values_from = total_relative_counts
  ) %>%
  mutate(delta_light_dark = lights_on - lights_off)

t.test(
  delta_light_dark ~ BPA_EXPOSURE,
  data = deltaf
)

deltam <- phase_summary %>%
  filter(SEX == "M") %>%
  pivot_wider(
    names_from = lights_phase,
    values_from = total_relative_counts
  ) %>%
  mutate(delta_light_dark = lights_on - lights_off)

t.test(
  delta_light_dark ~ BPA_EXPOSURE,
  data = deltam
)

delta_all <- phase_summary %>%
  pivot_wider(
    names_from = lights_phase,
    values_from = total_relative_counts
  ) %>%
  mutate(delta_light_dark = lights_on - lights_off)

delta_summary <- delta_all %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  summarise(
    mean_delta = mean(delta_light_dark, na.rm = TRUE),
    sem_delta  = sd(delta_light_dark, na.rm = TRUE) / sqrt(n()),
    n_ID = n(),
    .groups = "drop"
  )

ggplot(delta_summary, aes(x = BPA_EXPOSURE, y = mean_delta, fill = BPA_EXPOSURE)) +
  geom_col(alpha = 0.7, width = 0.6) +  # bars
  geom_errorbar(
    aes(ymin = mean_delta - sem_delta, ymax = mean_delta + sem_delta),
    width = 0.2
  ) +
  geom_jitter(
    data = delta_all,
    aes(x = BPA_EXPOSURE, y = delta_light_dark, color = BPA_EXPOSURE),
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
    y = "Delta (Lights ON – Lights OFF)",
    fill = "BPA Exposure",
    color = "BPA Exposure"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "top"
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

ggplot(
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

#CORRELATION ----

echoMRI_data_cor <- echoMRI_data_comparisons %>%
  filter(n_measurement %in% c("0 days with O.D", "133 days with O.D")) %>%
  select(ID, SEX, BPA_EXPOSURE, n_measurement, adiposity_index,COHORT) %>%
  pivot_wider(
    names_from  = n_measurement,
    values_from = adiposity_index
  ) %>%
  mutate(
    delta_ai = `133 days with O.D` - `0 days with O.D`
  ) %>% 
  drop_na()


delta_loco <- delta_all %>%
  select(ID, SEX, BPA_EXPOSURE, delta_light_dark) %>%
  drop_na(delta_light_dark)

delta_loco %>% count(SEX, BPA_EXPOSURE)

cor_df <- delta_loco %>%
  left_join(
    echoMRI_data_cor %>%
      select(ID, delta_ai),
    by = "ID"
  ) %>%
  drop_na(delta_ai)

cor_df %>%
  summarise(
    n_ID = n(),
    missing_ai = sum(is.na(delta_ai)),
    missing_loco = sum(is.na(delta_light_dark))
  )

cor_fem <- cor_df %>%
  filter(SEX == "F")

cor.test(
  cor_fem$delta_light_dark,
  cor_fem$delta_ai,
  method = "pearson"
)

# Compute correlation stats once
cor_res <- cor.test(
  cor_fem$delta_light_dark,
  cor_fem$delta_ai,
  method = "pearson"
)

label_text <- paste0(
  "Pearson r = ", round(cor_res$estimate, 2), "\n",
  "p = ", signif(cor_res$p.value, 2), "\n",
  "n = ", nrow(cor_fem)
)

ggplot(cor_fem, aes(x = delta_light_dark, y = delta_ai, color = BPA_EXPOSURE)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.1, vjust = 1.1,
    size = 4
  ) +
  scale_color_manual(values = c("NO" = "black", "YES" = "gray")) +
  labs(
    x = "Δ locomotion (Lights ON – Lights OFF)",
    y = "Δ adiposity index",
    color = "BPA exposure"
  ) +
  theme_minimal()

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

