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
library(emmeans)
library(pracma)


 # BW over time data----

METABPA <- read_csv("~/Documents/GitHub/data/data/METABPA.csv")

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT==15) %>% # CD1 mice lack DOB
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


BW_summary <- BW_data %>%
  group_by(day_rel,BPA_EXPOSURE,DIET_FORMULA,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ggplot(BW_summary,
       aes(x = day_rel,
           y = mean_BW,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_BW - sem_BW,
                  ymax = mean_BW + sem_BW),
              alpha = 0.25,
              color = NA) +
  facet_wrap(DIET_FORMULA ~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

# OGTT----

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

## two way anova for AUC (SEX x BPA)----

auc_df$BPA_EXPOSURE <- factor(auc_df$BPA_EXPOSURE)
auc_df$SEX <- factor(auc_df$SEX)

model <- aov(AUC ~ BPA_EXPOSURE * SEX, data = auc_df)
summary(model)

#Males clearly have higher AUCs than females (F) overall
#There is no statistically supported difference between
#BPA-exposed and non-exposed females (p=0.52) non-significant BPA × SEX interaction

## three way anova for AUC (SEX x BPA x DIET_FORMULA)----
auc_df <- auc_df %>%
  mutate(
    BPA_EXPOSURE = factor(BPA_EXPOSURE),
    SEX = factor(SEX),
    DIET_FORMULA = factor(DIET_FORMULA)
  )

model3 <- aov(AUC ~ BPA_EXPOSURE * SEX * DIET_FORMULA, data = auc_df)
summary(model3)

#A significant main effect of sex was observed (F₁,₁₆ = 6.85, p = 0.019)
#whereas no main effect of BPA exposure was detected (p = 0.31).
#Diet showed a trend toward an effect on AUC (F₁,₁₆ = 3.95, p = 0.064). 
#No significant two- or three-way interactions among BPA exposure, 
#sex, and diet were observed


# Body comp over time ----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16)) %>% 
  group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index,COHORT) %>% 
  left_join(METABPA, by= "ID") 

bodycomp_summary <- echoMRI_data %>%
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
  )

## adiposity index ----
ggplot(bodycomp_summary,
       aes(x = n_measurement,
           y = mean_AI,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_AI - sem_AI,
                  ymax = mean_AI + sem_AI),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "adiposity index (fat/lean mass)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

#here two weird things: 
#first males has higher basal adiposity index than females
#females exposed to BPA started with higher adiposity index 
#than females non exposed to BPA and the trend is 
#to normalizing adiposity index overtime 

#now the question is, after how many weeks of HFD or LFD 
#the first measured were done?

#is these effects dependent of DIET_FORMULA? YES, IS HIGHER IN FEMALES
#what if we collapsed by SEX and split by DIET_FORMULA? THE EFFECT IS HIGHER IN HFD
#is these effects dependent of COHORT? YES, IT SEEMS FOR COHORT 16 THIS IS NOT TRUE BUT I BELIEVE IS BECAUSE LACK OF STAT POWER

## fat mass ----

ggplot(bodycomp_summary,
       aes(x = n_measurement,
           y = mean_fat,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_fat - sem_fat,
                  ymax = mean_fat + sem_fat),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "fat mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

#here the same trend than for adiposity index: 
#first males has higher basal fat than females
#females exposed to BPA started with higher fat mass 
#than females non exposed to BPA and the trend is 
#to normalizing fat mass overtime 

## lean mass ----

ggplot(bodycomp_summary,
       aes(x = n_measurement,
           y = mean_lean,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_lean - sem_lean,
                  ymax = mean_lean + sem_lean),
              alpha = 0.25,
              color = NA) +
  facet_wrap(~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "lean mass (g)",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

#wow this is crazy! 
#females exposed to BPA started with higher lean mass 
#than females non exposed to BPA and the trend is not trending 
#to normalizing overtime 

