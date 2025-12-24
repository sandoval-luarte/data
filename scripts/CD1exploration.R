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
library(lubridate)
library(broom)
library(stringr)
library(forcats)


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
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 


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

# BW gain over time data----

BW_gain <- BW_data %>% 
  select(ID,day_rel,bw_rel,BPA_EXPOSURE,DIET_FORMULA,SEX)

BW_gain %>% 
  group_by(SEX,BPA_EXPOSURE,DIET_FORMULA) %>%
  summarise(n_ID = n_distinct(ID)) 

BW_gainsummary <- BW_gain  %>%
  group_by(day_rel,BPA_EXPOSURE,SEX,DIET_FORMULA) %>%
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
  facet_wrap(DIET_FORMULA ~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight gain",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

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

ggplot(BW_speed,
       aes(x = BPA_EXPOSURE,
           y = speed_g_per_day,
           fill = BPA_EXPOSURE)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_grid(DIET_FORMULA~ SEX) +
  labs(
    y = "BW gain speed (g of BW/day)",
    x = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

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
# BUT t = 1.36 so this is not statistically significant.

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
 facet_wrap(DIET_FORMULA~ SEX) +
  labs(
    x = "BPA exposure",
    y = "Glucose AUC (0–90 min)"
  ) +
  theme_classic()

## two way anova for AUC (SEX x BPA)----

auc_df$BPA_EXPOSURE <- factor(auc_df$BPA_EXPOSURE)
auc_df$SEX <- factor(auc_df$SEX)

model <- aov(AUC ~ BPA_EXPOSURE * SEX *DIET_FORMULA, data = auc_df)
summary(model)

#Males clearly have higher AUCs than females (F) overall
#There is no statistically supported difference between
#BPA-exposed and non-exposed females (p=0.52) non-significant BPA × SEX interaction
#there is a trend (p=0.06) for HFD to increased AUC

# Body comp over time ----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(15,16)) %>% 
  filter(!ID %in% c(9354, 9367, 9368, 9372, 9414)) %>%  # 9367 and 9368 have confused data from 8/1/25 echoMRI 
                                                        # 9354 and 9372 lack of complete schedule of measurements in echoMRI
                                                        # 9414 lack of complete schedule of measurements (lack of basal) in echoMRI
  group_by(ID) %>%
  arrange(Date) %>% 
  select(ID, Date, Fat, Lean, Weight, adiposity_index,COHORT) %>% 
  left_join(METABPA, by= "ID") %>% 
  ungroup()

#cohort 15 are 24 animals originally but if I eliminate 4 animals so total animals are 20 per date (9367, 9368, 9354, 9372)
#cohort 16 are 15 animals originally but if I eliminate 1 animal so total animals are 14 per date (9414)

echoMRI_data %>% 
  group_by(COHORT,SEX,BPA_EXPOSURE) %>%
  summarise(n_ID = n_distinct(ID)) 


echoMRI_data_comparisons <- echoMRI_data %>% 
  mutate(
    n_measurement= case_when(
      COHORT == 15 & Date == "2025-08-01" ~ 1,
      COHORT == 16 & Date == "2025-11-17" ~ 1,
      COHORT == 15 & Date == "2025-09-09" ~ 2,
      COHORT == 16 & Date == "2025-12-18" ~ 2)) %>% 
  drop_na()
                                               

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
  theme_classic(base_size = 14)+
  scale_x_continuous(breaks = c(1, 2))

#here two weird things: 
#first males has higher basal adiposity index than females (betweeb 16 to 18 week old)
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
  theme_classic(base_size = 14)+
  scale_x_continuous(breaks = c(1, 2))

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
  theme_classic(base_size = 14)+
  scale_x_continuous(breaks = c(1, 2))

#wow this is crazy! 
#females exposed to BPA started with higher lean mass 
#than females non exposed to BPA and the trend is not trending 
#to normalizing overtime 

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
    -DIET_FORMULA.x
  ) %>% 
  rename(
    SEX = SEX.x
  ) %>% 
  mutate(
   FIcumulative = cumsum(corrected_intake_kcal))


FI_summary <- FI_data %>%
  group_by(day_rel,BPA_EXPOSURE,SEX) %>%
  summarise(
    mean_FI = mean(FIcumulative, na.rm = TRUE),
    sem_FI  = sd(FIcumulative, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ggplot(FI_summary,
       aes(x = day_rel,
           y = mean_FI,
           color = BPA_EXPOSURE,
           fill  = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_FI - sem_FI,
                  ymax = mean_FI + sem_FI),
              alpha = 0.25,
              color = NA) +
  facet_wrap( ~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "food intake (kcal) over time",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

#independent of the diet females ate less than males, curious

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

ical_long15 <- ical_long15 %>% 
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
  group_by(datetime_day) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great

ical_long_all  <- ical_long_all %>% 
  select(ID, BW, count, BPA_EXPOSURE, SEX, DIET_FORMULA, day, datetime_day, cohort) %>% 
  mutate(datetime_day = str_remove(datetime_day, "^day \\d+\\s+"))

ical_long_all <- ical_long_all %>%
  mutate(
    datetime_parsed = as.POSIXct(paste("2000-01-01", datetime_day), format="%Y-%m-%d %H:%M"),
    hour = as.integer(format(datetime_parsed, "%H")),
    minute = as.integer(format(datetime_parsed, "%M")),
    daytime = case_when(
      (day == 1 & hour >= 20) | (day == 2 & hour < 6) ~ "night",
      (day == 2 & hour >= 7) ~ "light",
      TRUE ~ "light"
    )
  ) %>%
  select(-datetime_parsed, -hour, -minute) %>% 
  ungroup()

# Parse daytime_day to proper time
ical_long_all <- ical_long_all %>%
  mutate(datetime_parsed = as.POSIXct(paste("2000-01-01", datetime_day), format="%Y-%m-%d %H:%M"))


#  Calculate mean and SEM per hour per group
ical_hourly_summary <- ical_long_all %>%
  group_by(SEX, BPA_EXPOSURE, datetime_day,daytime) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    sem_count = sd(count, na.rm = TRUE) / sqrt(n_distinct(ID)),
    .groups = "drop"
  ) %>%
  # parse hour for proper ordering (20:00 → 19:00 next day)
  mutate(
    hour_numeric = as.integer(format(as.POSIXct(paste("2000-01-01", datetime_day), format="%Y-%m-%d %H:%M"), "%H")),
    hour_order = ifelse(hour_numeric >= 20, hour_numeric, hour_numeric + 24)
  ) %>%
  arrange(SEX, BPA_EXPOSURE, hour_order) %>%
  mutate(datetime_day = forcats::fct_reorder(datetime_day, hour_order))

#  Compute cumulative sum of mean counts over hours
ical_hourly_summary2 <- ical_hourly_summary %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  mutate(cumsum_mean = cumsum(mean_count),
         cumsum_sem = sqrt(cumsum(sem_count^2))) %>%  # propagate SEM for cumulative sum
  ungroup()

#  Plot cumulative mean ± SEM
ggplot(ical_hourly_summary2, aes(x = datetime_day, y = cumsum_mean, color = BPA_EXPOSURE, group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = cumsum_mean - cumsum_sem, ymax = cumsum_mean + cumsum_sem, fill = BPA_EXPOSURE),
              alpha = 0.2, color = NA) +
  facet_wrap(~SEX) +
  labs(x = "Time of Day", y = "Cumulative Mean Count ± SEM", color = "BPA Exposure", fill = "BPA Exposure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ical_hourly_relative <- ical_hourly_summary2 %>%
  group_by(SEX, BPA_EXPOSURE) %>%
  # subtract first hour's cumulative mean to start at 0
  mutate(
    rel_cumsum_mean = cumsum_mean - first(cumsum_mean),
    rel_cumsum_sem  = cumsum_sem  # SEM remains the same, optionally can recalc relative
  ) %>%
  ungroup()

# Plot relative cumulative mean ± SEM
ggplot(ical_hourly_relative, aes(x = datetime_day, y = rel_cumsum_mean, color = BPA_EXPOSURE, group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = rel_cumsum_mean - rel_cumsum_sem, ymax = rel_cumsum_mean + rel_cumsum_sem, fill = BPA_EXPOSURE),
              alpha = 0.2, color = NA) +
  facet_wrap(~SEX) +
  labs(x = "Time of Day", y = "Relative Cumulative Mean Count ± SEM", color = "BPA Exposure", fill = "BPA Exposure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Assuming you have your hourly summary:
# datetime_day, BPA_EXPOSURE, SEX, mean_count, sem_count, daytime

ical_daytime_relcumsum <- ical_hourly_summary %>%
  group_by(SEX, BPA_EXPOSURE, daytime) %>%
  arrange(datetime_day) %>%
  mutate(
    # Cumulative sum relative to the first hour of the daytime period
    rel_cumsum_mean = cumsum(mean_count) - first(mean_count),
    rel_cumsum_sem  = sqrt(cumsum(sem_count^2))  # SEM of cumulative sum
  ) %>%
  ungroup()

ggplot(ical_daytime_relcumsum, aes(x = datetime_day, y = rel_cumsum_mean,
                                   color = BPA_EXPOSURE, group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = rel_cumsum_mean - rel_cumsum_sem,
                  ymax = rel_cumsum_mean + rel_cumsum_sem,
                  fill = BPA_EXPOSURE), alpha = 0.2, color = NA) +
  facet_wrap(~SEX + daytime, scales = "free_x") +
  labs(x = "Hour", y = "Relative Cumulative Count ± SEM",
       color = "BPA Exposure", fill = "BPA Exposure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ical_auc_per_ID <- ical_long_all %>%
  arrange(ID, daytime, datetime_day) %>%
  group_by(ID, SEX, BPA_EXPOSURE, daytime) %>%
  summarise(
    AUC = trapz(
      x = as.numeric(factor(datetime_day, levels = unique(datetime_day))),
      y = count
    ),
    .groups = "drop"
  )

table(ical_auc_per_ID$SEX,
      ical_auc_per_ID$BPA_EXPOSURE,
      ical_auc_per_ID$daytime)

anova_night_auc <- aov(
  AUC ~ SEX * BPA_EXPOSURE,
  data = filter(ical_auc_per_ID, daytime == "night")
)

summary(anova_night_auc)

anova_light_auc <- aov(
  AUC ~ SEX * BPA_EXPOSURE,
  data = filter(ical_auc_per_ID, daytime == "light")
)

summary(anova_light_auc)

# Summarize max relative cumulative counts per daytime
ical_daytime_max <- ical_daytime_relcumsum %>%
  group_by(SEX, BPA_EXPOSURE, daytime) %>%
  summarise(
    max_rel_cumsum = max(rel_cumsum_mean, na.rm = TRUE),
    sem_at_max = sem_count[which.max(rel_cumsum_mean)], # optional SEM at max
    .groups = "drop"
  )

# Column plot
ggplot(ical_daytime_max, aes(x = daytime, y = max_rel_cumsum,
                             fill = BPA_EXPOSURE)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = max_rel_cumsum - sem_at_max,
                    ymax = max_rel_cumsum + sem_at_max),
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_wrap(~SEX) +
  labs(
    x = "Daytime",
    y = "Max Relative Cumulative Count",
    fill = "BPA Exposure"
  ) +
  theme_minimal()


# 1. Compute per-animal cumulative relative count per daytime
ical_daytime_per_ID <- ical_long_all %>%
  group_by(ID, SEX, BPA_EXPOSURE, daytime) %>%
  arrange(datetime_day) %>%
  mutate(
    rel_cumsum = cumsum(count) - first(count)
  ) %>%
  summarise(
    max_rel_cumsum = max(rel_cumsum, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Check the data
table(ical_daytime_per_ID$SEX, ical_daytime_per_ID$BPA_EXPOSURE)
# Now you should have multiple IDs per group

# 3. Run ANOVA
anova_night <- aov(max_rel_cumsum ~ SEX * BPA_EXPOSURE, 
                   data = filter(ical_daytime_per_ID, daytime == "night"))
summary(anova_night)

anova_light <- aov(max_rel_cumsum ~ SEX * BPA_EXPOSURE, 
                   data = filter(ical_daytime_per_ID, daytime == "light"))
summary(anova_light)



