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
  group_by(SEX,BPA_EXPOSURE) %>%
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

BW_gainsummary <- BW_gain  %>%
  group_by(day_rel,BPA_EXPOSURE,SEX) %>%
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
  facet_wrap( ~ SEX) +
  labs(
    x = "Days relative to first measurement",
    y = "Body weight gain",
    color = "BPA exposure",
    fill  = "BPA exposure"
  ) +
  theme_classic(base_size = 14)

# Speed = rate of change of body weight over time data----

BW_speed <- BW_data %>%
  group_by(ID, SEX, BPA_EXPOSURE) %>%
  do(tidy(lm(BW ~ day_rel, data = .))) %>%
  filter(term == "day_rel") %>%
  rename(
    speed_g_per_day = estimate,
    se = std.error
  )

ggplot(BW_speed,
       aes(x = BPA_EXPOSURE,
           y = speed_g_per_day,
           fill = BPA_EXPOSURE)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_grid(~ SEX) +
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

echoMRI_data %>% 
  group_by(SEX,BPA_EXPOSURE,n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) 

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

library(dplyr)
library(ggplot2)
library(forcats)

# Parse daytime_day to proper time
ical_long_all <- ical_long_all %>%
  mutate(datetime_parsed = as.POSIXct(paste("2000-01-01", datetime_day), format="%Y-%m-%d %H:%M"))

# Summarize counts per hour and group
ical_cumsum_summary <- ical_long_all %>%
  group_by(datetime_day, BPA_EXPOSURE, SEX) %>%
  summarise(count_hour = sum(count, na.rm = TRUE), .groups = "drop") %>%
  # Reorder daytime_day to start at 20:00 and go to 19:00
  mutate(
    hour_numeric = as.integer(format(as.POSIXct(paste("2000-01-01", datetime_day), format="%Y-%m-%d %H:%M"), "%H")),
    # Custom order: hours 20-23, then 0-19
    hour_order = case_when(
      hour_numeric >= 20 ~ hour_numeric,
      TRUE ~ hour_numeric + 24
    )
  ) %>%
  arrange(BPA_EXPOSURE, SEX, hour_order) %>%
  group_by(BPA_EXPOSURE, SEX) %>%
  mutate(count_cumsum = cumsum(count_hour)) %>%
  ungroup()

# Optional: factor for proper plotting
ical_cumsum_summary <- ical_cumsum_summary %>%
  mutate(datetime_day = fct_reorder(datetime_day, hour_order))

# Plot cumulative counts over “daytime_day” starting at 20:00
ggplot(ical_cumsum_summary, aes(x = datetime_day, y = count_cumsum, color = BPA_EXPOSURE, group = BPA_EXPOSURE)) +
  geom_line(size = 1) +
  facet_wrap(~SEX) +
  labs(x = "Time of Day", y = "Cumulative Counts", color = "BPA Exposure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



##heat####
ical_data_heat<- read_csv("~/Documents/GitHub/data/data/ical_results_Kotz_083024_heat.csv", skip = 2) 
colnames(ical_data_heat)[1:2] <- c("ID", "BW")

ical_long_heat <- ical_data_heat %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "kcal_hr") %>% 
  filter(!ID %in% c(7858, 7859)) %>%  #we want to eliminate the ID were used for antibody titration
  left_join(META, by = "ID")

ical_long_heat <- ical_long_heat %>%
  mutate(
    datetime_clean = str_remove(datetime_raw, "^x"),  # remove 'x' prefix
    datetime_clean = str_replace_all(datetime_clean, "_", "/"), # convert underscores to slashes
    datetime_parsed = lubridate::mdy_hm(datetime_clean)
  ) %>%
  select(ID, BW, kcal_hr, datetime_parsed,DIET_FORMULA) %>%
  mutate(
    date = as.Date(datetime_parsed),
    time = format(datetime_parsed, format = "%H:%M:%S"),
    hour = lubridate::hour(datetime_parsed),
    daycycle = ifelse(hour >= 20 | hour < 6, "dark", "light")
  ) %>%
  drop_na()

ical_long_heat <- ical_long_heat %>% 
  group_by(ID,daycycle) %>% 
  mutate(cumsumkcal_hr = cumsum(kcal_hr)) 
max_kcal_hr <- ical_long_heat %>%
  group_by(ID, daycycle,DIET_FORMULA) %>%
  summarise(max_kcal_hr = max(kcal_hr, na.rm = TRUE), .groups = "drop")

ggplot(max_kcal_hr, aes(x = daycycle, y = max_kcal_hr, group = ID)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA) +
  theme_minimal()

##RER####
ical_data_RER<- read_csv("~/Documents/GitHub/data/data/ical_results_Kotz_083024_RER.csv", skip = 1) 
colnames(ical_data_RER)[1:2] <- c("ID", "BW")

ical_long_RER <- ical_data_RER %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "RER") %>% 
  filter(!ID %in% c(7858, 7859)) %>%   #we want to eliminate the ID were used for antibody titration
  left_join(META, by = "ID")


ical_long_RER <- ical_long_RER %>%
  mutate(
    datetime_clean = str_remove(datetime_raw, "^x"),  # remove 'x' prefix
    datetime_clean = str_replace_all(datetime_clean, "_", "/"), # convert underscores to slashes
    datetime_parsed = lubridate::mdy_hm(datetime_clean)) %>%
  select(ID, BW, RER, datetime_parsed,DIET_FORMULA) %>%
  mutate(
    date = as.Date(datetime_parsed),
    time = format(datetime_parsed, format = "%H:%M:%S"),
    hour = lubridate::hour(datetime_parsed),
    daycycle = ifelse(hour >= 20 | hour < 6, "dark", "light")
  ) %>%
  drop_na()

ical_long_RER <- ical_long_RER %>% 
  group_by(ID,daycycle) %>% 
  mutate(meanRER = mean(RER)) 

maxRER <- ical_long_RER %>%
  group_by(ID, daycycle,DIET_FORMULA) %>%
  summarise(maxRER = max(RER, na.rm = TRUE), .groups = "drop")
ggplot(maxRER, aes(daycycle, maxRER,group = ID)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA)+
  theme_minimal()
# contextual object recognition test data ----

