#This script aim to explore adiposity index (fat/lean mass) in wt mice
#Our aim is to test no differences in adiposity index among 3 experimental group
#group 1 n=3: vehice (DMSO 2% and 98% corn oil) IP x 5 days and 2 days off
#group 2 n=3: RTIOXA-43 30 mg/kg IP (Yanan)
#group 3 n=3: RTIOXA-43 30 mg/kg IP (medchemexpress)

#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(car)
library(broom)
library(tidyr)
library(purrr)
library(ggpubr)

#body comp analysis####

##data import####
echomri_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import 

echomri_data <- echomri_data %>% 
  filter(COHORT ==9) %>% 
  select(ID, adiposity_index, n_measurement) %>% 
  mutate(GROUP = case_when(
    ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
    ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
    ID %in% c(8095, 8100, 8103) ~ "CONTROL")) 
 
 echomri_filtered <- echomri_data[echomri_data$n_measurement %in% c(2, 3), ]
  
## Summarize ####
summary_data <- echomri_filtered %>%
  group_by(n_measurement,GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

##plot before and after injection####

plot <- ggplot() +
  geom_point(data = echomri_filtered, aes(x = as.factor(n_measurement), y = adiposity_index), alpha = 0.7) +  # Optional: raw data points
  geom_point(data = summary_data, aes(x = as.factor(n_measurement), y = mean_adiposity), 
             color = "red", size = 3) +
  geom_line(data = echomri_filtered, aes(x = as.factor(n_measurement), y = adiposity_index, group = ID), alpha = 0.3) +
  geom_errorbar(data = summary_data, aes(x = as.factor(n_measurement), 
                                         ymin = mean_adiposity - se_adiposity, 
                                         ymax = mean_adiposity + se_adiposity),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("Adiposity Index") +
  xlab("Group") +
  scale_x_discrete(labels = c("2" = "Before", "3" = "After")) +
  facet_grid(~GROUP)
plot

##delta data####
delta_data <- echomri_filtered %>%
  filter(n_measurement %in% c(2, 3)) %>%
  pivot_wider(names_from = n_measurement, values_from = adiposity_index, names_prefix = "T") %>%
  mutate(delta_adiposity = T3 - T2)  # After - Before

## ANOVA test assumptions####
# Fit ANOVA model
anova_model <- aov(delta_adiposity ~ GROUP, data = delta_data)
# QQ plot
qqnorm(residuals(anova_model))
qqline(residuals(anova_model), col = "blue")

# Histogram
hist(residuals(anova_model), main = "Histogram of Residuals", xlab = "Residuals")

shapiro.test(residuals(anova_model))

plot(anova_model, which = 1)  # This is the default residuals vs fitted plot

leveneTest(delta_adiposity ~ GROUP, data = delta_data)

## One-way ANOVA####
anova_model <- aov(delta_adiposity ~ GROUP, data = delta_data)
summary(anova_model)


##plot of the deltas####

ggplot(delta_data, aes(x = GROUP, y = delta_adiposity, fill = GROUP)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  ylab("Δ Adiposity Index (After - Before)") +
  xlab("Group")


#food intake analysis####
##data import####
fi_data <-read_csv("~/Documents/GitHub/data/data/FI.csv") #data import 

fi_data <- fi_data %>% 
  filter(COHORT ==9) %>% 
  select(ID, DATE, corrected_intake_kcal, corrected_intake_gr) %>% 
  mutate(GROUP = case_when(
    ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
    ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
    ID %in% c(8095, 8100, 8103) ~ "CONTROL")) %>% 
  group_by(ID) %>% 
  mutate(cumulative_fi =cumsum(corrected_intake_kcal)) %>% 
  ungroup() %>% 
  drop_na()

##plot FI slopes longitudinal#### 
plot <- ggplot() +
  geom_point(data = fi_data, aes(x = DATE, y = cumulative_fi), alpha = 0.7) +  # Optional: raw data points
  geom_line(data = fi_data, aes(x = DATE, y = cumulative_fi, group = ID), alpha = 0.3) +
  theme_minimal() +
  geom_smooth(data = fi_data, aes(x = DATE, y = cumulative_fi, group = ID), 
              method = "lm", se = FALSE, color = "blue", linewidth = 0.6, linetype = "dashed") +  # Slopes
  theme_minimal() +
  ylab("cumulativeFood intake (kcal)") +
  xlab("DATE") +
  facet_grid(~GROUP)
plot

# make sure DATE is in numeric format (e.g., days since start)
fi_data <- fi_data %>%
  mutate(DATE_num = as.numeric(DATE - min(DATE)))

# Fit linear model for each animal
slopes <- fi_data %>%
  group_by(GROUP, ID) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(cumulative_fi ~ DATE_num, data = .x)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "DATE_num") %>%
  select(GROUP, ID, slope = estimate)

# View slopes
print(slopes)

##boxplot slopes####
ggplot(slopes, aes(x = GROUP, y = slope, fill = GROUP)) +
  geom_boxplot(alpha = 0.4) +
  geom_point(width = 0.2, size = 3, alpha = 0.7) +
  labs(
    y = "cumulative food intake over time (kcal/day)"
  ) +
  theme_minimal()
##one way ANOVA####
anova_result <- aov(slope ~ GROUP, data = slopes)
summary(anova_result) #the daily rate of food intake among the groups is not different

#bw analysis####
##data import####
bw_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import 

bw_data <- bw_data %>% 
  filter(COHORT ==9) %>% 
  select(ID, DATE, BW) %>% 
  mutate(GROUP = case_when(
    ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
    ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
    ID %in% c(8095, 8100, 8103) ~ "CONTROL")) %>% 
  filter(DATE %in% as.Date(c("2025-04-16", "2025-04-21"))) %>% 
  mutate(TIMEPOINT = case_when(
    DATE == as.Date("2025-04-16") ~ "Before",
    DATE == as.Date("2025-04-21") ~ "After"
  ))

## Summarize ####
summary_data_bw <- bw_data %>%
  group_by(TIMEPOINT,GROUP) %>%
  summarise(
    mean_bw = mean(BW, na.rm = TRUE),
    se_bw = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

##plot bw before and after injection####

plot <- ggplot() +
  geom_point(data = bw_data, aes(x = as.factor(TIMEPOINT), y = BW), alpha = 0.7) +  # Optional: raw data points
  geom_point(data = summary_data_bw, aes(x = as.factor(TIMEPOINT), y = mean_bw), 
             color = "red", size = 3) +
  geom_line(data = bw_data, aes(x = as.factor(TIMEPOINT), y = BW, group = ID), alpha = 0.3) +
  geom_errorbar(data = summary_data_bw, aes(x = as.factor(TIMEPOINT), 
                                         ymin = mean_bw - se_bw, 
                                         ymax = mean_bw + se_bw),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("BW (g)") +
  xlab("Group") +
  scale_x_discrete(labels = c("2" = "Before", "3" = "After")) +
  facet_grid(~GROUP)
plot

delta_data_bw <- bw_data %>%
  pivot_wider(names_from = TIMEPOINT, values_from = BW) %>%
  mutate(delta_bw = After - Before)

