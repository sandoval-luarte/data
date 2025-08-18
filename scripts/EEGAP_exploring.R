# This script aims to explore changes in body weight, food intake and body comp
# in middle age NZO and C57 after different stages of feeding:
#1: Basal, 2: peak obesity, 3: Acute body weight loss, 4: BW maintenance, 5: chronic RTIOXA-47 injections
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()
library(ggpubr)


##body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
         ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726,7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
                   7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
         ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 46 animals,
                        # n=22 NZO, n=24 C57


highlight_data <- BW_data %>%
  filter(COMMENTS %in% c("DAY_1_INJECTIONS")) #indirect calorimetry measurement


plot <- BW_data %>%
  filter(STRAIN =="NZO/HlLtJ") %>% 
  ggplot(aes(DATE, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
 #  geom_smooth()+
 # geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "day_rel",
    y = "BW (grams)")+
  geom_vline(data = highlight_data, 
             aes(xintercept = as.numeric(DATE)), 
             linetype = "dashed", color = "red", alpha = 0.7)+
  format.plot

plot

#body comp analysis####
##echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT > 1 & COHORT < 6) %>% #Just NZO females
  arrange(Date) %>% 
  filter(!ID %in% c(3712, 3715)) %>% #died during study 
 mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726,7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
              7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted"))%>% 
  group_by(ID) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index,GROUP) %>% 
  mutate(day_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean)) %>% 
  left_join(BW_data %>% select(ID, DIET_FORMULA,SEX,STRAIN), by = "ID") 

NZO_injections <-echoMRI_data %>%
  filter(STRAIN =="NZO/HlLtJ") %>% 
  slice_tail(n = 1) %>% 
  mutate(INJECTION = case_when(
    ID %in% c(3714,
              3728,
              3725,
              3724,
              3727,
              3720,
              3707,
              3709,
              3711,
              3713,
              3706) ~ "Vehicle",
    ID %in% c(3710,
              3729,
              3708,
              3723,
              3721,
              3722,
              3726,
              3719,
              3718,
              3716,
              3717) ~ "rtioxa_47"))
  


C57_injections <-echoMRI_data %>%
  filter(STRAIN =="C57BL6/J") %>% 
  slice_tail(n = 1) %>% 
  mutate(INJECTION = case_when(
    ID %in% c(7880,
              7881,
              7882,
              7883,
              7869,
              7870,
              7871,
              7868,
              7872,
              7878,
              7861,
              7863,
              7875,
              7876,
              7867,
              7864)  ~ "Vehicle",
    ID %in% c(7874,
              7877,
              7865,
              7866,
              7873,
              7879,
              7860,
              7862) ~ "rtioxa_47"))

###adiposity index####
plot_echo <- NZO_injections %>%
  ggplot(aes(INJECTION, adiposity_index, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
   geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label
  
plot_echo

###injection assignation NZO####
plot_inject <- echoMRI_data %>%
  ggplot(aes(day_rel, adiposity_index, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP*SEX) +
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label

plot_inject

###adiposity index####
plot_echo <- C57_injections %>%
  filter(DIET_FORMULA=="D12451i") %>% 
  ggplot(aes(INJECTION, adiposity_index, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP*SEX) +
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label

plot_echo

library(dplyr)
library(purrr)
library(broom)

# Filter relevant data
anova_data <- C57_injections %>%
  filter(DIET_FORMULA == "D12451i", !is.na(INJECTION))

# Split data by GROUP and SEX
split_data <- anova_data %>%
  group_by(GROUP, SEX) %>%
  group_split()

# Create informative names for each split
group_labels <- anova_data %>%
  distinct(GROUP, SEX) %>%
  mutate(label = paste(GROUP, SEX, sep = "_")) %>%
  pull(label)

# Run ANOVA on each subgroup
anova_results <- map(split_data, ~ {
  if (n_distinct(.x$INJECTION) > 1) {
    # only run ANOVA if both injection groups are present
    model <- aov(adiposity_index ~ INJECTION, data = .x)
    broom::tidy(model)
  } else {
    tibble(term = NA, p.value = NA, warning = "Only one injection group present")
  }
})

# Set names
names(anova_results) <- group_labels

# Print results
anova_results

##SLOPE AFTER DAY 1 OF INJECTIONS FOR NZO MICE
# Load necessary libraries
library(dplyr)
library(readr)
library(purrr)
library(broom)

# Read and clean the data
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>%  # Only NZO females
  filter(!ID %in% c(3712, 3715)) %>%   # Exclude animals that died
  arrange(ID, DATE) %>% 
  group_by(ID) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"),
         DRUG = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
           ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"),
         day_rel = DATE - first(DATE)) %>% 
  ungroup()

# Step 1: Get injection day per ID
injection_days <- BW_data %>%
  filter(COMMENTS == "DAY_1_INJECTIONS") %>%
  select(ID, inj_date = DATE)

# Step 2: Join to original data and filter post-injection
BW_post_injection <- BW_data %>%
  left_join(injection_days, by = "ID") %>%
  filter(DATE >= inj_date)


plot <- BW_post_injection %>%
  ggplot(aes(DATE, BW, group = ID)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, aes(group = ID), color = "blue")+ 
  facet_wrap(~GROUP*DRUG) +
  labs(
    x = "Date",
    y = "BW (grams)")
plot

# Step 3: Compute slope (BW ~ DATE) per ID
slopes <- BW_post_injection %>%
  group_by(ID) %>%
  filter(n() >= 2) %>%  # Make sure there's enough data
  do(tidy(lm(BW ~ as.numeric(DATE), data = .))) %>%
  filter(term == "as.numeric(DATE)") %>%
  select(ID, slope = estimate)

# Preview the result
print(slopes)

# 1. Merge slope with group info (just once per ID)
slopes_grouped <- slopes %>%
  left_join(BW_data %>% select(ID, DRUG,GROUP) %>% distinct(), by = "ID")

# 2. Plot: Boxplot + individual points
library(ggplot2)

ggplot(slopes_grouped, aes(x = GROUP, y = slope)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Boxplot without outlier dots
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) + # Add individual slopes
  labs(title = "Slope of BW Change After DAY_1_INJECTIONS",
       y = "Slope (g/day)",
       x = "Group") +
  theme_minimal() +
  theme(legend.position = "none")+
  facet_grid(GROUP~DRUG) +
  format.plot

# Subset the data
adlib_data <- slopes_grouped %>% filter(GROUP == "ad lib")
restricted_data <- slopes_grouped %>% filter(GROUP == "restricted")

# Compare vehicle vs RTIOXA_47 within each diet group using Wilcoxon rank-sum test
adlib_test <- wilcox.test(slope ~ DRUG, data = adlib_data)
restricted_test <- wilcox.test(slope ~ DRUG, data = restricted_data)

# Show results
adlib_test
restricted_test

# Compare means between vehicle and RTIOXA_47 within the ad lib group
t.test(slope ~ DRUG, data = adlib_data)

# Compare means between vehicle and RTIOXA_47 within the restricted group
t.test(slope ~ DRUG, data = restricted_data)



