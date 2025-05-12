#Description####
#This script aims to explore the relationship between orexin neuron activation and diet.
#Our hypothesis is that HFD mice will have less orexin neuron activation within LHA and locomotion 

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(janitor) #to clean columns


#metadata import####
META<- read_csv("~/Documents/GitHub/data/data/META.csv") %>% 
  filter(COHORT==0) %>% #We just want clarity guys
  filter(!ID %in% c(7858, 7859)) %>% #we want to eliminate the ID were used for antibody titration
  select(ID,DIET_FORMULA)
  
#bw data import####
BW_data <- read_csv("~/Documents/GitHub/data/data/BW.csv") %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = BW - first(BW),
    date_rel = DATE - first(DATE)
  ) %>% 
  filter(COHORT==0) %>% #We just want clarity guys
  filter(!ID %in% c(7858, 7859)) #we want to eliminate the ID were used for antibody titration

ggplot(BW_data, aes(date_rel, bw_rel)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA) +
  geom_smooth() +
  theme_minimal()

#fi data import####
FI_data <- read_csv("~/Documents/GitHub/data/data/FI.csv") %>% 
  filter(COHORT==0) %>% #We just want clarity guys
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  group_by(ID) %>%
 mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0)) %>% #we jump the days in which there is no data such as the day of the echoMRI
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
    date_rel = DATE - first(DATE),
    FI_cum =cumsum(corrected_intake_kcal)) %>% 
  filter(!ID %in% c(7858, 7859))  #we want to eliminate the ID were used for antibody titration

ggplot(FI_data, aes(date_rel, FI_cum)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  facet_wrap(~DIET_FORMULA) +
  geom_smooth() +
  theme_minimal()

#ical data import####
##counts####
ical_data_counts <- read_csv("~/Documents/GitHub/data/data/ical_results_Kotz_083024_counts.csv", skip = 1) 
colnames(ical_data_counts)[1:2] <- c("ID", "BW")

ical_long <- ical_data_counts %>%
  pivot_longer(cols = -c(ID, BW), 
               names_to = "datetime_raw", 
               values_to = "count") %>% 
  filter(!ID %in% c(7858, 7859))  %>% #we want to eliminate the ID were used for antibody titration
  left_join(META, by = "ID")

ical_long <- ical_long %>%
  mutate(
    datetime_clean = str_remove(datetime_raw, "^x"),  # remove 'x' prefix
    datetime_clean = str_replace_all(datetime_clean, "_", "/"), # convert underscores to slashes
    datetime_parsed = lubridate::mdy_hm(datetime_clean)
  ) %>%
  select(ID, BW, count, datetime_parsed,DIET_FORMULA) %>%
  mutate(
    date = as.Date(datetime_parsed),
    time = format(datetime_parsed, format = "%H:%M:%S"),
    hour = lubridate::hour(datetime_parsed),
    daycycle = ifelse(hour >= 20 | hour < 6, "dark", "light")
  ) %>%
  drop_na()

ical_counts <- ical_long %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(period= cumsum(as.numeric(replace_na(daycycle != lag(daycycle),0)))) %>% 
  ungroup() %>% 
  group_by(ID,daycycle) %>% 
  mutate(period_rel = scales::rescale(period,to= c(0,1))) %>% 
  ungroup() %>% 
  group_by(ID,daycycle,period_rel,DIET_FORMULA) %>% 
summarise(total_movement =sum(count))

counts_model <- ical_counts %>% 
  ungroup() %>% 
  group_by(ID,daycycle, DIET_FORMULA) %>% 
  summarise(mean_counts = mean(total_movement),
            sd_counts = sd(total_movement))


ggplot(counts_model, aes(daycycle, mean_counts, group = ID)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  geom_line() +
  facet_wrap(~DIET_FORMULA) +
  theme_minimal()


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

#echoMRI data import####
echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>% 
  filter(COHORT==0) %>% 
  filter(!ID %in% c(7858, 7859), n_measurement ==2) %>% #I only consider the second measurement because measurement 1 and 2 has just 7 weeks apart
  group_by(ID)
# Now run the unpaired t-test
t_test_result <- t.test(adiposity_index ~ DIET_FORMULA, data = echoMRI_data,var.equal = TRUE)
# View the result
print(t_test_result)

ggplot(echoMRI_data, aes(DIET_FORMULA,adiposity_index)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
  theme_minimal()

#clarity data import####
cell_count <- read_csv("~/Documents/GitHub/data/data/Kotz_Cellcounts.csv") %>% 
  rename(ID = `Sample Name`, cellcounts =`Cells Detected Total (Curated)`) %>%
  mutate(ID = str_remove(ID, "Kotz_"), ID= as.numeric(ID)) %>% 
  select(ID, cellcounts) %>% 
  left_join(META, by = "ID") %>% 
  #filter(!ID == 4021) %>%  
  drop_na() %>% 
  left_join(counts_model, by = c("ID","DIET_FORMULA")) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(rel_movement = (mean_counts[daycycle=="dark"] - mean_counts[daycycle=="light"])/
           mean_counts[daycycle=="light"])


# Now run the unpaired t-test
t_test_result <- t.test(cellcounts ~ DIET_FORMULA, data = cell_count,var.equal = TRUE)
# View the result
print(t_test_result)

ggplot(cell_count, aes(DIET_FORMULA,cellcounts,group = ID)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = ID), vjust = -0.5, size = 3) +  # add ID labels above points
theme_minimal()

neuron_locomotion_model <- lm(data=cell_count %>% 
                                filter(daycycle=="dark"), rel_movement ~ cellcounts*DIET_FORMULA)
summary(neuron_locomotion_model)
