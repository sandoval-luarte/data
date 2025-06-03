# This script aims to explore changes in body weight, food intake,
# body composition and locomotion in middle age NZO female mice after 
# different stages of feeding:1: Basal, 2: peak obesity, 3: Acute body 
#weight loss 4: BW maintenance
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()
library(ggpubr)

#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

#NZO####
##body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
         ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
         ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals


plot <- BW_data %>%
 # filter(DATE <= "2025-02-24") %>% # just raw BW over dates
  ggplot(aes(day_rel, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
   geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "day_rel",
    y = "BW (grams)")+
  format.plot

 plot
 
 # Subset rows to identify when sable measurements occurs 
 highlight_data <- BW_data %>%
   filter(COMMENTS == "RESTRICTED_DAY_1")
 
 plot <- BW_data %>% 
   ggplot(aes(DATE, bw_rel, group = ID)) +
   geom_point() +
   geom_line() +
   facet_wrap(~GROUP) +
   geom_smooth()+
   geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
   labs(
     x = "Date",
     y = "% BW gain")+
   # Highlight using vertical lines
   geom_vline(data = highlight_data, 
              aes(xintercept = as.numeric(DATE)), 
              linetype = "dashed", color = "red", alpha = 0.7)+
   format.plot
 
 plot
 
 plot <- BW_data %>%
   filter(DATE <= 2025-24-03) %>% #weight change over relative days 
   ggplot(aes(DATE, bw_rel, group = ID)) +
   geom_point() +
   geom_line() +
   facet_wrap(~GROUP) +
   geom_smooth()+
   geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
   labs(
     x = "Relative day",
     y = "% BW gain")+
   format.plot
 plot

#food intake kcal ####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 )) %>% #died during study
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
   group_by(ID) %>%
  arrange(DATE) %>% 
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         date_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) 

n_distinct(BW_data$ID) #here we know there is 22 animals

plot <- FI_data %>% #daily food intake
 filter(corrected_intake_kcal <= 40) %>% 
  ggplot(aes(date_rel,corrected_intake_kcal , group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_smooth()+
  labs(
    x = "Relative day",
    y = "daily kcal food intake")+
  format.plot
plot

#body comp analysis####
##echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  group_by(ID) %>%
  arrange(Date) %>% 
  filter(!ID %in% c(3712, 3715)) %>% #died during study 
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"))%>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index,GROUP) %>% 
  mutate(day_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean))

###adiposity index####
plot_echo <- echoMRI_data %>% #Adiposity index
  ggplot(aes(Date, adiposity_index, group = ID)) +
  geom_line() +
#  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)+
  format.plot
plot_echo

###fat mass####
plot_echo <- echoMRI_data %>% #Fat mass
  ggplot(aes(day_rel, Fat, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)+
  format.plot
plot_echo

### Lean mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, Lean, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)+
  format.plot
plot_echo #Interesting, ad lib guys still build muscle. It seem's like restricted ones just mantain muscle mass

### % change adiposity index####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, adiposity_index_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)+
  format.plot
plot_echo 

### % change fat mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, fat_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)+
  format.plot
plot_echo 

### % change lean mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, lean_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo 

#sable####
sable_data <- readRDS("../data/sable_downsampled_data.rds") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 ))   #died during study

#still exploring: locomotion####

sable_data_locomotion <- sable_data %>%  #Here I should add light and night per day
  filter(parameter =="AllMeters") %>%
  arrange(sable_idx)%>% 
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(hr=lubridate::hour(DateTime),
         day = lubridate::day(DateTime),
         day_night = if_else(hr %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),"light","dark")) %>% 
  filter(ID %in% c(3706, 3708)) %>% 
  ungroup() %>% 
  group_by(ID,day_night,day,GROUP) %>% 
  summarise(meters = mean(value))

plot_locomotion <- sable_data_locomotion %>% #I'm not sure what I am seeing
  ggplot(aes(day, meters)) +
  geom_point() + 
  facet_wrap(~GROUP*day_night)+ #LACKS STAGE HERE PEAK OBESITY AND BASAL AND MAINTENCE
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label

plot_locomotion

#C57BL/6J####
#body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT ==2) %>% #Just C57
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
           ID %in% c(7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
                     7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
           ID %in% c(7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals


plot <- BW_data %>% # just raw BW over dates
  ggplot(aes(day_rel, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP*DIET_FORMULA*SEX) +
  geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "Relative day",
    y = " % BW gain")

plot

# Subset rows to ide# Subset rows to ide# Subset rows to identify when sable measurements occurs 
highlight_data <- BW_data %>%
  filter(DATE == "2025-02-24")

plot <- BW_data %>% 
  ggplot(aes(DATE, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP*DIET_FORMULA) +
  geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "Date",
    y = "% BW gain")+
  # Highlight using vertical lines
  geom_vline(data = highlight_data, 
             aes(xintercept = as.numeric(DATE)), 
             linetype = "dashed", color = "red", alpha = 0.7)

plot

plot <- BW_data %>%  #weight change over relative days 
  ggplot(aes(day_rel, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_smooth()+
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + #ID label
  labs(
    x = "Relative day",
    y = "% BW gain")
plot

#food intake kcal ####
FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 )) %>% #died during study
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>%
  filter(!is.na(corrected_intake_gr)) %>% 
  mutate(corrected_intake_kcal = replace_na(corrected_intake_kcal, 0),) %>% 
  mutate(FI_rel = corrected_intake_kcal - first(corrected_intake_kcal),
         date_rel = DATE - first(DATE),
         FI_cum =cumsum(corrected_intake_kcal),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) 

n_distinct(BW_data$ID) #here we know there is 22 animals

plot <- FI_data %>% #daily food intake
  filter(corrected_intake_kcal <= 40) %>% 
  ggplot(aes(date_rel,corrected_intake_kcal , group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_smooth()+
  labs(
    x = "Relative day",
    y = "daily kcal food intake")
plot

#body comp analysis####
##echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT ==2) %>% #Just c57
  group_by(ID) %>% 
  arrange(Date) %>% 
  mutate(GROUP = case_when(
    ID %in% c(7876, 7873, 7875, 7864, 7862, 7867, 7869, 7870, 7871,
              7879, 7860, 7880, 7881, 7882, 7883, 7868) ~ "ad lib",
    ID %in% c(7872, 7874, 7878, 7865, 7861, 7863, 7877, 7866) ~ "restricted")) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index,GROUP) %>% 
  mutate(day_rel = Date - first(Date),
         adiposity_index_rel = 100 * (adiposity_index - first(adiposity_index)) / first(adiposity_index),
         fat_rel = 100 * (Fat - first(Fat)) / first(Fat),
         lean_rel = 100 * (Lean - first(Lean)) / first(Lean)) %>% 
  left_join(BW_data %>% select(ID, DIET_FORMULA,SEX), by = "ID") %>% 
  group_by(ID, day_rel) %>% 
  slice(1) 

###adiposity index####
plot_echo <- echoMRI_data %>%
  filter(GROUP == "restricted") 

plot_echo <- plot_echo %>%  
  ggplot(aes(day_rel, adiposity_index, group = ID)) +
  geom_line() +
  facet_wrap(~SEX) +
  geom_text(data = plot_echo,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo

###fat mass####
plot_echo <- echoMRI_data %>% #Fat mass
  ggplot(aes(day_rel, Fat, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP*DIET_FORMULA*SEX) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo

### Lean mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, Lean, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo #Interesting, ad lib guys still build muscle. It seem's like restricted ones just mantain muscle mass

### % change adiposity index####
plot_echo <- echoMRI_data %>% 
  ggplot(aes(day_rel, adiposity_index_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP*DIET_FORMULA*SEX) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo 

### % change fat mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, fat_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo 

### % change lean mass####
plot_echo <- echoMRI_data %>% #Lean mass
  ggplot(aes(day_rel, lean_rel, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo 

#sable####
sable_data <- readRDS("../data/sable_downsampled_data.rds") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715 ))   #died during study

#still exploring: locomotion####

sable_data_locomotion <- sable_data %>%  #Here I should add light and night per day
  filter(parameter =="AllMeters") %>%
  arrange(sable_idx)%>% 
  mutate(GROUP = case_when(
    ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
    ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(hr=lubridate::hour(DateTime),
         day = lubridate::day(DateTime),
         day_night = if_else(hr %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),"light","dark")) %>% 
  filter(ID %in% c(3706, 3708)) %>% 
  ungroup() %>% 
  group_by(ID,day_night,day,GROUP) %>% 
  summarise(meters = mean(value))

plot_locomotion <- sable_data_locomotion %>% #I'm not sure what I am seeing
  ggplot(aes(day, meters)) +
  geom_point() + 
  facet_wrap(~GROUP*day_night)+ #LACKS STAGE HERE PEAK OBESITY AND BASAL AND MAINTENCE
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6)  #ID label

plot_locomotion

