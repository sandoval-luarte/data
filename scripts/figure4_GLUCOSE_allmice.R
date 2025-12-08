#we aim to evaluate changes in blood glucose levels in NZO and C57 mice during baseline, peak obesity
#BW loss, BW maintenance and BW regain phases. During BW regain RTIOXA 47 or vehicle was injected
#for 5 weeks

#libraries

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(car)
library(patchwork)
library(ggrepel)
library(tidyverse)


#format plot
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

# general data import ----

# List of CSV files
files <- c("COHORT_2.csv","COHORT_3.csv", "COHORT_4.csv", "COHORT_5.csv")

# Function to read and fix SEX
fix_sex <- function(file){
  read_csv(file.path("~/Documents/GitHub/data/data", file),
           col_types = cols(.default = "c")) %>%  # read all as character to avoid FALSE/TRUE issue
    mutate(SEX = ifelse(SEX %in% c("FALSE", FALSE), "F",
                        ifelse(SEX %in% c("TRUE", TRUE), "M", SEX)))
}

# Apply to all files and combine if needed
all_cohorts_fixed <- map_df(files, fix_sex)

# Optional: save fixed files
walk(files, ~ {
  df_fixed <- fix_sex(.x)
  write_csv(df_fixed, paste0("~/Documents/GitHub/data/data/", tools::file_path_sans_ext(.x), "_fixed.csv"))
})

# Check the result
unique(all_cohorts_fixed$SEX)


gludata <- all_cohorts_fixed %>% 
  mutate(DATE = lubridate::mdy(DATE)) %>% 
  mutate(FASTED_GLU_mg_dL= as.numeric(FASTED_GLU_mg_dL)) %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  filter(!DIET == "CHOW") %>% 
  select(ID,STRAIN,SEX,BODY_WEIGHT_G, DATE, COMMENTS, DIET_FORMULA,FASTED_GLU_mg_dL)  %>% 
  drop_na(FASTED_GLU_mg_dL) 

gludata <-gludata %>% 
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 
                7876, 7879, 7880, 7881,7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,
                7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    )) %>% 
  filter(!ID %in% c(3712, 3715)) %>% # died during study 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(
    n_measurement = row_number(),      # 1 for first DATE, 2 for second, etc.
    delta_measurement = DATE - lag(DATE)
  ) %>% 
  mutate(
    STATUS = case_when(
      n_measurement == 1 ~ "peak obesity",
      n_measurement == 2 ~ "BW loss",
      TRUE ~ "other"  # everything else
    )
  ) %>% 
  mutate(STATUS= factor(STATUS, levels = c("peak obesity", "BW loss"))) %>% 
  filter(STATUS != "other")


gludata <- gludata %>%
  mutate(
    STATUS = case_when(
      STATUS == "BW loss" & GROUP == "ad lib" ~ "BW maintenance",
      TRUE ~ as.character(STATUS)  # keep everything else the same
    ),
    STATUS = factor(STATUS, levels = c("peak obesity", "BW maintenance", "BW loss"))  # reorder factor if needed
  )

#NZO analysis ----

gludata_NZO <-gludata %>% 
filter(STRAIN == "NZO/HlLtJ")

gludata_NZO %>% 
  group_by(n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) 

# Find IDs missing n_measurement == 2
missing_measurement2 <- gludata_NZO %>%
  group_by(ID) %>%
  summarise(has_measure2 = any(n_measurement == 2)) %>%
  filter(!has_measure2) %>%
  pull(ID)

missing_measurement2 #check why we measured just once blood glucose in ID 3727. Did she was sick or somethinhg?

gludata_NZO <- gludata %>%
  filter(STRAIN == "NZO/HlLtJ") %>%
  filter(ID != 3727) %>%  # remove animals with only one measurement
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    n_measurement = row_number(),
    STATUS = case_when(
      n_measurement == 1 ~ "peak obesity",
      n_measurement == 2 & GROUP == "ad lib" ~ "BW maintenance",  # correctly rename for ad lib
      n_measurement == 2 & GROUP == "restricted" ~ "BW loss",
      TRUE ~ "other"
    ),
    STATUS = factor(STATUS, levels = c("peak obesity", "BW maintenance", "BW loss"))
  ) %>%
  ungroup() %>%
  mutate(
    highlight = case_when(
      # Restricted group: glucose increased
      GROUP == "restricted" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) > 0 ~ "Restricted Increase",
      # Ad lib group: glucose decreased
      GROUP == "ad lib" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) < 0 ~ "Ad lib Decrease",
      TRUE ~ "Other"
    )
  )


gludata_NZO <- gludata_NZO %>%
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    highlight = case_when(
      # Restricted group: glucose increased
      GROUP == "restricted" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) > 0 ~ "Restricted Increase",
      # Ad lib group: glucose decreased
      GROUP == "ad lib" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) < 0 ~ "Ad lib Decrease",
      TRUE ~ "Other"
    )
  ) %>%
  ungroup()

# Remove unused STATUS levels per GROUP
gludata_NZO_plot <- gludata_NZO %>%
  group_by(GROUP) %>%
  filter(STATUS %in% unique(STATUS)) %>%   # keep only STATUS that exist in this group
  ungroup() %>%
  mutate(STATUS = droplevels(STATUS))      # drop unused factor levels

# Align points in the center
gludata_NZO_plot_centered <- gludata_NZO_plot %>%
  mutate(x_pos = as.numeric(STATUS))  # convert factor to numeric for exact positions

# Then plot 
plot <- ggplot(gludata_NZO_plot, aes(x = STATUS, y = FASTED_GLU_mg_dL, fill = STATUS)) + 
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) + 
  geom_point(aes(color = highlight), alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) + 
  geom_text(aes(label = ID), size = 3, vjust = -0.5)+
scale_color_manual(values = c( "Restricted Increase" = "red", "Ad lib Decrease" = "blue", "Other" = "black" )) + 
  theme_minimal() + labs(x = NULL, y = "Fasted glucose (mg/dL)") +
  facet_wrap(~GROUP, scales = "free_x") + 
  theme(legend.position = "none") 
plot
  

#C57 analysis ----
  
gludata_C57 <-gludata %>% 
  filter(STRAIN == "C57BL/6J")

gludata_C57 %>% 
  group_by(n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) 

# Find IDs missing n_measurement == 2
missing_measurement2 <- gludata_C57 %>%
  group_by(ID) %>%
  summarise(has_measure2 = any(n_measurement == 2)) %>%
  filter(!has_measure2) %>%
  pull(ID)
missing_measurement2

gludata_C57 <- gludata_C57 %>%
  filter(ID != 7860) %>%  # remove animals with only one measurement
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    n_measurement = row_number(),
    STATUS = case_when(
      n_measurement == 1 ~ "peak obesity",
      n_measurement == 2 & GROUP == "ad lib" ~ "BW maintenance",  # correctly rename for ad lib
      n_measurement == 2 & GROUP == "restricted" ~ "BW loss",
      TRUE ~ "other"
    ),
    STATUS = factor(STATUS, levels = c("peak obesity", "BW maintenance", "BW loss"))
  ) %>%
  ungroup() %>%
  mutate(
    highlight = case_when(
      # Restricted group: glucose increased
      GROUP == "restricted" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) > 0 ~ "Restricted Increase",
      # Ad lib group: glucose decreased
      GROUP == "ad lib" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) < 0 ~ "Ad lib Decrease",
      TRUE ~ "Other"
    )
  )


gludata_C57 <- gludata_C57 %>%
  group_by(ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    highlight = case_when(
      # Restricted group: glucose increased
      GROUP == "restricted" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) > 0 ~ "Restricted Increase",
      # Ad lib group: glucose decreased
      GROUP == "ad lib" & FASTED_GLU_mg_dL - lag(FASTED_GLU_mg_dL) < 0 ~ "Ad lib Decrease",
      TRUE ~ "Other"
    )
  ) %>%
  ungroup()

# Remove unused STATUS levels per GROUP
gludata_C57_plot <- gludata_C57 %>%
  group_by(GROUP) %>%
  filter(STATUS %in% unique(STATUS)) %>%   # keep only STATUS that exist in this group
  ungroup() %>%
  mutate(STATUS = droplevels(STATUS))      # drop unused factor levels

# Align points in the center
gludata_C57_plot_centered <- gludata_C57_plot %>%
  mutate(x_pos = as.numeric(STATUS))  # convert factor to numeric for exact positions

# Then plot 
plot <- ggplot(gludata_C57_plot, aes(x = STATUS, y = FASTED_GLU_mg_dL, fill = STATUS)) + 
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) + 
  geom_point(aes(color = highlight), alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) + 
  geom_text(aes(label = ID), size = 3, vjust = -0.5)+
  scale_color_manual(values = c( "Restricted Increase" = "red", "Ad lib Decrease" = "blue", "Other" = "black" )) + 
  theme_minimal() + labs(x = NULL, y = "Fasted glucose (mg/dL)") +
  facet_wrap(~GROUP, scales = "free_x") + 
  theme(legend.position = "none") 
plot

#conclusion
#it seems fasted glucose in blood did not stabilizes before animals loss weight
#this happens in both NZO and C57 mice males and females
#we can split the population in some "responders to food restriction " animals 
#(those which decreased fasted blood glucose after restriction) and "non responders to food restriction"
#(those animals which did not decreased fasted blood glucose after restriction)

#non responders NZO (F): 3708, 3720, 3721, 3722, 3725, 3728, 3729
#non responders c57 (M): 7861, 7863
#non responders c57 (F): 7872, 7877, 7878

#there is some animals in the ad lib group that during the BW maintenance it seems decreased the fasted glucose
#which means that probably they are already diabetic, we can also divide the population in 

#ad lib diabetic: 3718, 3719
  