#we aim to evaluate changes in blood glucose levels in NZO and C57 mice during peak obesity
#and BW loss
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
  mutate(STATUS= factor(STATUS, levels = c("peak obesity", "BW loss"))) #%>% 
 # filter(STATUS != "other")


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

gludata_NZO %>% 
  group_by(n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) 

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

# Then plot of fasted glucose in peak obesity and bw maintenance/loss 
plot <- ggplot(gludata_NZO_plot, aes(x = STATUS, y = FASTED_GLU_mg_dL, fill = STATUS)) + 
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) + 
  geom_point(aes(color = highlight), alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) + 
 # geom_text(aes(label = ID), size = 3, vjust = -0.5)+
scale_color_manual(values = c( "Restricted Increase" = "red", "Ad lib Decrease" = "blue", "Other" = "black" )) + 
  theme_minimal() + labs(x = NULL, y = "Fasted glucose (mg/dL)") +
  facet_wrap(~GROUP, scales = "free_x") + 
  theme(legend.position = "none") 
plot

# Ensure numeric columns
gludata_NZO_plot <- gludata_NZO_plot %>%
  mutate(
    BODY_WEIGHT_G = as.numeric(BODY_WEIGHT_G),
    FASTED_GLU_mg_dL = as.numeric(FASTED_GLU_mg_dL)
  )

# Plot Body weight
plot_BW <- ggplot(gludata_NZO_plot, aes(x = STATUS, y = BODY_WEIGHT_G, fill = STATUS)) + 
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) + 
  geom_point(aes(color = highlight), alpha = 0.7, size = 2, 
             position = position_jitter(width = 0.15)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) + 
  scale_color_manual(values = c(
    "Restricted Increase" = "red", 
    "Ad lib Decrease" = "blue", 
    "Other" = "black"
  )) + 
  theme_minimal() + 
  labs(x = NULL, y = "Body weight (grams)") +
  facet_wrap(~GROUP, scales = "free_x") + 
  theme(legend.position = "none")

plot_BW

# Select paired data
adlib_data <- gludata_NZO_plot %>%
  filter(GROUP == "ad lib", STATUS %in% c("peak obesity", "BW maintenance")) %>%
  select(ID, STATUS, BODY_WEIGHT_G) %>%
  pivot_wider(names_from = STATUS, values_from = BODY_WEIGHT_G)

# Paired t-test
t.test(adlib_data$`BW maintenance`, adlib_data$`peak obesity`, paired = TRUE)

restricted_data <- gludata_NZO_plot %>%
  filter(GROUP == "restricted", STATUS %in% c("peak obesity", "BW loss")) %>%
  select(ID, STATUS, BODY_WEIGHT_G) %>%
  pivot_wider(names_from = STATUS, values_from = BODY_WEIGHT_G)

t.test(restricted_data$`BW loss`, restricted_data$`peak obesity`, paired = TRUE)


adlib_glu <- gludata_NZO_plot %>%
  filter(GROUP == "ad lib", STATUS %in% c("peak obesity", "BW maintenance")) %>%
  select(ID, STATUS, FASTED_GLU_mg_dL) %>%
  pivot_wider(names_from = STATUS, values_from = FASTED_GLU_mg_dL)

# Paired t-test
t.test(adlib_glu$`BW maintenance`, adlib_glu$`peak obesity`, paired = TRUE)

restricted_glu <- gludata_NZO_plot %>%
  filter(GROUP == "restricted", STATUS %in% c("peak obesity", "BW loss")) %>%
  select(ID, STATUS, FASTED_GLU_mg_dL) %>%
  pivot_wider(names_from = STATUS, values_from = FASTED_GLU_mg_dL)

# Paired t-test
t.test(restricted_glu$`BW loss`, restricted_glu$`peak obesity`, paired = TRUE)


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
  ) %>% 
  filter(!(STRAIN == "C57BL6/J" & DIET_FORMULA == "D12450Ki"))


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

# Ensure numeric columns
gludata_C57_plot <- gludata_C57_plot %>%
  mutate(
    BODY_WEIGHT_G = as.numeric(BODY_WEIGHT_G),
    FASTED_GLU_mg_dL = as.numeric(FASTED_GLU_mg_dL)
  )

# Plot Body weight
plot_BW <- ggplot(gludata_C57_plot, aes(x = STATUS, y = BODY_WEIGHT_G, fill = STATUS)) + 
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) + 
  geom_point(aes(color = highlight), alpha = 0.7, size = 2, 
             position = position_jitter(width = 0.15)) +
  geom_line(aes(group = ID), color = "gray50", alpha = 0.5) + 
  scale_color_manual(values = c(
    "Restricted Increase" = "red", 
    "Ad lib Decrease" = "blue", 
    "Other" = "black"
  )) + 
  theme_minimal() + 
  labs(x = NULL, y = "Body weight (grams)") +
  facet_wrap(~GROUP, scales = "free_x") + 
  theme(legend.position = "none")+
  geom_text(aes(label = ID), size = 3, vjust = -0.5)

plot_BW

# Select paired data
adlib_data <- gludata_NZO_plot %>%
  filter(GROUP == "ad lib", STATUS %in% c("peak obesity", "BW maintenance")) %>%
  select(ID, STATUS, BODY_WEIGHT_G) %>%
  pivot_wider(names_from = STATUS, values_from = BODY_WEIGHT_G)

# Paired t-test
t.test(adlib_data$`BW maintenance`, adlib_data$`peak obesity`, paired = TRUE)

restricted_data <- gludata_NZO_plot %>%
  filter(GROUP == "restricted", STATUS %in% c("peak obesity", "BW loss")) %>%
  select(ID, STATUS, BODY_WEIGHT_G) %>%
  pivot_wider(names_from = STATUS, values_from = BODY_WEIGHT_G)

t.test(restricted_data$`BW loss`, restricted_data$`peak obesity`, paired = TRUE)


adlib_glu <- gludata_NZO_plot %>%
  filter(GROUP == "ad lib", STATUS %in% c("peak obesity", "BW maintenance")) %>%
  select(ID, STATUS, FASTED_GLU_mg_dL) %>%
  pivot_wider(names_from = STATUS, values_from = FASTED_GLU_mg_dL)

# Paired t-test
t.test(adlib_glu$`BW maintenance`, adlib_glu$`peak obesity`, paired = TRUE)

restricted_glu <- gludata_NZO_plot %>%
  filter(GROUP == "restricted", STATUS %in% c("peak obesity", "BW loss")) %>%
  select(ID, STATUS, FASTED_GLU_mg_dL) %>%
  pivot_wider(names_from = STATUS, values_from = FASTED_GLU_mg_dL)

# Paired t-test
t.test(restricted_glu$`BW loss`, restricted_glu$`peak obesity`, paired = TRUE)


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


#I would like to know now if there is a correlation between the delta fasted glucose between peak obesity 
#and bw loss with delta BW between the two stages

gludata_wide <- gludata_NZO %>% 
  distinct(ID, STATUS, .keep_all = TRUE) %>% 
  pivot_wider( names_from = STATUS, values_from = c(BODY_WEIGHT_G, FASTED_GLU_mg_dL) ) 

write_csv(gludata_wide, "../data/gludata_wide.csv")

gludata_wide_modified <- read_csv("~/Documents/GitHub/data/data/gludata_wide_modified.csv") %>% 
  mutate(DATE = lubridate::mdy(DATE)) 
View(gludata_wide_modified)

gludata_wide_modified <- gludata_wide_modified %>% 
  mutate(across(starts_with("BODY_WEIGHT_G"), ~ as.numeric(.)))
  
gludata_wide_modified <- gludata_wide_modified %>% 
    mutate(across(starts_with("FASTED_GLU_mg_dL"), ~ as.numeric(.)))

sapply(gludata_wide_modified[, grepl("BODY_WEIGHT_G", names(gludata_wide_modified))], class)
sapply(gludata_wide_modified[, grepl("FASTED_GLU_mg_dL", names(gludata_wide_modified))], class)

gludata_NZO_delta <- gludata_wide_modified %>%
  mutate(
    delta_BW = case_when(
      GROUP == "restricted" ~ `BODY_WEIGHT_G_BW loss` - `BODY_WEIGHT_G_peak obesity`,
      GROUP == "ad lib"     ~ `BODY_WEIGHT_G_BW maintenance` - `BODY_WEIGHT_G_peak obesity`
    ),
    delta_fasted_glu = case_when(
      GROUP == "restricted" ~ `FASTED_GLU_mg_dL_BW loss` - `FASTED_GLU_mg_dL_peak obesity`,
      GROUP == "ad lib"     ~ `FASTED_GLU_mg_dL_BW maintenance` - `FASTED_GLU_mg_dL_peak obesity`
    ),
    SEX = ifelse(SEX %in% c("FALSE", FALSE), "F", SEX)
  ) %>% 
  select(ID,SEX,DATE,DIET_FORMULA,GROUP,highlight, delta_BW,delta_fasted_glu)

#correlation

gludata_NZO_delta %>%
  group_by(GROUP) %>%
  summarise(
    r = cor(delta_fasted_glu, delta_BW, use = "complete.obs")
  )


# Compute R² per group
r2_per_group <- gludata_NZO_delta %>%
  group_by(GROUP) %>%
  summarise(
    r = cor(delta_fasted_glu, delta_BW, use = "complete.obs"),
    R2 = r^2,
    # Determine label position for each group
    x_pos = max(delta_fasted_glu, na.rm = TRUE) - 5,  # slightly left from max
    y_pos = max(delta_BW, na.rm = TRUE) - 5           # slightly below max
  )

# Scatter plot with regression lines and R² labels
ggplot(gludata_NZO_delta, aes(x = delta_fasted_glu, y = delta_BW, color = GROUP)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(
    data = r2_per_group,
    aes(x = x_pos, y = y_pos, label = paste0("R² = ", round(R2, 2))),
    inherit.aes = FALSE,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    x = "Delta Fasted Glucose (mg/dL)",
    y = "Delta Body Weight (g)",
    title = "Relationship between Delta BW and Delta Fasted Glucose by GROUP"
  )

#In NZO females, BW loss under restriction reduces body weight significantly but does not
#strongly affect fasting glucose. In contrast, ad lib animals show a stronger BW–glucose
#link. This suggests that obesity-induced metabolic dysfunction is partially resistant 
#to weight loss, highlighting the concept of "metabolic memory" in these mice.

