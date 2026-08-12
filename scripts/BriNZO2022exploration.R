# This script aims to explore changes in body composition and behavior in females NZO after exposure to RTIOXA47 for 4 weeks

#libraries----
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr) 
library(purrr)
library(Hmisc)
library(lme4)
library(lmerTest)
library(emmeans)
library(pracma)
library(lubridate)
library(stringr)
library(forcats)
library(patchwork)
library(ggpattern)
library(car)
library(broom) 
library(rstatix)

#BODY WEIGHT (BW) ANALYSIS----
## These NZO mice were on chow ----

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT ==21)  %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(50,51,52,53,56,57,61,66,67,68,70) ~ "vehicle", 
      ID %in% c(49,54,55,58,59,60,62,63,64,65,69,71) ~ "RTI_47")
  ) %>%
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE) %>% 
  mutate(
    STATUS = case_when(
       DATE == as.Date("2022-10-28") ~ "baseline",
        DATE == as.Date("2022-12-02") ~ "end"
  ) )

BW_data  %>% 
  group_by(SEX,STRAIN,DRUG) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf)
  
  
BW_data_2 <- BW_data %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) 

BW_data_2 %>% 
  group_by(day_rel) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great this is consistent


BW_summary <- BW_data_2 %>%
  group_by(day_rel,SEX,STRAIN) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) 


# BODY COMPOSITION ANALYSIS----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT ==21)  %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(50,51,52,53,56,57,61,66,67,68,70) ~ "vehicle", 
      ID %in% c(49,54,55,58,59,60,62,63,64,65,69,71) ~ "RTI_47")
  )
  
echoMRI_data <- echoMRI_data %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(fat_perc = (Fat/Weight)*100,
         lean_perc = (Lean/Weight)*100) %>% 
  mutate(
    STATUS = case_when(
      n_measurement == 1 ~ "start",
      n_measurement == 2 ~ "end"
    ) )

echoMRI_data %>% 
  group_by(SEX,DRUG,STRAIN,n_measurement) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

delta_bodycomp <- echoMRI_data %>%
  select(
    ID, SEX, STRAIN, DRUG, STATUS,
    adiposity_index, Fat, Lean, fat_perc, lean_perc
  ) %>%
  filter(STATUS %in% c("start", "end")) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = c(
      adiposity_index,
      Fat,
      Lean,
      fat_perc,
      lean_perc
    ),
    names_glue = "{.value}_{STATUS}"
  ) %>%
  mutate(
    delta_ai        = adiposity_index_end - adiposity_index_start,
    delta_fat       = Fat_end - Fat_start,
    delta_lean      = Lean_end - Lean_start,
    delta_fat_perc  = fat_perc_end - fat_perc_start,
    delta_lean_perc = lean_perc_end - lean_perc_start
  )

delta_bodycomp

delta_bodycomp %>% 
  group_by(SEX,DRUG,STRAIN) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 


### ADIPOSITY INDEX----
### plot A: Adiposity before and after chronic injections of RTIOXA 47 ----

AI_summary <- echoMRI_data %>%
  group_by(STATUS, STRAIN, SEX, DRUG) %>%
  summarise(
    mean_ai = mean(adiposity_index, na.rm = TRUE),
    sem_ai  = sd(adiposity_index, na.rm = TRUE) / sqrt(sum(!is.na(adiposity_index))),
    n = sum(!is.na(adiposity_index)),
    .groups = "drop"
  ) %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    GROUP = paste(DRUG, STATUS, sep = "_")
  )

plot_ai <- ggplot(
  AI_summary,
  aes(
    x = DRUG,
    y = mean_ai,
    fill = GROUP,
    group = STATUS
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(
      ymin = mean_ai - sem_ai,
      ymax = mean_ai + sem_ai
    ),
    position = position_dodge(width = 0.8),
    width = 0.15
  ) +
  scale_fill_manual(
    values = c(
      "vehicle_start" = "white",
      "vehicle_end" = "grey85",
      "RTI_47_start" = "#FDD0A2",
      "RTI_47_end" = "#E67E22"
    ),
    labels = c(
      "vehicle_start" = "Vehicle start",
      "vehicle_end" = "Vehicle end",
      "RTI_47_start" = "RTI-47 start",
      "RTI_47_end" = "RTI-47 end"
    )
  ) +
  facet_wrap(~ STRAIN) +
  labs(
    x = NULL,
    y = "Adiposity index (fat/lean mass)",
    fill = NULL
  ) +
  theme_classic(base_size = 14)

plot_ai


### plot B: Delta adiposity after chronic injections of RTIOXA 47 ----

AIdelta_summary <- delta_bodycomp %>%
  group_by(STRAIN, DRUG) %>%
  summarise(
    mean_aidelta = mean(delta_ai, na.rm = TRUE),
    sem_aidelta  = sd(delta_ai, na.rm = TRUE) / sqrt(sum(!is.na(delta_ai))),
    n = sum(!is.na(delta_ai)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )

plot_aidelta <- ggplot(
  AIdelta_summary,
  aes(
    x = DRUG,
    y = mean_aidelta,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(
      ymin = mean_aidelta - sem_aidelta,
      ymax = mean_aidelta + sem_aidelta
    ),
    width = 0.15
  ) +
  scale_fill_manual(
    values = c(
      "vehicle" = "white",
      "RTI_47" = "#E67E22"
    ),
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  facet_wrap(~ STRAIN) +
  labs(
    x = NULL,
    y = "Change in adiposity index (end - start)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )

plot_aidelta

# Final igure adiposity index ----
# Find common y-axis limits across both plots
y_min <- min(
  AI_summary$mean_ai - AI_summary$sem_ai,
  AIdelta_summary$mean_aidelta - AIdelta_summary$sem_aidelta,
  na.rm = TRUE
)

y_max <- max(
  AI_summary$mean_ai + AI_summary$sem_ai,
  AIdelta_summary$mean_aidelta + AIdelta_summary$sem_aidelta,
  na.rm = TRUE
)

# Apply the same y-axis limits
plot_ai <- plot_ai +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "A")

plot_aidelta <- plot_aidelta +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "B")

combined_plot <- plot_ai | plot_aidelta

combined_plot

##STATS

AI_delta_models <- delta_bodycomp %>%
  filter(!is.na(delta_ai)) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_ai ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

AI_delta_contrasts <- AI_delta_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

AI_delta_contrasts

