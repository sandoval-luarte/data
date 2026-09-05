# This script aims to explore changes in body composition and behavior in C57BL6J and NZO mice after exposure to RTIOXA47 for 4 weeks

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
## These c57bl6j and nzo mice were on chow ----

# asignation to the drugs based on food intake calculation sheet csv from Bri original data

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(8, 20, 21)) %>%
#  filter(STRAIN == "C57BL/6J") %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(
        1,3,5,7,9,11,13,15,17,19,21,23, # cohort 8
        25,27,30,31,33,35,37,39,41,43,45,48, # cohort 20
        49,51,54,55,57,59,61,63,65,67,69 # cohort 21
      ) ~ "vehicle",
      
      ID %in% c(
        2,4,6,8,10,12,14,16,18,20,22,24, # cohort 8
        26,28,29,32,34,36,38,40,42,44,46,47, # cohort 20
        50,52,53,56,58,60,62,64,66,68,70,71 # cohort 21
      ) ~ "RTI_47"
    )
  ) %>%                        
  mutate(
    DATE = ymd(DATE)
  ) %>% 
  arrange(DATE) %>% 
  mutate(
    STATUS = case_when(
      COHORT == 8 & DATE == as.Date("2022-03-04") ~ "start",
      COHORT == 8 & DATE == as.Date("2022-04-08") ~ "end",
      COHORT == 20 & DATE == as.Date("2022-04-22") ~ "start",
      COHORT == 20 & DATE == as.Date("2022-05-27") ~ "end",
      COHORT == 21 & DATE == as.Date("2022-10-28") ~ "start",
      COHORT == 21 & DATE == as.Date("2022-12-02") ~ "end",
      TRUE ~ NA_character_
    )
  ) %>% 
ungroup()
   
BW_data  %>% 
  group_by(SEX,COHORT,STATUS,STRAIN) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf)
  
  
BW_data_2 <- BW_data %>% 
  group_by(COHORT, ID) %>% 
  mutate(
    start_date = DATE[STATUS %in% "start"][1],
    day_rel = as.numeric(DATE - start_date),
    bw_start = BW[STATUS %in% "start"][1],
    bw_rel = 100 * (BW - bw_start) / bw_start
  ) %>% 
  ungroup() %>% 
  filter(day_rel >= 0, day_rel < 40)

BW_data_2 %>% 
  group_by(day_rel,COHORT,STATUS,SEX,STRAIN) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great this is consistent we have 5 weeks of RTIOXA-47 injections

BW_data_2 <- BW_data_2 %>%
  mutate(
    day_rel_factor = factor(day_rel, levels = sort(unique(day_rel)))
  )

BW_summary <- BW_data_2 %>%
  group_by(day_rel, day_rel_factor, SEX, STRAIN, DRUG) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW  = sd(BW, na.rm = TRUE) / sqrt(sum(!is.na(BW))),
    n = sum(!is.na(BW)),
    .groups = "drop"
  )


plot_BW <- ggplot() +
  geom_line(data = BW_data_2, #follow the IDs
            aes(x = day_rel_factor,y = BW,group = ID,color = DRUG),alpha = 0.3,linewidth = 0.5) +
  geom_point(data = BW_data_2,#Individual IDs
    aes(x = day_rel_factor,y = BW,group = ID,color = DRUG),alpha = 0.5,size = 1.5) +
  geom_line(data = BW_summary,  # Group mean
    aes(x = day_rel_factor,y = mean_BW,group = DRUG,
      color = DRUG),linewidth = 1.2) +
  geom_point( # Mean points
    data = BW_summary,
    aes(x = day_rel_factor,
      y = mean_BW,
      color = DRUG
    ),
    size = 3
  ) +
  
  # SEM
  geom_errorbar(
    data = BW_summary,
    aes(
      x = day_rel_factor,
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW,
      color = DRUG
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  facet_wrap(~SEX) +
  
  labs(
    x = "Day",
    y = "Body weight (g)",
    color = NULL
  ) +
  geom_text(
    data = BW_data_2,
    aes(
      x = day_rel_factor,
      y = BW,
      label = ID,
      color = DRUG
    ),
    vjust = -0.7,
    size = 3,
    show.legend = FALSE
  ) 

plot_BW

#stats----
BW_start_end <- BW_data %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, SEX, DRUG, STATUS, BW)

BW_start_end %>% 
  group_by(SEX,STATUS) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) #ok great this is consistent we have 12 males and 24 females in this study

BW_ttest <- BW_start_end %>% 
  group_by(SEX, STATUS) %>% 
  t_test(
    BW ~ DRUG,
    paired = FALSE
  ) %>% 
  add_significance()

BW_ttest

#Does RTI_47 reduce BW during the 4-week treatment?----
BW_change <- BW_data %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, SEX, DRUG, STATUS, BW) %>% 
  pivot_wider(
    names_from = STATUS,
    values_from = BW
  ) %>% 
  mutate(
    BW_change = end - start,
    BW_change_percent = 100 * (end - start) / start
  )

BW_change_ttest <- BW_change %>% 
  group_by(SEX) %>% 
  t_test(
    BW_change ~ DRUG,
    paired = FALSE
  ) %>% 
  add_significance()

BW_change_ttest #there was a non-significant trend toward reduced BW gain in RTI_47-treated males.


BW_ANCOVA <- BW_data %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, SEX, DRUG, STATUS, BW) %>% 
  pivot_wider(
    names_from = STATUS,
    values_from = BW
  ) %>% 
  drop_na(start, end)

BW_ANCOVA %>% 
  group_by(SEX, DRUG) %>% 
  summarise(
    mean_start = mean(start),
    sem_start = sd(start) / sqrt(n()),
    mean_end = mean(end),
    sem_end = sd(end) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ANCOVA_F <- lm(
  end ~ start + DRUG,
  data = BW_ANCOVA %>% filter(SEX == "F")
)

anova(ANCOVA_F)
summary(ANCOVA_F)

ANCOVA_F_interaction <- lm(
  end ~ start * DRUG,
  data = BW_ANCOVA %>% filter(SEX == "F")
)

anova(ANCOVA_F_interaction)

summary(ANCOVA_F_interaction)


ANCOVA_M <- lm(
  end ~ start + DRUG,
  data = BW_ANCOVA %>% filter(SEX == "M")
)

anova(ANCOVA_M)
summary(ANCOVA_M)

#FOOD INTAKE----

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT ==21)  %>% 
  filter(STRAIN =="C57BL/6J")  %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(50,51,52,53,56,57,61,66,67,68,70) ~ "vehicle", 
      ID %in% c(49,54,55,58,59,60,62,63,64,65,69,71) ~ "RTI_47")
  ) %>%
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    day_rel = as.integer(as.Date(DATE) - as.Date(first(DATE)))
  ) %>% 
  left_join(METABPA, by= "ID") %>% 
  select(
    -SEX.y,
    -DIET_FORMULA.y,
    -DIET_FORMULA.x,
    -COHORT.y
  ) %>% 
  rename(
    SEX = SEX.x,
    COHORT = COHORT.x
  ) %>% 
  mutate(
    FIcumulative = cumsum(corrected_intake_kcal)) %>% 
  ungroup() %>% 
  # filter(!ID ==9406) %>% #9406 has a  weird pattern in locomotion
  mutate(week_rel = day_rel / 7) %>%  
  mutate(week_rel = round( week_rel)) %>%  #Mice did not get fed on 2/9, staff misread calendar
  filter(week_rel<=18) 

FI_data   %>% 
  group_by(week_rel) %>%
  summarise(n_ID = n_distinct(ID)) #we have consistently 56 animals across all weeks.


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
  
  # Individual animals
  geom_point(
    data = echoMRI_data %>%
      filter(STATUS %in% c("start", "end")) %>%
      mutate(
        STATUS = factor(STATUS, levels = c("start", "end")),
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = adiposity_index,
      group = STATUS
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
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
  facet_wrap(~ STRAIN*SEX) +
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
  
  # Individual animals
  geom_jitter(
    data = delta_bodycomp,
    aes(
      x = DRUG,
      y = delta_ai
    ),
    inherit.aes = FALSE,
    width = 0.08,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
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
  facet_wrap(~ STRAIN*SEX) +
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

### plot C: Adiposity at baseline within each strain ----


AI_summary_bs <- echoMRI_data %>%
  filter(STATUS == "start") %>%
  group_by(STRAIN, SEX, DRUG) %>%
  summarise(
    mean_ai = mean(adiposity_index, na.rm = TRUE),
    sem_ai  = sd(adiposity_index, na.rm = TRUE) / 
      sqrt(sum(!is.na(adiposity_index))),
    n = sum(!is.na(adiposity_index)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )


plot_bs <- ggplot(
  AI_summary_bs,
  aes(
    x = DRUG,
    y = mean_ai,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals at baseline only
  geom_point(
    data = echoMRI_data %>%
      filter(
        STATUS == "start",
        !is.na(adiposity_index)
      ) %>%
      mutate(
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = adiposity_index
    ),
    position = position_jitter(
      width = 0.08
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_ai - sem_ai,
      ymax = mean_ai + sem_ai
    ),
    width = 0.15
  ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "white",
      "RTI_47" = "#FDD0A2"
    ),
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  
  facet_wrap(~ STRAIN * SEX) +
  
  labs(
    x = NULL,
    y = "Adiposity index at start",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14)

plot_bs

# Final figure adiposity index ----
# Find common y-axis limits across both plots
y_min <- min(
  AI_summary$mean_ai - AI_summary$sem_ai,
  AIdelta_summary$mean_aidelta - AIdelta_summary$sem_aidelta,
  AI_summary_bs$mean_ai - AI_summary_bs$sem_ai,
  na.rm = TRUE
)

y_max <- max(
  AI_summary$mean_ai + AI_summary$sem_ai,
  AIdelta_summary$mean_aidelta + AIdelta_summary$sem_aidelta,
  AI_summary_bs$mean_ai + AI_summary_bs$sem_ai,
  na.rm = TRUE
)

# Apply the same y-axis limits
plot_ai <- plot_ai +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "A")

plot_aidelta <- plot_aidelta +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "B")

plot_bs <- plot_bs +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "C")


combined_plot <- plot_ai | plot_aidelta| plot_bs

combined_plot

##STATS----

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

### Baseline comparison: Vehicle start vs RTI-47 start within each strain ----

AI_start_models <- echoMRI_data %>%
  filter(
    STATUS == "start",
    !is.na(adiposity_index),
    !is.na(DRUG)
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(adiposity_index ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

AI_start_contrasts <- AI_start_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

AI_start_contrasts
#Great. Based on these results, there is no statistically significant difference in baseline adiposiy index between Vehicle and RTI-47 animals within either strain


# LEAN % (lean mass/BW)----
### plot A: % lean mass before and after chronic injections of RTIOXA 47 ----

leanperc_summary <- echoMRI_data %>%
  group_by(STATUS, STRAIN, SEX, DRUG) %>%
  summarise(
    mean_leanperc = mean(lean_perc, na.rm = TRUE),
    sem_leanperc  = sd(lean_perc, na.rm = TRUE) / sqrt(sum(!is.na(lean_perc))),
    n = sum(!is.na(lean_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    GROUP = paste(DRUG, STATUS, sep = "_")
  )

plot_leanperc <- ggplot(
  leanperc_summary,
  aes(
    x = DRUG,
    y = mean_leanperc,
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
  
  # Individual animals
  geom_point(
    data = echoMRI_data %>%
      filter(STATUS %in% c("start", "end")) %>%
      mutate(
        STATUS = factor(STATUS, levels = c("start", "end")),
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = lean_perc,
      group = STATUS
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_leanperc - sem_leanperc,
      ymax = mean_leanperc + sem_leanperc
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
  facet_wrap(~ STRAIN*SEX) +
  labs(
    x = NULL,
    y = "Lean % (lean mass/BW)",
    fill = NULL
  ) +
  theme_classic(base_size = 14)

plot_leanperc

### plot B: Delta lean % after chronic injections of RTIOXA 47 ----

leanpercdelta_summary <- delta_bodycomp %>%
  group_by(STRAIN, DRUG) %>%
  summarise(
    mean_leanpercdelta = mean(delta_lean_perc, na.rm = TRUE),
    sem_leanpercdelta  = sd(delta_lean_perc, na.rm = TRUE) / sqrt(sum(!is.na(delta_lean_perc))),
    n = sum(!is.na(delta_lean_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )

plot_leanpercdelta <- ggplot(
  leanpercdelta_summary,
  aes(
    x = DRUG,
    y = mean_leanpercdelta,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals
  geom_jitter(
    data = delta_bodycomp,
    aes(
      x = DRUG,
      y = delta_lean_perc
    ),
    inherit.aes = FALSE,
    width = 0.08,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_leanpercdelta - sem_leanpercdelta,
      ymax = mean_leanpercdelta + sem_leanpercdelta
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
  facet_wrap(~ STRAIN*SEX) +
  labs(
    x = NULL,
    y = "Change in lean% (end - start)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )

plot_leanpercdelta

### plot C: lean% at baseline within each strain ----

leanperc_summary_bs <- echoMRI_data %>%
  filter(STATUS == "start") %>%
  group_by(STRAIN, SEX, DRUG) %>%
  summarise(
    mean_leanperc = mean(lean_perc, na.rm = TRUE),
    sem_leanperc  = sd(lean_perc, na.rm = TRUE) / 
      sqrt(sum(!is.na(lean_perc))),
    n = sum(!is.na(lean_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )


plot_leanpercbs <- ggplot(
  leanperc_summary_bs,
  aes(
    x = DRUG,
    y = mean_leanperc,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals at baseline only
  geom_point(
    data = echoMRI_data %>%
      filter(
        STATUS == "start",
        !is.na(lean_perc)
      ) %>%
      mutate(
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = lean_perc
    ),
    position = position_jitter(
      width = 0.08
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_leanperc - sem_leanperc,
      ymax = mean_leanperc + sem_leanperc
    ),
    width = 0.15
  ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "white",
      "RTI_47" = "#FDD0A2"
    ),
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  
  facet_wrap(~ STRAIN * SEX) +
  
  labs(
    x = NULL,
    y = "Lean % at start",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14)

plot_leanpercbs

# Final figure lean % ----
# Find common y-axis limits across both plots
y_min <- min(
  leanperc_summary$mean_leanperc - leanperc_summary$sem_leanperc,
  leanpercdelta_summary$mean_leanpercdelta - leanpercdelta_summary$sem_leanpercdelta,
  leanperc_summary_bs$mean_leanperc - leanperc_summary_bs$sem_leanperc,
  na.rm = TRUE
)

y_max <- max(
  leanperc_summary$mean_leanperc + leanperc_summary$sem_leanperc,
  leanpercdelta_summary$mean_leanpercdelta + leanpercdelta_summary$sem_leanpercdelta,
  leanperc_summary_bs$mean_leanperc + leanperc_summary_bs$sem_leanperc,
  na.rm = TRUE
)

# Apply the same y-axis limits
plot_leanperc  <- plot_leanperc  +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "A")

plot_leanpercdelta <- plot_leanpercdelta +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "B")

plot_leanpercbs  <- plot_leanpercbs  +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "C")


combined_plot <- plot_leanperc | plot_leanpercdelta | plot_leanpercbs

combined_plot

##STATS----

leanperc_delta_models <- delta_bodycomp %>%
  filter(!is.na(delta_lean_perc)) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_lean_perc ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

leanperc_delta_contrasts <- leanperc_delta_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

leanperc_delta_contrasts

### Baseline comparison: Vehicle start vs RTI-47 start within each strain ----

leanperc_start_models <- echoMRI_data %>%
  filter(
    STATUS == "start",
    !is.na(lean_perc),
    !is.na(DRUG)
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(lean_perc ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

leanperc_start_contrasts <- leanperc_start_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

leanperc_start_contrasts
#Great. Based on these results, there is no statistically significant difference in baseline % lean mass between Vehicle and RTI-47 animals within either strain

#now the interesting question is for panel B (i.e is the % lean mass chanfe within the RTI-47 group within each strain different from zero?)

### RTI-47: Is the change in lean % significantly different from zero? ----

leanperc_RTI_delta_models <- delta_bodycomp %>%
  filter(
    DRUG == "RTI_47",
    !is.na(delta_lean_perc)
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_lean_perc ~ 1, data = .x)),
    emmean = map(model, ~ emmeans(.x, ~ 1)),
    test_zero = map(emmean, ~ test(.x, null = 0))
  )

leanperc_RTI_delta_zero <- leanperc_RTI_delta_models %>%
  select(STRAIN, test_zero) %>%
  mutate(
    test_zero = map(test_zero, ~ as.data.frame(.x))
  ) %>%
  unnest(test_zero)

leanperc_RTI_delta_zero

### Vehicle: Is the change in lean % significantly different from zero? ----

leanperc_vehicle_delta_models <- delta_bodycomp %>%
  filter(
    DRUG == "vehicle",
    !is.na(delta_lean_perc)
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_lean_perc ~ 1, data = .x)),
    emmean = map(model, ~ emmeans(.x, ~ 1)),
    test_zero = map(emmean, ~ test(.x, null = 0))
  )

leanperc_vehicle_delta_zero <- leanperc_vehicle_delta_models %>%
  select(STRAIN, test_zero) %>%
  mutate(
    test_zero = map(test_zero, ~ as.data.frame(.x))
  ) %>%
  unnest(test_zero)

leanperc_vehicle_delta_zero


#Vehicle: lean % decreased by 5.26 percentage points, and this change was significantly different from zero (p = 0.0073).
#RTI-47: lean % decreased by only 1.55 percentage points, and this change was not significantly different from zero (p = 0.124).
#So you can say:
#In NZO/HlLtJ mice, lean mass as a percentage of body weight significantly decreased from baseline in the vehicle group, whereas no significant change from baseline was detected in the RTI-47 group.
#In C57BL/6J mice, lean mass as a percentage of body weight significantly decreased from baseline in both the Vehicle and RTI-47 groups. However, the numerical decrease was smaller in RTI-47-treated mice (−2.96 percentage points) than in Vehicle-treated mice (−7.14 percentage points)


# FAT % (fat mass/BW)----
### plot A: % fat mass before and after chronic injections of RTIOXA 47 ----

fatperc_summary <- echoMRI_data %>%
  group_by(STATUS, STRAIN, SEX, DRUG) %>%
  summarise(
    mean_fatperc = mean(fat_perc, na.rm = TRUE),
    sem_fatperc  = sd(fat_perc, na.rm = TRUE) / sqrt(sum(!is.na(fat_perc))),
    n = sum(!is.na(fat_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    GROUP = paste(DRUG, STATUS, sep = "_")
  )

plot_fatperc <- ggplot(
  fatperc_summary,
  aes(
    x = DRUG,
    y = mean_fatperc,
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
  
  # Individual animals
  geom_point(
    data = echoMRI_data %>%
      filter(STATUS %in% c("start", "end")) %>%
      mutate(
        STATUS = factor(STATUS, levels = c("start", "end")),
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = fat_perc,
      group = STATUS
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_fatperc - sem_fatperc,
      ymax = mean_fatperc + sem_fatperc
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
  facet_wrap(~ STRAIN*SEX) +
  labs(
    x = NULL,
    y = "Fat % (fat mass/BW)",
    fill = NULL
  ) +
  theme_classic(base_size = 14)

plot_fatperc

### plot B: Delta fat % after chronic injections of RTIOXA 47 ----

fatpercdelta_summary <- delta_bodycomp %>%
  group_by(STRAIN, DRUG) %>%
  summarise(
    mean_fatpercdelta = mean(delta_fat_perc, na.rm = TRUE),
    sem_fatpercdelta  = sd(delta_fat_perc, na.rm = TRUE) / sqrt(sum(!is.na(delta_fat_perc))),
    n = sum(!is.na(delta_fat_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )

plot_fatpercdelta <- ggplot(
  fatpercdelta_summary,
  aes(
    x = DRUG,
    y = mean_fatpercdelta,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals
  geom_jitter(
    data = delta_bodycomp,
    aes(
      x = DRUG,
      y = delta_fat_perc
    ),
    inherit.aes = FALSE,
    width = 0.08,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_fatpercdelta - sem_fatpercdelta,
      ymax = mean_fatpercdelta + sem_fatpercdelta
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
  facet_wrap(~ STRAIN*SEX) +
  labs(
    x = NULL,
    y = "Change in fat% (end - start)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )

plot_fatpercdelta

### plot C: fat% at baseline within each strain ----

fatperc_summary_bs <- echoMRI_data %>%
  filter(STATUS == "start") %>%
  group_by(STRAIN, SEX, DRUG) %>%
  summarise(
    mean_fatperc = mean(fat_perc, na.rm = TRUE),
    sem_fatperc  = sd(fat_perc, na.rm = TRUE) / 
      sqrt(sum(!is.na(fat_perc))),
    n = sum(!is.na(fat_perc)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )


plot_fatpercbs <- ggplot(
  fatperc_summary_bs,
  aes(
    x = DRUG,
    y = mean_fatperc,
    fill = DRUG
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals at baseline only
  geom_point(
    data = echoMRI_data %>%
      filter(
        STATUS == "start",
        !is.na(fat_perc)
      ) %>%
      mutate(
        DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
      ),
    aes(
      x = DRUG,
      y = fat_perc
    ),
    position = position_jitter(
      width = 0.08
    ),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_fatperc - sem_fatperc,
      ymax = mean_fatperc + sem_fatperc
    ),
    width = 0.15
  ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "white",
      "RTI_47" = "#FDD0A2"
    ),
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  
  facet_wrap(~ STRAIN * SEX) +
  
  labs(
    x = NULL,
    y = "fat % at start",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14)

plot_fatpercbs

# Final figure fat % ----
# Find common y-axis limits across both plots
y_min <- min(
  fatperc_summary$mean_fatperc - fatperc_summary$sem_fatperc,
  fatpercdelta_summary$mean_fatpercdelta - fatpercdelta_summary$sem_fatpercdelta,
  fatperc_summary_bs$mean_fatperc - fatperc_summary_bs$sem_fatperc,
  na.rm = TRUE
)

y_max <- max(
  fatperc_summary$mean_fatperc + fatperc_summary$sem_fatperc,
  fatpercdelta_summary$mean_fatpercdelta + fatpercdelta_summary$sem_fatpercdelta,
  fatperc_summary_bs$mean_fatperc + fatperc_summary_bs$sem_fatperc,
  na.rm = TRUE
)

# Apply the same y-axis limits
plot_fatperc  <- plot_fatperc  +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "A")

plot_fatpercdelta <- plot_fatpercdelta +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "B")

plot_fatpercbs  <- plot_fatpercbs  +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(tag = "C")


combined_plot <- plot_fatperc | plot_fatpercdelta | plot_fatpercbs

combined_plot

##STATS----

fatperc_delta_models <- delta_bodycomp %>%
  filter(!is.na(delta_fat_perc)) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_fat_perc ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

fatperc_delta_contrasts <- fatperc_delta_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

fatperc_delta_contrasts

### Baseline comparison: Vehicle start vs RTI-47 start within each strain ----

fatperc_start_models <- echoMRI_data %>%
  filter(
    STATUS == "start",
    !is.na(fat_perc),
    !is.na(DRUG)
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(fat_perc ~ DRUG, data = .x)),
    emmeans = map(model, ~ emmeans(.x, ~ DRUG)),
    contrast = map(emmeans, ~ pairs(.x))
  )

fatperc_start_contrasts <- fatperc_start_models %>%
  select(STRAIN, contrast) %>%
  mutate(
    contrast = map(contrast, ~ as.data.frame(.x))
  ) %>%
  unnest(contrast)

fatperc_start_contrasts
#Great. Based on these results, there is no statistically significant difference in baseline % fat mass between Vehicle and RTI-47 animals within either strain

#now the interesting question is for panel B (i.e is the % fat mass change within the RTI-47 group within each strain different from zero?)

### RTI-47: Is the change in fat % significantly different from zero? ----

fatperc_RTI_delta_models <- delta_bodycomp %>%
  filter(
    DRUG == "RTI_47",
    !is.na(delta_fat_perc)
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_fat_perc ~ 1, data = .x)),
    emmean = map(model, ~ emmeans(.x, ~ 1)),
    test_zero = map(emmean, ~ test(.x, null = 0))
  )

fatperc_RTI_delta_zero <- fatperc_RTI_delta_models %>%
  select(STRAIN, test_zero) %>%
  mutate(
    test_zero = map(test_zero, ~ as.data.frame(.x))
  ) %>%
  unnest(test_zero)

fatperc_RTI_delta_zero

### Vehicle: Is the change in fat % significantly different from zero? ----

fatperc_vehicle_delta_models <- delta_bodycomp %>%
  filter(
    DRUG == "vehicle",
    !is.na(delta_fat_perc)
  ) %>%
  group_by(STRAIN) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(delta_fat_perc ~ 1, data = .x)),
    emmean = map(model, ~ emmeans(.x, ~ 1)),
    test_zero = map(emmean, ~ test(.x, null = 0))
  )
fatperc_vehicle_delta_zero <- fatperc_vehicle_delta_models %>%
  select(STRAIN, test_zero) %>%
  mutate(
    test_zero = map(test_zero, ~ as.data.frame(.x))
  ) %>%
  unnest(test_zero)
fatperc_vehicle_delta_zero

#So descriptively, both strains accumulated less relative fat under RTI-OXA-47 than under vehicle, particularly NZO, where the RTI group was essentially unchanged

