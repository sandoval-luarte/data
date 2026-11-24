
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
library(ARTool)
library(ggrepel)

# Color palettes ----
drug_colors <- c(
  "vehicle" = "gray60",
  "RTI_47" = "#E67E22"
)

drug_labels <- c(
  "vehicle" = "Vehicle",
  "RTI_47" = "RTI-47"
)

#format plot ----
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

#RTIOXA 47 for 5 weeks 40 mg/kg ----
## BODY COMPOSITION----

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(8, 20, 21))  %>%
  mutate(
    DRUG = case_when(
      ID %in% c(
        1,3,5,7,9,11,13,15,17,19,21,23, # cohort 8
        25,27,30,31,33,35,37,39,41,43,45,48, # cohort 20
        50,51,52,53,56,57,61,66,67,68,70 # cohort 21
      ) ~ "vehicle",
      
      ID %in% c(
        2,4,6,8,10,12,14,16,18,20,22,24, # cohort 8
        26,28,29,32,34,36,38,40,42,44,46,47, # cohort 20
        49,54,55,58,59,60,62,63,64,65,69,71 # cohort 21
      ) ~ "RTI_47"
    )
  )    %>% 
  filter(COHORT %in% c(8, 21)) #seems like cohort 20 has issues with solubility of the drug but in Bri's notes nothing is noticed


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
  group_by(DRUG) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

delta_bodycomp <- echoMRI_data %>%
  select(
    ID, SEX, STRAIN, DRUG, STATUS,
    adiposity_index, Fat, Lean, fat_perc, lean_perc,COHORT
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
    delta_lean_perc = lean_perc_end - lean_perc_start,
    delta_ai_standar= (adiposity_index_end - adiposity_index_start)/adiposity_index_start,
    delta_fat_standar = (Fat_end - Fat_start)/Fat_start,
    delta_lean_standar = (Lean_end - Lean_start)/Lean_start,
    delta_fat_perc_standar  = (fat_perc_end - fat_perc_start)/fat_perc_start,
    delta_lean_perc_standar = (lean_perc_end - lean_perc_start)/lean_perc_start) %>% 
  ungroup()


delta_bodycomp

delta_bodycomp %>% 
  group_by(DRUG) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 


### 1. Fat mass ----

fatdelta_summary <- delta_bodycomp %>%
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    # Absolute fat change
    mean_fatdelta = mean(delta_fat, na.rm = TRUE),
    sem_fatdelta = sd(delta_fat, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_fat))),
    n_fatdelta = sum(!is.na(delta_fat)),
    
    # Fat percentage change
    mean_delta_fat_perc = mean(delta_fat_perc, na.rm = TRUE),
    sem_delta_fat_perc = sd(delta_fat_perc, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_fat_perc))),
    n_delta_fat_perc = sum(!is.na(delta_fat_perc)),
    
    # Standardized absolute fat change
    mean_delta_fat_standar = mean(delta_fat_standar, na.rm = TRUE),
    sem_delta_fat_standar = sd(delta_fat_standar, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_fat_standar))),
    n_delta_fat_standar = sum(!is.na(delta_fat_standar)),
    
    # Standardized fat percentage change
    mean_delta_fat_perc_standar = mean(delta_fat_perc_standar, na.rm = TRUE),
    sem_delta_fat_perc_standar = sd(delta_fat_perc_standar, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_fat_perc_standar))),
    n_delta_fat_perc_standar = sum(!is.na(delta_fat_perc_standar)),
    
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )


#### Fat mass plot ----

plot_fatdelta <- ggplot(
  fatdelta_summary,
  aes(
    x = DRUG,
    y = mean_fatdelta,
    fill = DRUG
  )
) +
  
  # Mean bars
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Individual animals
  geom_jitter(
    data = delta_bodycomp,
    aes(
      x = DRUG,
      y = delta_fat
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.7
    ),
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_fatdelta - sem_fatdelta,
      ymax = mean_fatdelta + sem_fatdelta
    ),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  
  # Colors
  scale_fill_manual(
    values = drug_colors,
    labels = drug_labels
  ) +
  
  scale_x_discrete(
    labels = drug_labels
  ) +
  
  labs(
    x = NULL,
    y = "Change in fat mass (end-start)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  ) +
  
  facet_wrap(~SEX*STRAIN)

plot_fatdelta

####Fat mass stats----

fatdelta_ttest <- delta_bodycomp %>%
  filter(!is.na(delta_fat), !is.na(DRUG)) %>%
  group_by(SEX,STRAIN) %>%
  rstatix::t_test(delta_fat ~ DRUG) %>%
  rstatix::add_significance()

fatdelta_ttest

### 2. Adiposity index ----

adipositydelta_summary <- delta_bodycomp %>%
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    # Absolute fat change
    mean_delta_adiposity = mean(delta_ai, na.rm = TRUE),
    sem_delta_adiposity = sd(delta_ai, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_ai))),
    n_aidelta = sum(!is.na(delta_ai)),
    
    # Standardized absolute adiposity index change
    mean_delta_adiposity_standar = mean(delta_ai_standar, na.rm = TRUE),
    sem_delta_adiposity_standar = sd(delta_ai_standar, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_ai_standar))),
    n_delta_ai_standar = sum(!is.na(delta_ai_standar)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )

### adiposity index plot ----

plot_adipositydelta <- ggplot(
  adipositydelta_summary,
  aes(
    x = DRUG,
    y = mean_delta_adiposity,
    fill = DRUG
  )
) +
  
  # Mean bars
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
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
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.7
    ),
    size = 2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_delta_adiposity - sem_delta_adiposity,
      ymax = mean_delta_adiposity + sem_delta_adiposity
    ),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  
  # Colors
  scale_fill_manual(
    values = drug_colors,
    labels = drug_labels
  ) +
  
  scale_x_discrete(
    labels = drug_labels
  ) +
  
  labs(
    x = NULL,
    y = " Adiposity index change (end-start)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  ) +
  
  facet_wrap(~SEX*STRAIN)

plot_adipositydelta

###adiposity index stats----

adipositydelta_ttest <- delta_bodycomp %>%
  filter(!is.na(delta_ai), !is.na(DRUG)) %>%
  group_by(SEX,STRAIN) %>%
  rstatix::t_test(delta_ai ~ DRUG) %>%
  rstatix::add_significance()

adipositydelta_ttest

#now we want to test is this value is different from zero
#i.e do these animals lose lean mass?

leandelta_absolute_zero_ttest <- delta_bodycomp %>%
  filter(!is.na(delta_lean), !is.na(DRUG)) %>%
  group_by(DRUG, SEX) %>%
  rstatix::t_test(
    delta_lean ~ 1,
    mu = 0
  ) %>%
  rstatix::add_significance()

leandelta_absolute_zero_ttest


###3.  Lean mass----

leandelta_summary <- delta_bodycomp %>%
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    # Absolute fat change
    mean_delta_lean = mean(delta_lean, na.rm = TRUE),
    sem_delta_lean = sd(delta_lean, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_lean))),
    n_leandelta = sum(!is.na(delta_lean)),
    
    # Lean percentage change
    mean_delta_lean_perc = mean(delta_lean_perc, na.rm = TRUE),
    sem_delta_lean_perc = sd(delta_lean_perc, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_lean_perc))),
    n_delta_lean_perc = sum(!is.na(delta_lean_perc)),
    
    # Standardized absolute lean change
    mean_delta_lean_standar = mean(delta_lean_standar, na.rm = TRUE),
    sem_delta_lean_standar = sd(delta_lean_standar, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_lean_standar))),
    n_delta_lean_standar = sum(!is.na(delta_lean_standar)),
    
    # Standardized lean percentage change
    mean_delta_lean_perc_standar = mean(delta_lean_perc_standar, na.rm = TRUE),
    sem_delta_lean_perc_standar = sd(delta_lean_perc_standar, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_lean_perc_standar))),
    n_delta_lean_perc_standar = sum(!is.na(delta_lean_perc_standar)),
    
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )


### lean mass plot  ----
#to start practising coding 
