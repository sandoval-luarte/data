# This script aims to explore:
#1 changes in body composition and glucose managing in C57BL6J and NZO mice after 5 weeks of IP injections of RTIOXA47 40 mg-kg
#2 changes in body composition in C57BL6J mice after 5 days of IP injections of RTIOXA43 10 mg-kg

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
## BW ----
### These c57 and nzo mice were on chow over the course of the study
### assignation to the drugs was based on .csv file from Bri original data

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(8, 20, 21)) %>%
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
ungroup() %>% 
  filter(DATE >= as.Date("2022-03-04")) %>% #cohort 8 has acclimation, it seems like cohort 20 an 21 did not
  filter(COHORT %in% c(8, 21)) #seems like cohort 20 has issues with solubility of the drug but in Bri's notes nothing is noticed
   
BW_data  %>% 
#  filter(SEX=="F") %>% 
 # filter(STRAIN == "NZO/HlLtJ") %>% 
  group_by(DRUG,SEX,STRAIN,ID,COHORT) %>% #so cohort 21 are just females from both strains
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf)
  
###1. BW longitudinal ----  
BW_data_2 <- BW_data %>% 
  group_by(ID,DRUG,SEX,STRAIN) %>% 
  mutate(
    start_date = DATE[STATUS == "start"][1],
    day_rel = as.numeric(DATE - start_date),
    week_rel = floor(day_rel / 7),
    bw_start = BW[STATUS == "start"][1],
    bw_rel = 100 * (BW - bw_start) / bw_start
  ) %>% 
  ungroup() %>% 
  filter(day_rel >= 0, day_rel < 40) %>%
  mutate(
    week_rel = factor(week_rel, levels = 0:5)
  ) #%>% 
  #filter(STRAIN == "C57BL/6J") %>% 
  #filter(SEX=="F") 

BW_data_2 %>% 
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    n_ID = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  print(n = Inf)


BW_summary <- BW_data_2 %>%
  group_by(week_rel,DRUG,SEX,STRAIN) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(sum(!is.na(BW))),
    n = sum(!is.na(BW)),
    .groups = "drop"
  )

### plot BW longitudinal ----
plot_BW <- ggplot() +
  geom_line(   # Individual mouse trajectories
    data = BW_data_2,
    aes(
      x = week_rel,
      y = BW,
      group = ID,
      color = DRUG
    ),
    alpha = 0.25,
    linewidth = 0.5
  ) +
geom_point(   # Individual measurements
    data = BW_data_2,
    aes(
      x = week_rel,
      y = BW,
      color = DRUG
    ),
    alpha = 0.35,
    size = 1.5
  ) +
  geom_line(   # Group mean
    data = BW_summary,
    aes(
      x = week_rel,
      y = mean_BW,
      group = DRUG,
      color = DRUG
    ),
    linewidth = 1.2
  ) +
  geom_point(  # Mean points
    data = BW_summary,
    aes(
      x = week_rel,
      y = mean_BW,
      color = DRUG
    ),
    size = 3
  ) +
  geom_errorbar(  # SEM
    data = BW_summary,
    aes(
      x = week_rel,
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW,
      color = DRUG
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  facet_wrap(~ SEX*STRAIN) +
  
  scale_color_manual(
    values = c(
      "vehicle" = "gray60",
      "RTI_47" = "#E67E22"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "0" = "0",
      "1" = "1",
      "2" = "2",
      "3" = "3",
      "4" = "4",
      "5" = "5"
    )
  ) +
  
  labs(
    x = "Week",
    y = "Body weight (g)",
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )+
  scale_color_manual(
    values = drug_colors,
    labels = drug_labels
  )

plot_BW


### stats BW longitudinal ----
#At each week, is body weight significantly different between RTI-47 and Vehicle?
# Linear mixed-effects model collapsed by SEX and STRAIN
BW_model <- lmer(
  BW ~ DRUG * week_rel + (1 | ID),
  data = BW_data_2,
  REML = FALSE
)

anova(BW_model)

BW_emmeans <- emmeans(
  BW_model,
  ~ DRUG | week_rel
)

BW_emmeans

BW_weekly_comparisons <- pairs(
  BW_emmeans,
  adjust = "holm"
)

BW_weekly_comparisons


# Linear mixed-effects model separated by SEX and STRAIN
# Make sure categorical variables are factors
BW_data_2 <- BW_data_2 %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    week_rel = factor(week_rel),
    SEX = factor(SEX),
    STRAIN = factor(STRAIN)
  )

# Use sum-to-zero contrasts for Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

# Full factorial linear mixed-effects model
BW_model <- lmer(
  BW ~ DRUG * week_rel * SEX * STRAIN + (1 | ID),
  data = BW_data_2,
  REML = FALSE
)

# Overall tests of fixed effects
anova(BW_model, type = 3, ddf = "Satterthwaite")


# Estimated marginal means
BW_emmeans <- emmeans(
  BW_model,
  ~ DRUG | week_rel * SEX * STRAIN
)

# Compare RTI-47 vs Vehicle within each week,
# sex, and strain
BW_contrasts <- pairs(
  BW_emmeans,
  by = c("week_rel", "SEX", "STRAIN"),
  adjust = "none"
)

# Extract results
BW_contrasts_df <- as.data.frame(BW_contrasts)

# Apply Holm correction across ALL treatment comparisons
BW_contrasts_df$p_holm <- p.adjust(
  BW_contrasts_df$p.value,
  method = "holm"
)


# Display the results
BW_contrasts_df %>%
  dplyr::select(
    SEX, STRAIN, week_rel, contrast,
    estimate, SE, df, t.ratio, p.value, p_holm
  ) %>%
  dplyr::arrange(SEX, STRAIN, week_rel) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


### 2. BW comparison (end-start) ----
# Prepare summary data
BW_plot_summary <- BW_start_end %>%
  group_by(DRUG,STATUS) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(sum(!is.na(BW))),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    STATUS = factor(STATUS, levels = c("start", "end")),
    DRUG_STATUS = factor(
      paste(DRUG, STATUS, sep = "_"),
      levels = c(
        "vehicle_start",
        "vehicle_end",
        "RTI_47_start",
        "RTI_47_end"
      )
    )
  )


# Prepare individual mouse data
BW_points <- BW_start_end %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    STATUS = factor(STATUS, levels = c("start", "end")),
    DRUG_STATUS = factor(
      paste(DRUG, STATUS, sep = "_"),
      levels = c(
        "vehicle_start",
        "vehicle_end",
        "RTI_47_start",
        "RTI_47_end"
      )
    )
  )


### plot BW comparison (end-start) ----
plot_BW_start_end <- ggplot(
  BW_plot_summary,
  aes(
    x = DRUG_STATUS,
    y = mean_BW,
    fill = DRUG
  )
) +
  
  # Mean bars
  geom_col(
    width = 0.7
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW
    ),
    width = 0.2,
    linewidth = 0.7
  ) +
  
  # Paired individual mice
  geom_line(
    data = BW_points,
    aes(
      x = DRUG_STATUS,
      y = BW,
      group = ID
    ),
    color = "black",
    alpha = 0.35,
    linewidth = 0.6
  ) +
  
  # Individual mouse points
  geom_point(
    data = BW_points,
    aes(
      x = DRUG_STATUS,
      y = BW
    ),
    color = "black",
    size = 2.2
  ) +
  
  # Mouse IDs
  #geom_text(
   # data = BW_points,
    #aes(
     # x = DRUG_STATUS,
      #y = BW,
      #label = ID
  #  ),
   # vjust = -0.7,
  #  size = 3,
   # color = "black"
#  ) +
  
  # Drug colors
  scale_fill_manual(
    values = drug_colors,
    labels = drug_labels
  ) +
  
  # X-axis
  scale_x_discrete(
    labels = c(
      "vehicle_start" = "Start",
      "vehicle_end" = "End",
      "RTI_47_start" = "Start",
      "RTI_47_end" = "End"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Body weight (g)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  )#+
  #facet_wrap(~STRAIN)

plot_BW_start_end

### stats BW comparison (end-start) ----
BW_start_end <- BW_data %>% 
  #filter(STRAIN == "C57BL/6J") %>% 
  # filter(SEX=="F") %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, DRUG, STATUS, BW,STRAIN,SEX)

BW_start_end %>% 
  group_by(DRUG) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

BW_ttest <- BW_start_end %>% 
  group_by(DRUG) %>% 
  t_test(
    BW ~ STATUS,
    paired = TRUE
  ) %>% 
  add_significance()

BW_ttest


###3. ALTERNATIVE analysis by HIGH or LOW BW at start----

# 1. Extract each mouse's starting body weight
BW_start <- BW_data %>%
  filter(STATUS == "start") %>%
  select(ID, COHORT, STRAIN, SEX, BW,DRUG) %>%
  rename(BW_start = BW)

# 2. Calculate the mean starting BW for each strain and sex
BW_start <- BW_start %>%
#  group_by(SEX) %>%
  mutate(
    mean_BW_start = mean(BW_start, na.rm = TRUE),
    
    BW_group = case_when(
      BW_start > mean_BW_start ~ "HIGH",
      BW_start < mean_BW_start ~ "LOW",
      TRUE ~ "EQUAL"
    )
  ) %>%
  ungroup()

# 3. Check the classification
BW_start %>%
  group_by(BW_group) %>%
  summarise(
    n_mice = n_distinct(ID),
    mean_start_BW = mean(BW_start, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# 4. Add classification to every measurement

BW_data <- BW_data %>%
  select(-any_of(c("BW_start", "mean_BW_start", "BW_group"))) %>%
  left_join(
    BW_start %>%
      select(
        ID, COHORT, DRUG,
        BW_start, mean_BW_start, BW_group
      ),
    by = c("ID", "COHORT", "DRUG")
  )

# Check that the classification is consistent for each mouse
BW_data %>%
  distinct(ID, COHORT, STRAIN, SEX, BW_start, mean_BW_start, BW_group,DRUG) %>%
  arrange(STRAIN, SEX, BW_group, BW_start) %>%
  print(n = Inf)


### STATS: PAIRED START-END BW ANALYSIS BY BW GROUP ----

BW_start_end <- BW_data %>%
  filter(STATUS %in% c("start", "end")) %>%
  select(
    ID, COHORT, DRUG, STATUS, BW,
    STRAIN, SEX, BW_group
  )

# Check sample size in each group
BW_start_end %>%
  distinct(ID, COHORT, DRUG, STRAIN, SEX, BW_group) %>%
  group_by(STRAIN, SEX, DRUG, BW_group) %>%
  summarise(
    n_mice = n_distinct(ID),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Paired t-test comparing start vs end
BW_ttest <- BW_start_end %>%
  filter(BW_group %in% c("HIGH", "LOW")) %>%
  group_by(DRUG, BW_group) %>%
  t_test(
    BW ~ STATUS,
    paired = TRUE
  ) %>%
  add_significance()

BW_ttest


### PLOT: START VS END BW BY DRUG AND BW GROUP ----

# Summarize BW for plotting
BW_plot_summary <- BW_start_end %>%
  filter(BW_group %in% c("HIGH", "LOW")) %>%
  group_by(DRUG, BW_group, STATUS) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(sum(!is.na(BW))),
    n = sum(!is.na(BW)),
    .groups = "drop"
  )

# Set the order of the x-axis and facets
BW_start_end <- BW_start_end %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    BW_group = factor(BW_group, levels = c("LOW", "HIGH"))
  )

BW_plot_summary <- BW_plot_summary %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    BW_group = factor(BW_group, levels = c("LOW", "HIGH"))
  )


### PAIRED DOT PLOT: START VS END BW ----

# Prepare data
BW_start_end_plot <- BW_start_end %>%
  filter(BW_group %in% c("HIGH", "LOW")) %>%
  mutate(
    STATUS = factor(STATUS, levels = c("start", "end")),
    BW_group = factor(BW_group, levels = c("HIGH", "LOW")),
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )
# Calculate mean and SEM
BW_plot_summary <- BW_start_end_plot %>%
  group_by(DRUG, BW_group, STATUS) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) /
      sqrt(sum(!is.na(BW))),
    n = sum(!is.na(BW)),
    .groups = "drop"
  )

# Create paired dot plot
plot_BW_start_end <- ggplot(
  BW_start_end_plot,
  aes(x = STATUS, y = BW)
) +
  
  # Connect measurements from the same mouse
  geom_line(
    aes(
      group = interaction(ID, COHORT),
      color = DRUG
    ),
    alpha = 0.35,
    linewidth = 0.6
  ) +
  
  # Individual mouse measurements
  geom_point(
    aes(color = DRUG),
    size = 2,
    alpha = 0.7
  ) +
  
  # Group mean
  geom_line(
    data = BW_plot_summary,
    aes(
      x = STATUS,
      y = mean_BW,
      group = DRUG,
      color = DRUG
    ),
    inherit.aes = FALSE,
    linewidth = 1.3
  ) +
  
  geom_point(
    data = BW_plot_summary,
    aes(
      x = STATUS,
      y = mean_BW,
      color = DRUG
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  
  # SEM
  geom_errorbar(
    data = BW_plot_summary,
    aes(
      x = STATUS,
      ymin = mean_BW - sem_BW,
      ymax = mean_BW + sem_BW,
      color = DRUG
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.8
  ) +
  
  # Rows = BW category; columns = treatment
  facet_grid(
    BW_group ~ DRUG,
    labeller = labeller(DRUG = drug_labels)
  ) +
  
  scale_color_manual(
    values = drug_colors,
    labels = drug_labels
  ) +
  
  scale_x_discrete(
    labels = c(
      "start" = "Start",
      "end" = "End"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Body weight (g)",
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

plot_BW_start_end


###4. BW change (end-start)----

BW_change <- BW_data %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, SEX, DRUG, STATUS, BW,STRAIN) %>% 
  pivot_wider(
    names_from = STATUS,
    values_from = BW
  ) %>% 
  mutate(
    BW_change = end - start,
    BW_change_percent = 100 * (end - start) / start
  )

BW_change  %>% 
  group_by(DRUG,SEX,STRAIN) %>% #so cohort 21 are just females from both strains
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf)


# Summary by treatment group 

BW_change_summary <- BW_change %>%
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    mean_change = mean(BW_change, na.rm = TRUE),
    sem_change = sd(BW_change, na.rm = TRUE) /
      sqrt(sum(!is.na(BW_change))),
    n = sum(!is.na(BW_change)),
    mean_change_percent = mean(BW_change_percent, na.rm = TRUE),
    sem_change_percent = sd(BW_change_percent, na.rm = TRUE) /
      sqrt(sum(!is.na(BW_change_percent))),
    n_percent = sum(!is.na(BW_change_percent)),
    .groups = "drop"
  )

BW_change_summary

BW_change <- BW_change %>%
  mutate(
    DRUG = factor(
      DRUG,
      levels = c("vehicle", "RTI_47")
    )
  )

plot_BW_change <- ggplot(
  BW_change_summary,
  aes(x = DRUG, y =  mean_change_percent, fill = DRUG)
) +
  geom_col(width = 0.6, color = "black") +
  
  geom_errorbar(
    aes(
      ymin =  mean_change_percent -  sem_change_percent,
      ymax =  mean_change_percent + sem_change_percent
    ),
    width = 0.2,
    linewidth = 0.7
  ) +
geom_jitter(
  data = BW_change,
  aes(x = DRUG, y = BW_change_percent),
  width = 0.10,
  size = 2.5,
  alpha = 0.7,
  color = "black",
  inherit.aes = FALSE
) +
  
  
  #geom_text(
  #  data = BW_change %>% 
  #    filter(STRAIN == "NZO/HlLtJ", SEX == "F"),
  #  aes(
  #    x = DRUG,
  #    y = BW_change,
  #    label = ID
  #  ),
  #  position = position_jitter(width = 0.10, height = 0),
  #  hjust = -0.2,
  #  size = 3
 # ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "gray60",
      "RTI_47" = "#E67E22"
    ),
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
#  coord_cartesian(ylim = c(-5, 15)) +
 # scale_y_continuous(breaks = seq(-5, 15, 5)) +
  labs(
    x = NULL,
    y = "% body weight change (100 * (end - start) / start)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~ SEX*STRAIN)

plot_BW_change

### stats of BW change----

BW_change_ttest <- BW_change %>% 
#  filter(STRAIN == "NZO/HlLtJ") %>% 
  #filter(SEX=="F") %>% 
 group_by(SEX,STRAIN) %>% 
  t_test(
    BW_change_percent ~ DRUG,
    paired = FALSE
  ) %>% 
  add_significance()

BW_change_ttest 

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

drug_colors <- c(
  "vehicle" = "gray60",
  "RTI_47" = "#E67E22"
)

drug_labels <- c(
  "vehicle" = "Vehicle",
  "RTI_47" = "RTI-47"
)


#### Fat mass plot ----

plot_fatdelta <- ggplot(
  fatdelta_summary,
  aes(
    x = DRUG,
    y = mean_delta_fat,
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
      ymin = mean_delta_fat - sem_delta_fat,
      ymax = mean_delta_fat + sem_delta_fat
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

drug_colors <- c(
  "vehicle" = "gray60",
  "RTI_47" = "#E67E22"
)

drug_labels <- c(
  "vehicle" = "Vehicle",
  "RTI_47" = "RTI-47"
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

drug_colors <- c(
  "vehicle" = "gray60",
  "RTI_47" = "#E67E22"
)

drug_labels <- c(
  "vehicle" = "Vehicle",
  "RTI_47" = "RTI-47"
)


### lean mass plot ----

plot_leandelta <- ggplot(
  leandelta_summary,
  aes(
    x = DRUG,
    y = mean_delta_lean_perc,
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
      y = delta_lean_perc
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
      ymin = mean_delta_lean_perc - sem_delta_lean_perc,
      ymax = mean_delta_lean_perc + sem_delta_lean_perc
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
    y = "Change in lean percentage (end-start)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  ) +
  
  facet_wrap(~SEX*STRAIN)

plot_leandelta

###lean mass stats----

leandelta_ttest <- delta_bodycomp %>%
  filter(!is.na(delta_lean), !is.na(DRUG)) %>%
  group_by(SEX,STRAIN) %>%
  rstatix::t_test(delta_lean_perc ~ DRUG) %>%
  rstatix::add_significance()

leandelta_ttest

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

#RTIOXA 43 for 5 days 10 mg/kg ----
# group assignation 
# first we used the adiposity index as a criteria to assign the injection groups

echoMRI_data_43_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(9, 14)) %>% 
  filter(n_measurement == 1) %>% 
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8096, 8099, 8102, 41006, 41007, 41008, 41015, 41001, 41018, 41021) ~ "RTIOXA_43", 
      ID %in% c(8095, 8100, 8103, 41003, 41005, 41012, 41010, 41016, 41014, 41020) ~ "vehicle"
    )
  ) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

plot_43_5D_assignation <- ggplot(echoMRI_data_43_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_assignation

#pairwise comparison to evaluate if the adiposity index is the same among the groups at the baseline

#to check if variances differ between groups,
echoMRI_data_43_assignation <- echoMRI_data_43_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_43_assignation) #p = 0.47 > 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_43_assignation$adiposity_index,
  g = echoMRI_data_43_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups



echoMRI_data_43 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(9, 14)) %>% 
  filter(Date %in% as.Date(c("2025-04-14", "2025-04-21", "2025-12-18","2025-12-23"))) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102, 41006, 41007,
                41008, 41015, 41001, 41018, 41021) ~ "RTIOXA_43", 
      ID %in% c(8095, 8100, 8103, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle"),
    STATUS = case_when(
      COHORT == 9 & Date == as.Date("2025-04-14") ~ "baseline",
      COHORT == 9 & Date == as.Date("2025-04-21") ~ "end_drug",
      COHORT == 14 & Date == as.Date("2025-12-18") ~ "baseline",
      COHORT == 14 & Date == as.Date("2025-12-23") ~ "end_drug")
  ) %>%
  select(ID, Fat, Lean,adiposity_index,Weight, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    delta_adiposity_index = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"],
    delta_lean = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"],
    delta_fat = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"],
    delta_bw = (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"],
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

#plot bw changes

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_bw, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% of body weight change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_bw

#plot adiposity index changes

plot_43_5D_ai <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_adiposity_index, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% adiposity index change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_ai

#plot lean mass changes

plot_43_5D_lean <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_lean, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% lean mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_lean

#plot fat mass changes

plot_43_5D_fat <- ggplot(echoMRI_data_43, aes(x = DRUG, y = delta_fat, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "% fat mass change") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_fat

# combined plot for RTIOXA-43 x 5 days

combined_plot <- (plot_43_5D_bw | plot_43_5D_ai) /
  (plot_43_5D_fat | plot_43_5D_lean) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

combined_plot

#data analysis
#1 does rtioxa-43 decreases bw after 5 days of injection?
#to check if variances differ between groups,
echoMRI_data_43 <- echoMRI_data_43 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(delta_adiposity_index ~ DRUG, data = echoMRI_data_43) #p = 0.72 > 0.05, so variances can be considered equal
leveneTest(delta_fat ~ DRUG, data = echoMRI_data_43) #p = 0.65 > 0.05, so variances can be considered equal
leveneTest(delta_lean ~ DRUG, data = echoMRI_data_43) #p = 0.14 > 0.05, so variances can be considered equal
leveneTest(delta_bw ~ DRUG, data = echoMRI_data_43) #p = 0.42 > 0.05, so variances can be considered equal

summary(aov(delta_bw ~ DRUG, data = echoMRI_data_43))  #p = 0.4
summary(aov(delta_fat ~ DRUG, data = echoMRI_data_43))  #p = 0.18
summary(aov(delta_lean ~ DRUG, data = echoMRI_data_43))  #p = 0.59
summary(aov(delta_adiposity_index ~ DRUG, data = echoMRI_data_43))  #p = 0.16

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, 
# or adiposity after 5 days. 






