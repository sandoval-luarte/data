# This script aims to explore:
#1 changes in body composition and glucose managing in C57BL6J and NZO mice after 5 weeks of IP injections of RTIOXA47 40 mg-kg
#2 changes in body composition in C57BL6J mice after 5 days of IP injections of RTIOXA43 10 mg-kg
#2 changes in body composition in C57BL6J mice after 5 days of IP injections of RTIOXA47 10 mg-kg

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

 ##4.  GLUCOSE MANAGMENT----

glu_data <- read_csv("../data/Bri_glucose_2022.csv") %>% 
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
      ) ~ "RTI_47")) %>% 
  group_by(ID) %>% 
  mutate(glu_standar = (GLUCOSE_MG_DL/BW)) %>% 
  ungroup() %>% 
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47"))
  )

glu_datasummary<- glu_data %>%
  group_by(DRUG,SEX,STRAIN) %>%
  summarise(
    mean_glu = mean(GLUCOSE_MG_DL, na.rm = TRUE),
    sem_glu  = sd(GLUCOSE_MG_DL, na.rm = TRUE) / sqrt(n()),
    mean_glu_st = mean(glu_standar, na.rm = TRUE),
    sem_glu_st  = sd(glu_standar, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### plot glucose  ----
plot_glu <- ggplot(
  glu_datasummary,
  aes(
    x = DRUG,
    y = mean_glu,
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
      ymin = mean_glu - sem_glu,
      ymax = mean_glu + sem_glu
    ),
    width = 0.2,
    linewidth = 0.7
  )  +
  
  # Individual mouse points
  geom_point(
    data = glu_data,
    aes(
      x = DRUG,
      y = GLUCOSE_MG_DL
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
  
  labs(
    x = NULL,
    y = "Glucose (mg/dL)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  )+
  facet_wrap(~SEX*STRAIN)

plot_glu

### plot glucose standardized ----
plot_glu_st <- ggplot(
  glu_datasummary,
  aes(
    x = DRUG,
    y = mean_glu_st,
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
      ymin = mean_glu_st - sem_glu_st,
      ymax = mean_glu_st + sem_glu_st
    ),
    width = 0.2,
    linewidth = 0.7
  )  +
  
  # Individual mouse points
  geom_point(
    data = glu_data,
    aes(
      x = DRUG,
      y = glu_standar
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
  
  labs(
    x = NULL,
    y = "Glucose standardized (glu/BW)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  )+
facet_wrap(~SEX*STRAIN)

plot_glu_st 

###stats----

glu_ttest <- glu_data %>%
  group_by(SEX, STRAIN) %>%
  t_test(
    glu_standar ~ DRUG,
    paired = FALSE
  ) %>%
  add_significance()

glucose_ttest <- glu_data %>%
  group_by(SEX, STRAIN) %>%
  t_test(
    GLUCOSE_MG_DL ~ DRUG,
    paired = FALSE
  ) %>%
  add_significance()

glu_ttest 
glucose_ttest

##5. PREDICTIONS ----
####Does a higher initial body weight or initial fat mass predict a greater reduction in body weight or fat mass following RTIOXA-47?----

##### 5.1 Baseline phenotype as predictor of RTIOXA-47 response ----

baseline_predictors <- echoMRI_data %>%
  filter(STATUS %in% c("start", "end")) %>%
  select(
    ID, COHORT, SEX, STRAIN, DRUG,
    STATUS, Weight, Fat, adiposity_index
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = c(Weight, Fat,adiposity_index),
    names_glue = "{.value}_{STATUS}"
  ) %>%
  mutate(
    BW_change = Weight_end - Weight_start,
    fat_change = Fat_end - Fat_start,
    adiposity_index_change = adiposity_index_end - adiposity_index_start
  )

baseline_predictors <- baseline_predictors %>%
  mutate(
    STRAIN = factor(
      STRAIN,
      levels = c("C57BL/6J", "NZO/HlLtJ")
    )
  ) 
baseline_predictors 


baseline_predictors <- baseline_predictors %>%
  mutate(
    SEX = factor(SEX, levels = c("F", "M")),
    STRAIN = factor(
      STRAIN,
      levels = c("C57BL/6J", "NZO/HlLtJ")
    ),
    DRUG = factor(
      DRUG,
      levels = c("vehicle", "RTI_47")
    )
  )


##### 5.2 Does baseline body weight predict the magnitude of the RTIOXA-47 effect on body weight change?----

BW_response_model <- lm(
  BW_change ~ DRUG * Weight_start + SEX + STRAIN,
  data = baseline_predictors
)

anova(BW_response_model)
summary(BW_response_model) #No evidence that starting weight modifies RTIOXA47 effect on BW change

#There was no statistically significant evidence that baseline body weight modified 
#the association between RTIOXA-47 treatment and body-weight change (p = 0.540). 
#Sex was associated with body-weight change, while no significant strain difference was detected. 
#The relationship between baseline body weight and body-weight change should be 
#interpreted after centering baseline weight and clarifying the main-effect test.
### BODY WEIGHT: DRUG × STRAIN x SEX ----

BW_model_strain <- lm(
  BW_change ~ Weight_start + DRUG * STRAIN * SEX,
  data = baseline_predictors
)

anova(BW_model_strain)
summary(BW_model_strain) #No evidence that RTIOXA47 efficacy differs by strain


BW_treatment_by_group <- emmeans(
  BW_model_strain,
  pairwise ~ DRUG | STRAIN * SEX,
  at = list(
    Weight_start = mean(
      baseline_predictors$Weight_start,
      na.rm = TRUE
    )
  )
)

BW_treatment_by_group #the effect of RTI-47 on body-weight change differs between female and male NZO mice

#RTI-47-treated mice had an estimated 4.00 g lower body-weight change than vehicle-treated mice (p=0.0051)
#RTI-47-treated mice had an estimated 5.50 g higher body-weight change than vehicle-treated mice (p=0.0178)

##### 5.3 Does baseline fat mass predict the magnitude of the RTIOXA-47 effect on fat mass change?----

fat_response_model <- lm(
  fat_change ~ DRUG * Fat_start + SEX + STRAIN,
  data = baseline_predictors
)

anova(fat_response_model)
summary(fat_response_model)

#There is no statistically significant evidence that RTI-47's effect on fat-mass change varies according to baseline fat mass. 
#For example, our analysis does not demonstrate that mice starting with more fat lose more fat in response to RTI-47 than mice starting with less fat

##### 5.4 Does baseline adiposity index predict the magnitude of the RTIOXA-47 effect on adiposity index change?----

adiposity_index_response_model <- lm(
  adiposity_index_change ~ DRUG * adiposity_index_start + SEX + STRAIN,
  data = baseline_predictors
)

anova(adiposity_index_response_model)
summary(adiposity_index_response_model)

#we id not find statistically significant evidence that RTI-47 changes adiposity index overall or that its effect depends on baseline adiposity index
#we  did find a significant association between sex and adiposity index change (such as for fat mass)

##6. INDIRECT CALORIMETRY----
#I will assume like as always in the columbus room was # 20:00–05:59 DARK and # 06:00–19:59 LIGHT

####counts /locomotion ----
#####start-----

ical_data_counts_coh8 <- read_csv("~/Documents/GitHub/data/data/Coh 1 basal calorimetry analysis 03042022_counts.csv")%>% 
  rename(ID = Subject) %>% 
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
  )     %>% 
  mutate(
    STRAIN = case_when(
      ID %in% c(13,14,15,16,17,18,19,20,
                21,22,23,24) ~ "NZO/HlLtJ",
      ID %in% c(1,2,3,4,5,6,7,8,9,10,
                11,12) ~ "C57BL6/J")) %>% 
  mutate(
    SEX = case_when(
      ID %in% c(1,2,3,4,5,6,13,14,15,16,17,18) ~ "M",
      ID %in% c(7,8,9,10,11,12,19,20,21,22,23,24) ~ "F"))
  

cal_counts_long <- ical_data_counts_coh8 %>%
  pivot_longer(
    cols = -c(ID, BW,DRUG,STRAIN,SEX),
    names_to = "DateTime",
    values_to = "Count"
  ) %>%
  mutate(
    DateTime = as.POSIXct(DateTime, format = "%m/%d/%y %H:%M")
  )

light_periods <- cal_counts_long %>%
  mutate(
    date = as.Date(DateTime)
  ) %>%
  distinct(date) %>%
  mutate(
    light_off_start = as.POSIXct(
      paste(date, "20:00:00")
    ),
    light_off_end = as.POSIXct(
      paste(date + 1, "06:00:00")
    )
  )

count_summary <- cal_counts_long %>%
  group_by(DateTime,DRUG,SEX,STRAIN) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    sem_count = sd(Count, na.rm = TRUE) /
      sqrt(sum(!is.na(Count))),
    .groups = "drop"
  )

ggplot() +
  
  # Light-off background
  geom_rect(
    data = light_periods,
    aes(
      xmin = light_off_start,
      xmax = light_off_end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "grey85",
    inherit.aes = FALSE
  ) +
  
  # Individual animals
  geom_line(
    data = cal_counts_long,
    aes(
      x = DateTime,
      y = Count,
      group = ID
    ),
    alpha = 0.15,
    linewidth = 0.5
  ) +
  
  geom_point(
    data = cal_counts_long,
    aes(
      x = DateTime,
      y = Count,
      group = ID
    ),
    alpha = 0.15,
    size = 1
  ) +
  
  # SEM
  geom_ribbon(
    data = count_summary,
    aes(
      x = DateTime,
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count
    ),
    alpha = 0.25
  ) +
  
  # Mean
  geom_line(
    data = count_summary,
    aes(
      x = DateTime,
      y = mean_count
    ),
    linewidth = 1.2
  ) +
  
  labs(
    x = "Date and time",
    y = "Count"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )+
  facet_wrap(~DRUG*SEX*STRAIN)

##### end-----

ical_data_counts_coh8_end <- read_csv("~/Documents/GitHub/data/data/Coh 1 final calorimetry analysis 04082022_counts.csv")%>% 
  rename(ID = Subject) %>% 
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
  )     %>% 
  mutate(
    STRAIN = case_when(
      ID %in% c(13,14,15,16,17,18,19,20,
                21,22,23,24) ~ "NZO/HlLtJ",
      ID %in% c(1,2,3,4,5,6,7,8,9,10,
                11,12) ~ "C57BL6/J")) %>% 
  mutate(
    SEX = case_when(
      ID %in% c(1,2,3,4,5,6,13,14,15,16,17,18) ~ "M",
      ID %in% c(7,8,9,10,11,12,19,20,21,22,23,24) ~ "F"))

cal_counts_long_end <- ical_data_counts_coh8_end %>%
  pivot_longer(
    cols = -c(ID, BW,DRUG,SEX,STRAIN),
    names_to = "DateTime",
    values_to = "Count"
  ) %>%
  mutate(
    DateTime = as.POSIXct(DateTime, format = "%m/%d/%y %H:%M")
  )

light_periods_end <- cal_counts_long_end %>%
  mutate(
    date = as.Date(DateTime)
  ) %>%
  distinct(date) %>%
  mutate(
    light_off_start = as.POSIXct(
      paste(date, "20:00:00")
    ),
    light_off_end = as.POSIXct(
      paste(date + 1, "06:00:00")
    )
  )

count_summary_end <- cal_counts_long_end %>%
  group_by(DateTime,DRUG,SEX,STRAIN) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    sem_count = sd(Count, na.rm = TRUE) /
      sqrt(sum(!is.na(Count))),
    .groups = "drop"
  )

ggplot() +
  
  # Light-off background
  geom_rect(
    data = light_periods_end,
    aes(
      xmin = light_off_start,
      xmax = light_off_end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "grey85",
    inherit.aes = FALSE
  ) +
  
  # Individual animals
  geom_line(
    data = cal_counts_long_end,
    aes(
      x = DateTime,
      y = Count,
      group = ID
    ),
    alpha = 0.15,
    linewidth = 0.5
  ) +
  
  geom_point(
    data = cal_counts_long_end,
    aes(
      x = DateTime,
      y = Count,
      group = ID
    ),
    alpha = 0.15,
    size = 1
  ) +
  
  # SEM
  geom_ribbon(
    data = count_summary_end,
    aes(
      x = DateTime,
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count
    ),
    alpha = 0.25
  ) +
  
  # Mean
  geom_line(
    data = count_summary_end,
    aes(
      x = DateTime,
      y = mean_count
    ),
    linewidth = 1.2
  ) +
  
  labs(
    x = "Date and time",
    y = "Count at end"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )+
  facet_wrap(~DRUG*SEX*STRAIN)

#start----

cal_counts_long <- ical_data_counts_coh8 %>%
  pivot_longer(
    cols = -c(ID, BW, DRUG,SEX,STRAIN),
    names_to = "DateTime",
    values_to = "Count"
  ) %>%
  mutate(
    DateTime = as.POSIXct(
      DateTime,
      format = "%m/%d/%y %H:%M"
    ),
    
    # Experimental start
    experiment_start = as.POSIXct(
      "2022-02-28 20:00:00"
    ),
    
    # Assign experimental day
    day = if_else(
      DateTime < experiment_start,
      NA_integer_,
      floor(
        as.numeric(
          difftime(
            DateTime,
            experiment_start,
            units = "hours"
          )
        ) / 24
      ) + 1
    )
  )

cal_counts_long_complete <- cal_counts_long %>%
  filter(
    !is.na(day),
    day != 4
  )
cal_counts_long_complete <- cal_counts_long_complete %>%
  mutate(
    DATE = as.Date(DateTime),
    TIME = format(DateTime, "%H:%M:%S")
  )
cal_counts_long_complete <- cal_counts_long_complete %>%
  mutate(
    lights = case_when(
      TIME >= "20:00:00" | TIME < "06:00:00" ~ "OFF",
      TIME >= "06:00:00" & TIME < "20:00:00" ~ "ON",
      TRUE ~ NA_character_
    )
  )

#### Mean and SD per ID, day, and hour ####

count_hour_summary <- cal_counts_long_complete %>%
  mutate(
    HOUR = as.integer(format(DateTime, "%H"))
  ) %>%
  group_by(ID, DRUG, HOUR, lights,SEX,STRAIN) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    sd_count = sd(Count, na.rm = TRUE),
    n = sum(!is.na(Count)),
    .groups = "drop"
  )

#### Plot hourly locomotor activity averaged across 3 days ####

#### Reorder hours: 20:00 to 19:00 ####

count_hour_summary <- count_hour_summary %>%
  mutate(
    HOUR_ORDER = if_else(HOUR >= 20, HOUR - 20, HOUR + 4)
  )

count_group_summary <- count_hour_summary %>%
  group_by(DRUG, HOUR_ORDER, lights,SEX,STRAIN) %>%
  summarise(
    mean_count = mean(mean_count, na.rm = TRUE),
    sem_count = sd(mean_count, na.rm = TRUE) /
      sqrt(sum(!is.na(mean_count))),
    .groups = "drop"
  )

ggplot() +
  
  # Individual mice
  geom_line(
    data = count_hour_summary,
    aes(
      x = HOUR_ORDER,
      y = mean_count,
      group = ID
    ),
    alpha = 0.5,
    linewidth = 0.6,
    color = "grey50"
  ) +
  
  # Mean ± SEM
  geom_ribbon(
    data = count_group_summary,
    aes(
      x = HOUR_ORDER,
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count,
      group = DRUG
    ),
    alpha = 0.2
  ) +
  
  # Group mean
  geom_line(
    data = count_group_summary,
    aes(
      x = HOUR_ORDER,
      y = mean_count,
      color = DRUG,
      group = DRUG
    ),
    linewidth = 1.5
  ) +
  
  scale_x_continuous(
    breaks = 0:23,
    labels = c(20:23, 0:19)
  ) +
  
  scale_color_manual(
    values = c(
      "vehicle" = "grey40",
      "RTI_47" = "#E67E22"
    )
  ) +
  
  labs(
    x = "Hour of experimental day",
    y = "Locomotor counts",
    color = "Treatment"
  ) +
  
  theme_classic() +
  
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  facet_wrap(~DRUG*SEX*STRAIN)

#### Cumulative locomotor counts per ID, day, and lights ####

count_cumulative_lights <- cal_counts_long_complete %>%
  group_by(ID, DRUG, day, lights,SEX,STRAIN) %>%
  summarise(
    cumulative_count = sum(Count, na.rm = TRUE),
    n_measurements = sum(!is.na(Count)),
    .groups = "drop"
  ) %>%
  arrange(ID, day, lights)

#### Mean ± SEM of cumulative counts per day and light period ####

count_cumulative_summary <- count_cumulative_lights %>%
  group_by(DRUG, day, lights,SEX,STRAIN) %>%
  summarise(
    mean_count = mean(cumulative_count, na.rm = TRUE),
    sem_count = sd(cumulative_count, na.rm = TRUE) /
      sqrt(sum(!is.na(cumulative_count))),
    n = sum(!is.na(cumulative_count)),
    .groups = "drop"
  )

#### Bar plot: cumulative counts per light period ####

ggplot(
  count_cumulative_summary,
  aes(x = lights, y = mean_count, fill = DRUG)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.8
  ) +
  
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count
    ),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "grey60",
      "RTI_47" = "#E67E22"
    )
  ) +
  
  labs(
    x = "Light period",
    y = "Cumulative locomotor counts (mean ± SEM)",
    fill = "Treatment"
  ) +
  
  facet_wrap(~day*SEX*STRAIN, nrow = 1) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

ggplot(
  count_cumulative_summary,
  aes(x = lights, y = mean_count, fill = DRUG)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.8
  ) +
  
  # Individual mice
  geom_point(
    data = count_cumulative_lights,
    aes(x = lights, y = cumulative_count, color = DRUG),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count
    ),
    position = position_dodge(width = 0.8),
    width = 0.2,
    color = "black",
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      "vehicle" = "grey60",
      "RTI_47" = "#E67E22"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "vehicle" = "grey20",
      "RTI_47" = "#A84F00"
    ),
    guide = "none"
  ) +
  
  labs(
    x = "Light period",
    y = "Cumulative locomotor counts (mean ± SEM)",
    fill = "Treatment"
  ) +
  
  facet_wrap(~day*SEX*STRAIN, nrow = 1) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

ggplot(
  count_cumulative_summary,
  aes(x = lights, y = mean_count, fill = DRUG)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.8
  ) +
  
  # Individual mice
  geom_point(
    data = count_cumulative_lights,
    aes(x = lights, y = cumulative_count, color = DRUG),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  # Mouse IDs
  geom_text(
    data = count_cumulative_lights,
    aes(
      x = lights,
      y = cumulative_count,
      label = ID,
      color = DRUG
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.8
    ),
    size = 2.5,
    vjust = -0.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = mean_count - sem_count,
      ymax = mean_count + sem_count
    ),
    position = position_dodge(width = 0.8),
    width = 0.2,
    color = "black",
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c("vehicle" = "grey60", "RTI_47" = "#E67E22")
  ) +
  
  scale_color_manual(
    values = c("vehicle" = "grey20", "RTI_47" = "#A84F00"),
    guide = "none"
  ) +
  
  labs(
    x = "Light period",
    y = "Cumulative locomotor counts (mean ± SEM)",
    fill = "Treatment"
  ) +
  
  facet_wrap(~day*SEX*STRAIN, nrow = 1) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 11, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )



# Unpaired Welch's t-test: vehicle vs RTI-47
t_test_counts <- count_cumulative_lights %>%
  group_by(STRAIN, SEX, day, lights) %>%
  t_test(
    cumulative_count ~ DRUG,
    paired = FALSE,
    var.equal = FALSE
  ) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  ungroup()

# View results
t_test_counts %>%
  select(STRAIN, SEX, day, lights, group1, group2,
         n1, n2, statistic, df, p, p.adj, p.adj.signif) 

  print(t_test_counts, n = Inf)


  
  # Prepare the data
  cal_model_data <- count_cumulative_lights %>%
    mutate(
      ID = factor(ID),
      DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
      STRAIN = factor(STRAIN),
      SEX = factor(SEX),
      lights = factor(lights, levels = c("OFF", "ON")),
      day = factor(day),
      log_count = log1p(cumulative_count)
    )
  
  # Check that each mouse has one observation per day and light period
  cal_model_data %>%
    count(ID, day, lights) %>%
    filter(n > 1)
  
  # Fit the mixed-effects model
  model_counts <- lmer(
    log_count ~ DRUG * (STRAIN + SEX + lights) +
      day + (1 | ID),
    data = cal_model_data,
    REML = TRUE
  )
  
  # View the model results
  summary(model_counts)
  
  # Test the fixed effects
  anova(model_counts, type = 3)
  

  # Average over strain, sex, light period, and day
  drug_emm <- emmeans(
    model_counts,
    ~ DRUG,
    weights = "equal"
  )
  
  # Compare RTI-47 with vehicle
  pairs(drug_emm, adjust = "none")
  
  # Display estimated means and confidence intervals
  summary(drug_emm, infer = TRUE)
  
#Drug effect: No statistically significant overall difference was detected.
#Strain: Locomotor activity differed significantly between C57BL/6J and NZO/HlLtJ mice.
#Sex: Locomotor activity differed significantly between females and males.
#Light period: Activity differed significantly between ON and OFF periods.
#Drug interactions: No significant interactions with strain, sex, or light period were detected.

  
#end----

  
  # Start with the final calorimetry dataset
  cal_counts_long_end <- ical_data_counts_coh8_end %>%
    pivot_longer(
      cols = -c(ID, BW, DRUG, SEX, STRAIN),
      names_to = "DateTime",
      values_to = "Count"
    ) %>%
    mutate(
      DateTime = as.POSIXct(
        DateTime,
        format = "%m/%d/%y %H:%M"
      ),
      experiment_start = as.POSIXct("2022-04-05 10:00:00"),
      day = floor(
        as.numeric(
          difftime(DateTime, experiment_start, units = "hours")
        ) / 24
      ) + 1
    )
  
  # Check the time range and groups
  range(cal_counts_long_end$DateTime, na.rm = TRUE)
  
  ical_data_counts_coh8_end %>%
    count(DRUG, STRAIN, SEX)
  
  cal_counts_long_end <- cal_counts_long_end %>%
    mutate(
      TIME = format(DateTime, "%H:%M:%S"),
      lights = case_when(
        TIME >= "20:00:00" | TIME < "06:00:00" ~ "OFF",
        TIME >= "06:00:00" & TIME < "20:00:00" ~ "ON",
        TRUE ~ NA_character_
      )
    )
  
  # Keep only complete 24-hour experimental days
  count_cumulative_lights_end <- cal_counts_long_end %>%
    filter(day %in% 1:3, !is.na(lights)) %>%
    group_by(ID, DRUG, STRAIN, SEX, day, lights) %>%
    summarise(
      cumulative_count = sum(Count, na.rm = TRUE),
      n_measurements = sum(!is.na(Count)),
      .groups = "drop"
    )
  
  # Check observations per mouse, day, and light period
  count_cumulative_lights_end %>%
    count(ID, day, lights) %>%
    filter(n != 1)
  
  # Check recording completeness
  count_cumulative_lights_end %>%
    group_by(day, lights) %>%
    summarise(
      min_measurements = min(n_measurements),
      max_measurements = max(n_measurements),
      .groups = "drop"
    )
  
  cal_counts_long_end %>%
    distinct(ID, DRUG, STRAIN, SEX) %>%
    arrange(STRAIN, SEX, DRUG)
  
  cal_model_data_end <- count_cumulative_lights_end %>%
    mutate(
      ID = factor(ID),
      DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
      SEX = factor(SEX),
      lights = factor(lights, levels = c("OFF", "ON")),
      day = factor(day),
      log_count = log1p(cumulative_count)
    )
  
  model_counts_end <- lmer(
    log_count ~ DRUG * (SEX + lights) +
      day + (1 | ID),
    data = cal_model_data_end,
    REML = TRUE
  )
  
  summary(model_counts_end)
  
  anova(model_counts_end, type = 3)
  
  # Estimated treatment difference averaged across sex,
  # light period, and day
  drug_emm_end <- emmeans(
    model_counts_end,
    ~ DRUG,
    weights = "equal"
  )
  
  pairs(drug_emm_end, adjust = "none")
  summary(drug_emm_end, infer = TRUE)


   ggplot(
    count_cumulative_lights_end,
    aes(x = DRUG, y = cumulative_count, fill = DRUG)
  ) +
    geom_boxplot(
      alpha = 0.25,
      width = 0.55,
      outlier.shape = NA
    ) +
    geom_point(
      aes(color = DRUG),
      position = position_jitter(width = 0.08),
      size = 2.5,
      alpha = 0.8
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 3.5,
      color = "black"
    ) 
   
   
   ggplot(
     count_cumulative_lights_end,
     aes(x = DRUG, y = cumulative_count, fill = DRUG)
   ) +
     geom_boxplot(
       alpha = 0.20,
       width = 0.5,
       outlier.shape = NA
     ) +
     geom_point(
       aes(color = DRUG),
       position = position_jitter(width = 0.08, height = 0),
       size = 2.5,
       alpha = 0.8
     ) +
     geom_text_repel(
       aes(label = ID),
       size = 3,
       max.overlaps = Inf,
       show.legend = FALSE
     ) +
     stat_summary(
       fun = mean,
       geom = "point",
       shape = 18,
       size = 3.5,
       color = "black"
     ) +
     stat_summary(
       fun.data = mean_se,
       geom = "errorbar",
       width = 0.15,
       color = "black"
     ) +
     facet_grid(SEX ~ lights + day) +
     scale_fill_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     scale_color_manual(
       values = c("vehicle" = "gray40", "RTI_47" = "#E67E22")
     ) +
     labs(
       x = "Treatment",
       y = "Cumulative locomotor counts",
       title = "End calorimetry: individual locomotor activity",
       subtitle = "Individual mice labeled by ID; mean ± SEM"
     ) +
     theme_classic(base_family = "Arial") +
     theme(legend.position = "none")

   t_test_counts_end <- count_cumulative_lights_end %>%
     group_by(SEX, lights, day) %>%
     t_test(
       cumulative_count ~ DRUG,
       paired = FALSE,
       var.equal = FALSE
     ) %>%
     adjust_pvalue(method = "holm") %>%
     add_significance("p.adj") %>%
     ungroup()
   
   t_test_counts_end
   
   #These analyses did not detect a significant difference
   #in cumulative locomotor counts between vehicle- and RTI-47-treated 
   #NZO/HlLtJ mice in any sex, light period, or recording day.
   
   
# 4hours after injections----
   
   counts_11am_2pm_end <- cal_counts_long_end %>%
     filter(
       day %in% 1:3,
       format(DateTime, "%H:%M:%S") >= "11:00:00",
       format(DateTime, "%H:%M:%S") <= "14:00:00"
     ) %>%
     group_by(ID, DRUG, STRAIN, SEX, day) %>%
     summarise(
       count_11am_2pm = sum(Count, na.rm = TRUE),
       n_measurements = sum(!is.na(Count)),
       .groups = "drop"
     ) %>%
     arrange(ID, day)
   
   counts_11am_2pm_end
   
   
   ggplot(
     counts_11am_2pm_end,
     aes(x = DRUG, y = count_11am_2pm, fill = DRUG)
   ) +
     geom_boxplot(
       alpha = 0.20,
       width = 0.5,
       outlier.shape = NA
     ) +
     geom_point(
       aes(color = DRUG),
       position = position_jitter(width = 0.08, height = 0),
       size = 2.5,
       alpha = 0.8
     ) +
     geom_text_repel(
       aes(label = ID),
       size = 3,
       max.overlaps = Inf,
       show.legend = FALSE
     ) +
     stat_summary(
       fun = mean,
       geom = "point",
       shape = 18,
       size = 3.5,
       color = "black"
     ) +
     stat_summary(
       fun.data = mean_se,
       geom = "errorbar",
       width = 0.15,
       color = "black"
     ) +
     facet_grid(SEX ~ day) +
     scale_fill_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     scale_color_manual(
       values = c("vehicle" = "gray40", "RTI_47" = "#E67E22")
     ) +
     labs(
       x = "Treatment",
       y = "Cumulative locomotor counts (11 AM–2 PM)",
       title = "End calorimetry: locomotor activity",
       subtitle = "Individual mice; mean ± SEM"
     ) +
     theme_classic(base_family = "Arial") +
     theme(legend.position = "none")

   # Unpaired t-tests: 11 AM–2 PM
   t_test_11am_2pm_end <- counts_11am_2pm_end %>%
     group_by(SEX, day) %>%
     t_test(
       count_11am_2pm ~ DRUG,
       paired = FALSE,
       var.equal = FALSE
     ) %>%
     adjust_pvalue(method = "holm") %>%
     add_significance("p.adj") %>%
     ungroup()
   
   # Display results
   t_test_11am_2pm_end %>%
     select(
       SEX, day, group1, group2,
       n1, n2, statistic, df, p, p.adj, p.adj.signif
     )
  
  #heat also called tee----
   
   # End calorimetry: heat ----
   ical_data_heat_coh8_end <- read_csv(
     "~/Documents/GitHub/data/data/Coh 1 final calorimetry analysis 04082022_heat.csv"
   ) %>%
     rename(ID = Subject) %>%
     mutate(
       DRUG = case_when(
         ID %in% c(13, 15, 17, 19, 21, 23) ~ "vehicle",
         ID %in% c(14, 16, 18, 20, 22, 24) ~ "RTI_47",
         TRUE ~ NA_character_
       ),
       STRAIN = case_when(
         ID %in% 13:24 ~ "NZO/HlLtJ",
         TRUE ~ NA_character_
       ),
       SEX = case_when(
         ID %in% 13:18 ~ "M",
         ID %in% 19:24 ~ "F",
         TRUE ~ NA_character_
       )
     ) %>%
     mutate(
       DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
       STRAIN = factor(STRAIN),
       SEX = factor(SEX, levels = c("F", "M"))
     )
   
   # Verify animal assignments
   ical_data_heat_coh8_end %>%
     distinct(ID, DRUG, STRAIN, SEX) %>%
     arrange(SEX, DRUG, ID)
   
   ical_data_heat_coh8_end %>%
     count(DRUG, STRAIN, SEX)
   
   cal_heat_long_end <- ical_data_heat_coh8_end %>%
     pivot_longer(
       cols = -c(ID, BW, DRUG, STRAIN, SEX),
       names_to = "DateTime",
       values_to = "Heat"
     ) %>%
     mutate(
       DateTime = as.POSIXct(
         DateTime,
         format = "%m/%d/%y %H:%M"
       ),
       experiment_start = as.POSIXct("2022-04-05 10:00:00"),
       day = floor(
         as.numeric(
           difftime(DateTime, experiment_start, units = "hours")
         ) / 24
       ) + 1,
       TIME = format(DateTime, "%H:%M:%S"),
       lights = case_when(
         TIME >= "21:00:00" | TIME < "06:00:00" ~ "OFF",
         TIME >= "07:00:00" & TIME < "19:00:00" ~ "ON",
         TRUE ~ NA_character_
       )
     )
   
   # Check the recording
   range(cal_heat_long_end$DateTime, na.rm = TRUE)
   
   cal_heat_long_end %>%
     count(ID, day) %>%
     filter(day %in% 1:3)
   
   heat_cumulative_lights_end <- cal_heat_long_end %>%
     filter(day %in% 1:3, !is.na(lights)) %>%
     group_by(ID, DRUG, STRAIN, SEX, day, lights) %>%
     summarise(
       cumulative_heat = sum(Heat, na.rm = TRUE),
       n_measurements = sum(!is.na(Heat)),
       .groups = "drop"
     ) %>%
     arrange(ID, day, lights)
   
   # Check recording completeness
   heat_cumulative_lights_end %>%
     group_by(day, lights) %>%
     summarise(
       min_measurements = min(n_measurements),
       max_measurements = max(n_measurements),
       .groups = "drop"
     )
   
   ggplot(
     heat_cumulative_lights_end,
     aes(x = DRUG, y = cumulative_heat, fill = DRUG)
   ) +
     geom_boxplot(
       alpha = 0.20,
       width = 0.5,
       outlier.shape = NA
     ) +
     geom_point(
       aes(color = DRUG),
       position = position_jitter(width = 0.08, height = 0),
       size = 2.5,
       alpha = 0.8
     ) +
     geom_text_repel(
       aes(label = ID),
       size = 3,
       max.overlaps = Inf,
       show.legend = FALSE
     ) +
     stat_summary(
       fun = mean,
       geom = "point",
       shape = 18,
       size = 3.5,
       color = "black"
     ) +
     stat_summary(
       fun.data = mean_se,
       geom = "errorbar",
       width = 0.15,
       color = "black"
     ) +
     facet_grid(SEX ~ lights + day) +
     scale_fill_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     scale_color_manual(
       values = c("vehicle" = "gray40", "RTI_47" = "#E67E22")
     ) +
     labs(
       x = "Treatment",
       y = "Cumulative heat (sum of hourly values)",
       title = "End calorimetry: heat production",
       subtitle = "Individual mice; mean ± SEM"
     ) +
     theme_classic(base_family = "Arial") +
     theme(legend.position = "none")
   
   t_test_heat_lights_end <- heat_cumulative_lights_end %>%
     group_by(SEX, lights, day) %>%
     t_test(
       cumulative_heat ~ DRUG,
       paired = FALSE,
       var.equal = FALSE
     ) %>%
     adjust_pvalue(method = "holm") %>%
     add_significance("p.adj") %>%
     ungroup()
   
   t_test_heat_lights_end %>%
     select(
       SEX, lights, day,
       group1, group2, n1, n2,
       statistic, df, p, p.adj, p.adj.signif
     )
   
   heat_model_data_end <- heat_cumulative_lights_end %>%
     mutate(
       ID = factor(ID),
       DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
       SEX = factor(SEX),
       lights = factor(lights, levels = c("OFF", "ON")),
       day = factor(day),
       log_heat = log1p(cumulative_heat)
     )
   
   model_heat_end <- lmer(
     log_heat ~ DRUG * (SEX + lights) +
       day + (1 | ID),
     data = heat_model_data_end,
     REML = TRUE
   )
   
   summary(model_heat_end)
   anova(model_heat_end, type = 3)
   
   # Overall treatment comparison, averaged across sex, light period, and day
   heat_emm_end <- emmeans(
     model_heat_end,
     ~ DRUG,
     weights = "equal"
   )
   
   pairs(heat_emm_end, adjust = "none")
   summary(heat_emm_end, infer = TRUE)
   
#4 hours after injection heat----
   
   heat_11am_2pm_end <- cal_heat_long_end %>%
     filter(
       day %in% 1:3,
       format(DateTime, "%H:%M:%S") >= "11:00:00",
       format(DateTime, "%H:%M:%S") <= "14:00:00"
     ) %>%
     group_by(ID, DRUG, STRAIN, SEX, day) %>%
     summarise(
       heat_11am_2pm = sum(Heat, na.rm = TRUE),
       n_measurements = sum(!is.na(Heat)),
       .groups = "drop"
     ) %>%
     arrange(ID, day)
   
   # Check that four hourly values are included
   heat_11am_2pm_end %>%
     filter(n_measurements != 4)
   
   ggplot(
     heat_11am_2pm_end,
     aes(x = DRUG, y = heat_11am_2pm, fill = DRUG)
   ) +
     geom_boxplot(
       alpha = 0.20,
       width = 0.5,
       outlier.shape = NA
     ) +
     geom_point(
       aes(color = DRUG),
       position = position_jitter(width = 0.08, height = 0),
       size = 2.5,
       alpha = 0.8
     ) +
     geom_text_repel(
       aes(label = ID),
       size = 3,
       max.overlaps = Inf,
       show.legend = FALSE
     ) +
     stat_summary(
       fun = mean,
       geom = "point",
       shape = 18,
       size = 3.5,
       color = "black"
     ) +
     stat_summary(
       fun.data = mean_se,
       geom = "errorbar",
       width = 0.15,
       color = "black"
     ) +
     facet_grid(SEX ~ day) +
     scale_fill_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     scale_color_manual(
       values = c("vehicle" = "gray40", "RTI_47" = "#E67E22")
     ) +
     labs(
       x = "Treatment",
       y = "Heat (sum of hourly values, 11 AM–2 PM)",
       title = "End calorimetry: heat during the post-injection window",
       subtitle = "Individual mice; mean ± SEM"
     ) +
     theme_classic(base_family = "Arial") +
     theme(legend.position = "none")
   
   t_test_heat_11am_2pm_end <- heat_11am_2pm_end %>%
     group_by(SEX, day) %>%
     t_test(
       heat_11am_2pm ~ DRUG,
       paired = FALSE,
       var.equal = FALSE
     ) %>%
     adjust_pvalue(method = "holm") %>%
     add_significance("p.adj") %>%
     ungroup()
   
   t_test_heat_11am_2pm_end %>%
     select(
       SEX, day, group1, group2,
       n1, n2, statistic, df, p, p.adj, p.adj.signif
     )
   
#heat adjusted by BW----
   # Add body weight to the model dataset
   heat_model_data_end <- heat_cumulative_lights_end %>%
     left_join(
       ical_data_heat_coh8_end %>%
         select(ID, BW) %>%
         distinct(),
       by = "ID"
     ) %>%
     mutate(
       ID = factor(ID),
       DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
       SEX = factor(SEX),
       lights = factor(lights, levels = c("OFF", "ON")),
       day = factor(day),
       log_heat = log1p(cumulative_heat)
     )
   
   # Check that BW was added correctly
   heat_model_data_end %>%
     select(ID, BW, DRUG, SEX, day, lights, cumulative_heat) %>%
     head(12)
   
   # Check for missing body weights
   sum(is.na(heat_model_data_end$BW))

   model_heat_BW_end <- lmer(
     log_heat ~ DRUG * (SEX + lights) +
       day + BW + (1 | ID),
     data = heat_model_data_end,
     REML = TRUE
   )
   
   summary(model_heat_BW_end)
   anova(model_heat_BW_end, type = 3)
   
 ##In NZO/HlLtJ mice, absolute heat values were significantly higher during the light phase than the dark phase. 
#This effect remained significant after accounting for body weight. 
#No significant effect of RTI-47 treatment or treatment interactions with sex or light period was detected.
 
   # Center body weight around the overall mean
   heat_model_data_end <- heat_model_data_end %>%
     mutate(
       BW_centered = BW - mean(BW, na.rm = TRUE)
     )
   
   # Refit the model with centered BW
   model_heat_BW_end <- lmer(
     log_heat ~ DRUG * (SEX + lights) +
       day + BW_centered + (1 | ID),
     data = heat_model_data_end,
     REML = TRUE
   )
   
   # Adjust each observation to the mean BW using the estimated BW coefficient
   beta_BW <- fixef(model_heat_BW_end)["BW_centered"]
   
   heat_BW_adjusted <- heat_model_data_end %>%
     mutate(
       log_heat_adjusted = log_heat - beta_BW * BW_centered,
       heat_adjusted = expm1(log_heat_adjusted)
     ) %>%
     group_by(ID, DRUG, SEX, lights, day) %>%
     summarise(
       heat_adjusted = mean(heat_adjusted, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     group_by(ID, DRUG, SEX, lights) %>%
     summarise(
       heat_adjusted = mean(heat_adjusted, na.rm = TRUE),
       .groups = "drop"
     )
   
   # Summary: mean and SEM across individual mice
   heat_BW_summary <- heat_BW_adjusted %>%
     group_by(DRUG, SEX, lights) %>%
     summarise(
       mean_heat = mean(heat_adjusted, na.rm = TRUE),
       SEM = sd(heat_adjusted, na.rm = TRUE) / sqrt(n()),
       n = n(),
       .groups = "drop"
     )
   
   # Plot
   ggplot(
     heat_BW_adjusted,
     aes(x = DRUG, y = heat_adjusted, color = DRUG)
   ) +
     geom_boxplot(
       aes(fill = DRUG),
       alpha = 0.15,
       width = 0.55,
       outlier.shape = NA
     ) +
     geom_point(
       position = position_jitter(width = 0.07, height = 0),
       size = 2.5
     ) +
     geom_text_repel(
       aes(label = ID),
       size = 3,
       show.legend = FALSE,
       max.overlaps = Inf
     ) +
     geom_errorbar(
       data = heat_BW_summary,
       aes(
         x = DRUG,
         y = mean_heat,
         ymin = mean_heat - SEM,
         ymax = mean_heat + SEM
       ),
       inherit.aes = FALSE,
       color = "black",
       width = 0.15,
       linewidth = 0.7
     ) +
     geom_point(
       data = heat_BW_summary,
       aes(x = DRUG, y = mean_heat),
       inherit.aes = FALSE,
       color = "black",
       shape = 18,
       size = 3.5
     ) +
     facet_grid(SEX ~ lights) +
     scale_color_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     scale_fill_manual(
       values = c("vehicle" = "gray60", "RTI_47" = "#E67E22")
     ) +
     labs(
       x = NULL,
       y = "BW-adjusted cumulative heat",
       color = "Treatment",
       fill = "Treatment"
     ) +
     theme_classic(base_family = "Arial") +
     theme(
       strip.background = element_blank(),
       strip.text = element_text(size = 12),
       axis.text = element_text(size = 11),
       axis.title.y = element_text(size = 12),
       legend.position = "top"
     )
   
   
   
   #RTIOXA 43 for 5 days 10 mg/kg ----
##### group assignation ---------
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
    adiposity_index_change = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"]),
    lean_change = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"]),
    fat_change = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"]),
    BW_change =  (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"]),
    adiposity_index_change_percent = (100 * (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"]),
    lean_change_percent = (100 * (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"]),
    fat_change_percent =(100 * (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"]),
    BW_change_percent = (100 * (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"]),
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43"))) %>% 
  drop_na()

###1. BW change plot----

plot_43_5D_bw <- ggplot(echoMRI_data_43, aes(x = DRUG, y = BW_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "BW change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_bw

###2.  adiposity index change plot----

plot_43_5D_ai <- ggplot(echoMRI_data_43, aes(x = DRUG, y = adiposity_index_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_ai

###3. lean mass change plot----

plot_43_5D_lean <- ggplot(echoMRI_data_43, aes(x = DRUG, y = lean_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "lean mass change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_lean

###4. fat mass change plot----

plot_43_5D_fat <- ggplot(echoMRI_data_43, aes(x = DRUG, y =fat_change , fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "fat mass change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_43" = "#66C2A5")) +
  theme(legend.position = "none")

plot_43_5D_fat

###5. stats ----
#1 does rtioxa-43 decreases bw after 5 days of injection?
#to check if variances differ between groups,
echoMRI_data_43 <- echoMRI_data_43 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(adiposity_index_change ~ DRUG, data = echoMRI_data_43) #p = 0.8 > 0.05, so variances can be considered equal
leveneTest(fat_change ~ DRUG, data = echoMRI_data_43) #p = 0.8 > 0.05, so variances can be considered equal
leveneTest(lean_change ~ DRUG, data = echoMRI_data_43) #p = 0.2 > 0.05, so variances can be considered equal
leveneTest(BW_change ~ DRUG, data = echoMRI_data_43) #p = 0.3 > 0.05, so variances can be considered equal
leveneTest(adiposity_index_change_percent ~ DRUG, data = echoMRI_data_43) #p = 0.7 > 0.05, so variances can be considered equal
leveneTest(fat_change_percent ~ DRUG, data = echoMRI_data_43) #p = 0.6 > 0.05, so variances can be considered equal
leveneTest(lean_change_percent ~ DRUG, data = echoMRI_data_43) #p = 0.1 > 0.05, so variances can be considered equal
leveneTest(BW_change_percent ~ DRUG, data = echoMRI_data_43) #p = 0.42 > 0.05, so variances can be considered equal

summary(aov(BW_change ~ DRUG, data = echoMRI_data_43))  #p = 0.4
summary(aov(fat_change ~ DRUG, data = echoMRI_data_43))  #p = 0.4
summary(aov(lean_change ~ DRUG, data = echoMRI_data_43))  #p = 0.5
summary(aov(adiposity_index_change ~ DRUG, data = echoMRI_data_43))  #p = 0.3
summary(aov(BW_change_percent ~ DRUG, data = echoMRI_data_43))  #p = 0.4
summary(aov(fat_change_percent ~ DRUG, data = echoMRI_data_43))  #p = 0.1
summary(aov(lean_change_percent ~ DRUG, data = echoMRI_data_43))  #p = 0.5
summary(aov(adiposity_index_change_percent ~ DRUG, data = echoMRI_data_43))  #p = 0.1

#conclusion> There is no indication that RTIOXA-43 decreases body weight, fat mass, lean mass or adiposity after 5 days. 


#RTIOXA 47 for 5 days 10 mg/kg ----

#WE SHOULD CONFIRM IF THESE ANIMALS WERE ALSO EXPOSED TO RTI43.
##### group assignation ---------
echoMRI_data_47_assignation <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(12, 14)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    DRUG_ASSIGNATION = case_when(
      ID %in% c(8075, 8077, 8078, 41002, 41004, 
                41009, 41019, 41017, 41013, 41011) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle"
    )
  ) %>%
  filter(n_measurement == 1) %>% 
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION, levels = c("vehicle", "RTIOXA_47"))) %>% 
  drop_na()

plot_47_5D_assignation <- ggplot(echoMRI_data_47_assignation, aes(x = DRUG_ASSIGNATION, y = adiposity_index, fill = DRUG_ASSIGNATION)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index at baseline") + 
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange"))+
  theme(legend.position = "none") 
plot_47_5D_assignation

#pairwise comparison to evaluate if the adiposity index is the same between groups at the baseline

#to check if variances differ between groups,
echoMRI_data_47_assignation <- echoMRI_data_47_assignation %>%
  mutate(DRUG_ASSIGNATION = factor(DRUG_ASSIGNATION))

leveneTest(adiposity_index ~ DRUG_ASSIGNATION, data = echoMRI_data_47_assignation) #p = 0.95> 0.05, so variances can be considered equal

pairwise.t.test(
  x = echoMRI_data_47_assignation$adiposity_index,
  g = echoMRI_data_47_assignation$DRUG_ASSIGNATION,
  p.adjust.method = "bonferroni",  
  pool.sd = TRUE                     
) #No statistically significant differences in adiposity_index between groups


echoMRI_data_47 <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT %in% c(12, 14)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078, 41002, 41004, 
                41009, 41019, 41017, 41013, 41011) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle"),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      n_measurement == 2 ~ "end_drug")
  ) %>%
  select(ID, Fat, Lean,adiposity_index,Weight, DRUG, STATUS) %>%
  group_by(ID, DRUG) %>%
  summarise(
    adiposity_index_change = (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"]),
    lean_change = (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"]),
    fat_change = (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"]),
    BW_change =  (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"]),
    adiposity_index_change_percent = (100 * (adiposity_index[STATUS == "end_drug"] - adiposity_index[STATUS == "baseline"])/adiposity_index[STATUS == "baseline"]),
    lean_change_percent = (100 * (Lean[STATUS == "end_drug"] - Lean[STATUS == "baseline"])/Lean[STATUS == "baseline"]),
    fat_change_percent =(100 * (Fat[STATUS == "end_drug"] - Fat[STATUS == "baseline"])/Fat[STATUS == "baseline"]),
    BW_change_percent = (100 * (Weight[STATUS == "end_drug"] - Weight[STATUS == "baseline"])/Weight[STATUS == "baseline"]),
    .groups = "drop"
  ) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  drop_na()

###1. BW changes plot----

plot_47_5D_bw <- ggplot(echoMRI_data_47, aes(x = DRUG, y = BW_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "BW change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_bw

###2. adiposity index change plot ----

plot_47_5D_ai <- ggplot(echoMRI_data_47, aes(x = DRUG, y = adiposity_index_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "adiposity index change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_ai

###3. lean mass change plot----

plot_47_5D_lean <- ggplot(echoMRI_data_47, aes(x = DRUG, y = lean_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "lean mass change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_lean

###4. fat mass change plot ----

plot_47_5D_fat <- ggplot(echoMRI_data_47, aes(x = DRUG, y = fat_change, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "fat mass change (end-start)") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")

plot_47_5D_fat

###5. stats----
#1 does rtioxa-47 decreases bw after 5 days of injection?
#to check if variances differ between groups,
echoMRI_data_47 <- echoMRI_data_47 %>%
  mutate(DRUG = factor(DRUG))

leveneTest(adiposity_index_change ~ DRUG, data = echoMRI_data_47) #p = 0.3 > 0.05, so variances can be considered equal
leveneTest(fat_change ~ DRUG, data = echoMRI_data_47) #p = 0.3 > 0.05, so variances can be considered equal
leveneTest(lean_change ~ DRUG, data = echoMRI_data_47) #p = 0.009 < 0.05, so variances can not be considered equal
leveneTest(BW_change ~ DRUG, data = echoMRI_data_47) #p = 0.4 > 0.05, so variances can be considered equal
leveneTest(adiposity_index_change_percent ~ DRUG, data = echoMRI_data_47) #p = 0.3 > 0.05, so variances can be considered equal
leveneTest(fat_change_percent ~ DRUG, data = echoMRI_data_47) #p = 0.4 > 0.05, so variances can be considered equal
leveneTest(lean_change_percent ~ DRUG, data = echoMRI_data_47) #p = 0.01 < 0.05, so variances can not be considered equal
leveneTest(BW_change_percent ~ DRUG, data = echoMRI_data_47) #p = 0.4 > 0.05, so variances can be considered equal

summary(aov(BW_change ~ DRUG, data = echoMRI_data_47)) #p = 0.2
summary(aov(fat_change ~ DRUG, data = echoMRI_data_47)) #p = 0.5
oneway.test(lean_change ~ DRUG, data = echoMRI_data_47, var.equal = FALSE) #p = 0.5
summary(aov(adiposity_index_change ~ DRUG, data = echoMRI_data_47)) #p = 0.4
summary(aov(BW_change_percent ~ DRUG, data = echoMRI_data_47)) #p = 0.29
summary(aov(fat_change_percent ~ DRUG, data = echoMRI_data_47)) #p = 0.6
oneway.test(lean_change_percent ~ DRUG, data = echoMRI_data_47, var.equal = FALSE) #p = 0.4
summary(aov(adiposity_index_change_percent ~ DRUG, data = echoMRI_data_47)) #p = 0.6

#conclusion> There is no indication that RTIOXA-47 decreases changes in fat mass, lean mass or adiposity after 5 days there is a trend in decrease % bw change

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(12, 14)) %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  filter(str_detect(COMMENTS, regex("INJECTION_DAY_[1-5]", ignore_case = TRUE))) %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078, 41002, 41004, 
                41009, 41019, 41017, 41013, 41011) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle")) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  drop_na()

FI_data <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(12, 14)) %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  rename(DIET_FORMULA = DIET_FORMULA.x) %>% #There is no differences between columns x and y. 
  select(-DIET_FORMULA.y) %>% 
  filter(str_detect(COMMENTS, regex("INJECTION_DAY_[1-5]", ignore_case = TRUE))) %>% 
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078, 41002, 41004, 
                41009, 41019, 41017, 41013, 41011) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079, 41003, 41005, 
                41012, 41010, 41016, 41014, 41020) ~ "vehicle")) %>% 
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))) %>% 
  drop_na()

# left join FI_data with BW_data by ID and COMMENTS
FI_BW_joined <- FI_data %>%
  left_join(BW_data %>% select(ID, COMMENTS, BW), by = c("ID", "COMMENTS"))

# standarization of daily food intake per daily BW
FI_BW_joined <- FI_BW_joined %>%
  mutate(kcal_per_gBW = corrected_intake_kcal / BW)

# cumulative sum per id
FI_BW_joined_sum <- FI_BW_joined %>%
  group_by(ID, DRUG) %>%   # group by ID and DRUG
  summarise(
    total_kcal = sum(kcal_per_gBW, na.rm = TRUE),
    .groups = "drop"
  )

plot_47_foodintake <- ggplot(FI_BW_joined_sum, aes(x = DRUG, y = total_kcal, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(x = NULL, y = "Total kcal per g body weight over 5 days") +
  format.plot+
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "orange")) +
  theme(legend.position = "none")+
  geom_text_repel(aes(label = ID)) 
plot_47_foodintake

#clearlt both groups ate the same, cuek





