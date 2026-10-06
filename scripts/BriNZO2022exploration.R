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

# COLOR PALETTE ----
drug_colors <- c(
  "vehicle" = "gray60",
  "RTI_47" = "#E67E22"
)

drug_labels <- c(
  "vehicle" = "Vehicle",
  "RTI_47" = "RTI-47"
)

#BODY WEIGHT (BW) ANALYSIS----
## These c57bl6j and nzo mice were on chow ----

# assignation to the drugs based on food intake calculation sheet csv from Bri original data

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
  filter(DATE >= as.Date("2022-03-04")) #cohort 8 has acclimation, it seems like cohort 20 an 21 did not
   
BW_data  %>% 
  group_by(SEX,STRAIN,DRUG,COHORT) %>% #so cohort 21 are just females from both strains
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf)
  
  
BW_data_2 <- BW_data %>% 
  group_by( ID, STRAIN, SEX,DRUG) %>% 
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
  ) %>% 
  filter(STRAIN == "C57BL/6J") %>% 
  filter(SEX=="F") 

BW_data_2 %>% 
  group_by(DRUG) %>%
  summarise(
    n_ID = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  print(n = Inf)


BW_summary <- BW_data_2 %>%
  group_by(week_rel,DRUG) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(sum(!is.na(BW))),
    n = sum(!is.na(BW)),
    .groups = "drop"
  )


plot_BW <- ggplot() +
  
  # Individual mouse trajectories
  geom_line(
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
  
  # Individual measurements
  geom_point(
    data = BW_data_2,
    aes(
      x = week_rel,
      y = BW,
      color = DRUG
    ),
    alpha = 0.35,
    size = 1.5
  ) +
  
  # Group mean
  geom_line(
    data = BW_summary,
    aes(
      x = week_rel,
      y = mean_BW,
      group = DRUG,
      color = DRUG
    ),
    linewidth = 1.2
  ) +
  
  # Mean points
  geom_point(
    data = BW_summary,
    aes(
      x = week_rel,
      y = mean_BW,
      color = DRUG
    ),
    size = 3
  ) +
  
  # SEM
  geom_errorbar(
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
  
  #facet_wrap(~ COHORT) +
  
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


#stats----
BW_start_end <- BW_data %>% 
  filter(STRAIN == "C57BL/6J") %>% 
  filter(SEX=="F") %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, DRUG, STATUS, BW,STRAIN)

BW_start_end %>% 
  group_by(STATUS,DRUG) %>%
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


# Prepare summary data
BW_plot_summary <- BW_start_end %>%
  group_by(DRUG, STATUS) %>%
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


# Plot
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
#  facet_wrap(~COHORT)

plot_BW_start_end


#Does RTI_47 reduce BW during the 5-week treatment?----

BW_change <- BW_data %>% 
  filter(STATUS %in% c("start", "end")) %>% 
  select(ID, COHORT, SEX, DRUG, STATUS, BW,STRAIN,STATUS) %>% 
  pivot_wider(
    names_from = STATUS,
    values_from = BW
  ) %>% 
  mutate(
    BW_change = end - start,
    BW_change_percent = 100 * (end - start) / start
  )


# Summary
BW_change_summary <- BW_change %>%
  filter(
    STRAIN == "NZO/HlLtJ",
    SEX == "F"
  ) %>%
  group_by(DRUG,COHORT) %>%
  summarise(
    mean_change = mean(BW_change, na.rm = TRUE),
    sem_change = sd(BW_change, na.rm = TRUE) /
      sqrt(sum(!is.na(BW_change))),
    n = sum(!is.na(BW_change)),
    .groups = "drop"
  )

BW_change <- BW_change %>%
  mutate(
    DRUG = factor(
      DRUG,
      levels = c("vehicle", "RTI_47")
    )
  )

plot_BW_change <- ggplot(
  BW_change_summary,
  aes(x = DRUG, y = mean_change, fill = DRUG)
) +
  geom_col(width = 0.6, color = "black") +
  
  geom_errorbar(
    aes(
      ymin = mean_change - sem_change,
      ymax = mean_change + sem_change
    ),
    width = 0.2,
    linewidth = 0.7
  ) +
  
  geom_jitter(
    data = BW_change %>% 
      filter(STRAIN == "NZO/HlLtJ", SEX == "F"),
    aes(x = DRUG, y = BW_change),
    width = 0.10,
    size = 2.5,
    alpha = 0.7,
    color = "black"
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
  
  coord_cartesian(ylim = c(-5, 15)) +
  scale_y_continuous(breaks = seq(-5, 15, 5)) +
  
  labs(
    x = NULL,
    y = "Change in body weight (g)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )+
  facet_wrap(~COHORT)

plot_BW_change

BW_change_ttest <- BW_change %>% 
  filter(STRAIN == "NZO/HlLtJ") %>% 
  filter(SEX=="F") %>% 
  group_by(COHORT) %>% 
  t_test(
    BW_change ~ DRUG,
    paired = FALSE
  ) %>% 
  add_significance()

BW_change_ttest #there was a non-significant trend toward reduced BW gain in RTI_47-treated FEmales.


# BODY COMPOSITION ANALYSIS----

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
#  filter(STRAIN == "NZO/HlLtJ") %>% 
  filter(SEX=="F")        
  
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
  group_by(SEX,DRUG,n_measurement,STRAIN) %>%
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
    delta_lean_perc = lean_perc_end - lean_perc_start
  )

delta_bodycomp

delta_bodycomp %>% 
  group_by(SEX,DRUG,STRAIN) %>%
  summarise(n_ID = n_distinct(ID)) %>% 
  print(n = Inf) 

### fat mass CHANGE----

fatdelta_summary <- delta_bodycomp %>%
  filter(SEX == "F") %>%
  group_by(DRUG, STRAIN, SEX,COHORT) %>%
  summarise(
    mean_fatdelta = mean(delta_fat, na.rm = TRUE),
    sem_fatdelta  = sd(delta_fat, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_fat))),
    n = sum(!is.na(delta_fat)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    STRAIN = factor(
      STRAIN,
      levels = c("C57BL/6J", "NZO/HlLtJ")
    )
  )
plot_fatdelta <- ggplot(
  fatdelta_summary,
  aes(
    x = DRUG,
    y = mean_fatdelta,
    fill = STRAIN
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
    data = delta_bodycomp %>%
      filter(SEX == "F"),
    aes(
      x = DRUG,
      y = delta_fat,
      fill = STRAIN
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.7
    ),
    size = 2,
    shape = 21,
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
  
  scale_fill_manual(
    values = c(
      "C57BL/6J" = "gray60",
      "NZO/HlLtJ" = "pink"
    ),
    labels = c(
      "C57BL/6J" = "C57BL/6J",
      "NZO/HlLtJ" = "NZO/HlLtJ"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, by = 5)
  ) +
  
  labs(
    x = NULL,
    y = "Change in fat mass (g)\n(end - start)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top"
  )+
  facet_wrap(~COHORT)

plot_fatdelta  


fatdelta_ttest <- delta_bodycomp %>% 
  filter(SEX == "F") %>%
  group_by(DRUG,COHORT) %>% 
  t_test(delta_fat ~ STRAIN) %>% 
  add_significance()

fatdelta_ttest 

### fat mass at baseline and end ----

fat_start_end <- echoMRI_data %>%
  filter(
    STRAIN == "NZO/HlLtJ",
    SEX == "F",
    STATUS %in% c("start", "end")
  ) %>%
  mutate(
    x_position = case_when(
      DRUG == "vehicle" & STATUS == "start" ~ 1,
      DRUG == "vehicle" & STATUS == "end"   ~ 2,
      DRUG == "RTI_47"  & STATUS == "start" ~ 3,
      DRUG == "RTI_47"  & STATUS == "end"   ~ 4
    )
  )
fat_start_end_summary <- fat_start_end %>%
  group_by(DRUG, STATUS, x_position,COHORT) %>%
  summarise(
    mean_fat = mean(Fat, na.rm = TRUE),
    sem_fat = sd(Fat, na.rm = TRUE) /
      sqrt(sum(!is.na(Fat))),
    n = sum(!is.na(Fat)),
    .groups = "drop"
  )


plot_fat_start_end <- ggplot(
  fat_start_end_summary,
  aes(
    x = x_position,
    y = mean_fat,
    fill = DRUG
  )
) +
  
  # Mean bars
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Connect the same mouse from Start to End
  geom_line(
    data = fat_start_end,
    aes(
      x = x_position,
      y = Fat,
      group = ID
    ),
    color = "black",
    alpha = 0.35,
    linewidth = 0.6
  ) +
  
  # Individual mice
  geom_point(
    data = fat_start_end,
    aes(
      x = x_position,
      y = Fat
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_fat - sem_fat,
      ymax = mean_fat + sem_fat
    ),
    width = 0.15
  ) +
  
  # X-axis
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("Start", "End", "Start", "End"),
    limits = c(0.5, 4.5)
  ) +
  
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
  
  labs(
    x = NULL,
    y = "Fat mass (g)",
    fill = NULL
  ) +
  
  # Drug labels
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = "Vehicle",
    vjust = 2,
    size = 5
  ) +
  
  annotate(
    "text",
    x = 3.5,
    y = Inf,
    label = "RTI-47",
    vjust = 2,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    plot.margin = margin(
      t = 35,
      r = 10,
      b = 10,
      l = 10
    )
  )+
  facet_wrap(~COHORT)

plot_fat_start_end

### Statistics: Start vs End within each treatment group ----

fat_start_end_ttest <- fat_start_end %>%
  group_by(DRUG,COHORT) %>%
  t_test(
    Fat ~ STATUS,
    paired = TRUE
  ) %>%
  add_significance()

fat_start_end_ttest


### lean mass CHANGE----

leandelta_summary <- delta_bodycomp %>%
  filter(SEX == "F") %>%
  group_by(DRUG, STRAIN, SEX) %>%
  summarise(
    mean_leandelta = mean(delta_lean, na.rm = TRUE),
    sem_leandelta  = sd(delta_lean, na.rm = TRUE) /
      sqrt(sum(!is.na(delta_lean))),
    n = sum(!is.na(delta_lean)),
    .groups = "drop"
  ) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTI_47")),
    STRAIN = factor(
      STRAIN,
      levels = c("C57BL/6J", "NZO/HlLtJ")
    )
  )
plot_leandelta <- ggplot(
  leandelta_summary,
  aes(
    x = DRUG,
    y = mean_leandelta,
    fill = STRAIN
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
    data = delta_bodycomp %>%
      filter(SEX == "F"),
    aes(
      x = DRUG,
      y = delta_lean,
      fill = STRAIN
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.7
    ),
    size = 2,
    shape = 21,
    color = "black",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_leandelta - sem_leandelta,
      ymax = mean_leandelta + sem_leandelta
    ),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  
  scale_fill_manual(
    values = c(
      "C57BL/6J" = "gray60",
      "NZO/HlLtJ" = "pink"
    ),
    labels = c(
      "C57BL/6J" = "C57BL/6J",
      "NZO/HlLtJ" = "NZO/HlLtJ"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "vehicle" = "Vehicle",
      "RTI_47" = "RTI-47"
    )
  ) +
  
 # scale_y_continuous(
  #  limits = c(0, 15),
   # breaks = seq(0, 15, by = 5)
  #) +
  
  labs(
    x = NULL,
    y = "Change in lean mass (g)\n(end - start)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top"
  ) #+
  #facet_wrap(~COHORT)

plot_leandelta  


leandelta_ttest <- delta_bodycomp %>% 
  filter(SEX == "F") %>%
  group_by(DRUG) %>% 
  t_test(delta_lean ~ STRAIN) %>% 
  add_significance()

leandelta_ttest 

### lean mass at baseline and end ----

lean_start_end <- echoMRI_data %>%
  filter(
    STRAIN == "NZO/HlLtJ",
    SEX == "F",
    STATUS %in% c("start", "end")
  ) %>%
  mutate(
    x_position = case_when(
      DRUG == "vehicle" & STATUS == "start" ~ 1,
      DRUG == "vehicle" & STATUS == "end"   ~ 2,
      DRUG == "RTI_47"  & STATUS == "start" ~ 3,
      DRUG == "RTI_47"  & STATUS == "end"   ~ 4
    )
  )
lean_start_end_summary <- lean_start_end %>%
  group_by(DRUG, STATUS, x_position) %>%
  summarise(
    mean_lean = mean(Lean, na.rm = TRUE),
    sem_lean = sd(Lean, na.rm = TRUE) /
      sqrt(sum(!is.na(Lean))),
    n = sum(!is.na(Lean)),
    .groups = "drop"
  )


plot_lean_start_end <- ggplot(
  lean_start_end_summary,
  aes(
    x = x_position,
    y = mean_lean,
    fill = DRUG
  )
) +
  
  # Mean bars
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Connect the same mouse from Start to End
  geom_line(
    data = lean_start_end,
    aes(
      x = x_position,
      y = Lean,
      group = ID
    ),
    color = "black",
    alpha = 0.35,
    linewidth = 0.6
  ) +
  
  # Individual mice
  geom_point(
    data = lean_start_end,
    aes(
      x = x_position,
      y = Lean
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = mean_lean - sem_lean,
      ymax = mean_lean + sem_lean
    ),
    width = 0.15
  ) +
  
  # X-axis
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("Start", "End", "Start", "End"),
    limits = c(0.5, 4.5)
  ) +
  
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
  
  labs(
    x = NULL,
    y = "Lean mass (g)",
    fill = NULL
  ) +
  
  # Drug labels
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = "Vehicle",
    vjust = 2,
    size = 5
  ) +
  
  annotate(
    "text",
    x = 3.5,
    y = Inf,
    label = "RTI-47",
    vjust = 2,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    plot.margin = margin(
      t = 35,
      r = 10,
      b = 10,
      l = 10
    )
  )#+
  #facet_wrap(~COHORT)

plot_lean_start_end

### Statistics: Start vs End within each treatment group ----

lean_start_end_ttest <- lean_start_end %>%
  group_by(DRUG) %>%
  t_test(
    Lean ~ STATUS,
    paired = TRUE
  ) %>%
  add_significance()

lean_start_end_ttest


