# This script aims to explore changes in body weight in NZO and c57 after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#libraries----
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr)
library(purrr)
library(broom)
library(Hmisc)
library(lme4)
library(emmeans)
library(patchwork)


#BW CSV data import RTIOXA 47 ----
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 
                7876, 7879, 7880, 7881,7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,
                7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871,
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ),
    day_rel = DATE - first(DATE)
  ) %>%
  mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 202 & GROUP == "restricted" ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 183 & GROUP == "restricted"~ "peak obesity",
      STRAIN == "C57BL6/J" & day_rel == 204 & GROUP == "ad lib" ~ "peak obesity",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 241 & GROUP == "restricted"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 218 & GROUP == "restricted"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "M" & day_rel == 239 & GROUP == "ad lib"~ "BW loss",
      STRAIN == "C57BL6/J" & SEX == "F" & day_rel == 218 & GROUP == "ad lib"~ "BW loss",
      STRAIN == "NZO/HlLtJ" & day_rel == 161 ~ "BW loss",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance", 
      STRAIN == "NZO/HlLtJ" &  COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    ) ) #%>% 
  #filter(STRAIN == "NZO/HlLtJ" | (STRAIN == "C57BL6/J" & DIET_FORMULA == "D12451i"))

n_distinct(BW_data$ID) #9 NZO in total and 4 C57 in total

highlight_data <- BW_data %>%
  filter(STATUS == "peak obesity")

highlight_data_2 <- BW_data %>%
  filter(STATUS == "BW loss")

highlight_data_3 <- BW_data %>%
  filter(STATUS == "BW maintenance")

highlight_data_4 <- BW_data %>%
  filter(STATUS == "BW regain")

#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14))


#plot 1 description of BW over time per ID----
plot <- BW_data %>% 
  ggplot(aes(day_rel, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~STRAIN*GROUP) +
  geom_smooth() +
  labs(
    x = "Days relative to baseline",
    y = "BW (grams)"
  ) +
  geom_vline(data = highlight_data, aes(xintercept = day_rel),
             linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(data = highlight_data_2, aes(xintercept = day_rel),
             linetype = "dashed", color = "green", alpha = 0.7) +
  geom_vline(data = highlight_data_3, aes(xintercept = day_rel),
             linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_vline(data = highlight_data_4, aes(xintercept = day_rel),
             linetype = "dashed", color = "orange", alpha = 0.7) +
  format.plot

plot

# Make STATUS an ordered factor
BW_data <- BW_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain"))) %>% 
  filter(!is.na(STATUS))

BW_data %>%
  group_by(STATUS) %>%
  summarise(n_ID = n_distinct(ID)) #this is good


# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX, DIET_FORMULA
summary_data <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP,SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

#plot 2 BW in each STATUS collapsed by DIET_FORMULA----
ggplot() +
  # Individual points, smaller and more transparent
  geom_point(data = BW_data, 
             aes(x = STATUS, y = BW, color = SEX, group = ID),
             size = 1.2, alpha = 0.3) +
  geom_line(data = BW_data, 
            aes(x = STATUS, y = BW, color = DRUG, group = ID),
            alpha = 0.15) +
  
  # Mean ± SEM as shaded ribbon
  geom_ribbon(data = summary_data,
              aes(x = STATUS, ymin = mean_BW - sem_BW, ymax = mean_BW + sem_BW,
                  group = interaction(GROUP, SEX), fill = SEX),
              alpha = 0.2) +
  geom_line(data = summary_data,
            aes(x = STATUS, y = mean_BW, group = interaction(GROUP, SEX), color = SEX),
            size = 1) +
  
  # Facet by STRAIN (rows) × GROUP  (columns)
  facet_grid(~ STRAIN*GROUP, scales = "free_x") +
  
  scale_color_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  scale_fill_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  
  labs(
    x = "Status",
    y = "Body weight (grams)",
    color = "Sex",
    fill = "Sex"
  ) +
  format.plot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX using bw_rel
summary_data_rel <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP, SEX) %>%
  summarise(
    mean_bw_rel = mean(bw_rel, na.rm = TRUE),
    sem_bw_rel = sd(bw_rel, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

#plot 3 % BW gain in each STATUS collapsed by DIET_FORMULA----

ggplot() +
  # Individual points, smaller and more transparent
  geom_point(data = BW_data, 
             aes(x = STATUS, y = bw_rel, color = SEX, group = ID),
             size = 1.2, alpha = 0.3) +
  geom_line(data = BW_data, 
            aes(x = STATUS, y = bw_rel, group = ID),
            alpha = 0.15) +
  
  # Mean ± SEM as shaded ribbon
  geom_ribbon(data = summary_data_rel,
              aes(x = STATUS, ymin = mean_bw_rel - sem_bw_rel, ymax = mean_bw_rel + sem_bw_rel,
                  group = interaction(GROUP, SEX), fill = SEX),
              alpha = 0.2) +
  geom_line(data = summary_data_rel,
            aes(x = STATUS, y = mean_bw_rel, group = interaction(GROUP, SEX), color = SEX),
            size = 1) +
  
  # Facet by STRAIN (rows) × GROUP (columns)
  facet_grid(~ STRAIN*GROUP, scales = "free_x") +
  
  # Colors
  scale_color_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  scale_fill_manual(values = c("F" = "#C03830FF", "M" = "#317EC2FF")) +
  
  # Labels and theme
  labs(
    x = "Status",
    y = "Body weight gain (%)",
    color = "Sex",
    fill = "Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  format.plot


#plot 4 BW in all females together ----


# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX
summary_data_rel <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP, SEX) %>%
  summarise(
    mean_BW = mean(BW, na.rm = TRUE),
    sem_BW = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define custom colors per STRAIN (still harmonious with your sex colors)
strain_colors <- c(
  "C57BL6/J" = "#C03830FF",   # red tone
  "NZO/HlLtJ"   = "#317EC2FF"    # blue tone
)

# Plot
ggplot() +
  # Individual points
  geom_point(data = BW_data,
             aes(x = STATUS, y = BW, color = STRAIN, group = ID),
             size = 1.2, alpha = 0.3) +
  
  # Individual trajectories
  geom_line(data = BW_data,
            aes(x = STATUS, y = BW, group = ID, color = STRAIN),
            alpha = 0.15) +
  
  # Mean ± SEM ribbons per STRAIN
  geom_ribbon(data = summary_data_rel,
              aes(x = STATUS,
                  ymin = mean_BW - sem_BW,
                  ymax = mean_BW + sem_BW,
                  group = interaction(STRAIN, GROUP, SEX),
                  fill = STRAIN),
              alpha = 0.15) +
  
  # Mean lines per STRAIN
  geom_line(data = summary_data_rel,
            aes(x = STATUS, y = mean_BW,
                group = interaction(STRAIN, GROUP, SEX),
                color = STRAIN),
            size = 1) +
  
  # Facet by GROUP × SEX (not by strain)
  facet_grid(GROUP ~ SEX, scales = "free_x") +
  
  # Custom colors for strains
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  labs(
    x = "Status",
    y = "Body weight (grams)",
    color = "Strain",
    fill = "Strain"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  format.plot

#plot 5 % BW gain in all females together ----

# Compute mean ± SEM per STRAIN, STATUS, GROUP, SEX
summary_data_rel <- BW_data %>%
  group_by(STRAIN, STATUS, GROUP, SEX) %>%
  summarise(
    mean_bw_rel = mean(bw_rel, na.rm = TRUE),
    sem_bw_rel = sd(bw_rel, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define custom colors per STRAIN (still harmonious with your sex colors)
strain_colors <- c(
  "C57BL6/J" = "#C03830FF",   # red tone
  "NZO/HlLtJ"   = "#317EC2FF"    # blue tone
)

# Plot
ggplot() +
  # Individual points
  geom_point(data = BW_data,
             aes(x = STATUS, y = bw_rel, color = STRAIN, group = ID),
             size = 1.2, alpha = 0.3) +
  
  # Individual trajectories
  geom_line(data = BW_data,
            aes(x = STATUS, y = bw_rel, group = ID, color = STRAIN),
            alpha = 0.15) +
  
  # Mean ± SEM ribbons per STRAIN
  geom_ribbon(data = summary_data_rel,
              aes(x = STATUS,
                  ymin = mean_bw_rel - sem_bw_rel,
                  ymax = mean_bw_rel + sem_bw_rel,
                  group = interaction(STRAIN, GROUP, SEX),
                  fill = STRAIN),
              alpha = 0.15) +
  
  # Mean lines per STRAIN
  geom_line(data = summary_data_rel,
            aes(x = STATUS, y = mean_bw_rel,
                group = interaction(STRAIN, GROUP, SEX),
                color = STRAIN),
            size = 1) +
  
  # Facet by GROUP × SEX (not by strain)
  facet_grid(GROUP ~ SEX, scales = "free_x") +
  
  # Custom colors for strains
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  
  labs(
    x = "Status",
    y = "Body weight change (%)",
    color = "Strain",
    fill = "Strain"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  format.plot

# plot 6 RTI47 PLOT BW% GAIN ----

# Step 1. Compute % change from BW maintenance → BW regain


BW_delta <- BW_data %>%
  filter(STATUS %in% c("BW maintenance", "BW regain")) %>%
  filter(DRUG %in% c("vehicle", "RTIOXA_47")) %>%
  select(ID, DRUG, STATUS, BW,STRAIN,SEX) %>%
  pivot_wider(names_from = STATUS, values_from = BW) %>%
  filter(!is.na(`BW maintenance`) & !is.na(`BW regain`)) %>%
  group_by(ID,SEX,STRAIN) %>% 
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    bw_rel = 100 * (`BW regain` - `BW maintenance`) / `BW maintenance` # % change
  )

# Step 2. Plot mean ± SEM bars with individual mice
ggplot(BW_delta, aes(x = DRUG, y = bw_rel, fill = DRUG)) +
  # Bars for mean ± SEM
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  # Individual mice points
  geom_point(position = position_jitter(width = 0.1), size = 2, color = "black")  +
  # T-test statistics
  stat_compare_means(
    comparisons = list(c("vehicle", "RTIOXA_47")),
    method = "t.test",
    label = "p.signif",
    label.y = max(BW_delta$bw_rel) + 2
  ) +
  labs(
    x = "",
    y = "% BW change from BW maintenance",
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray70", "#8DA0CB"))+
  facet_grid(~STRAIN*SEX)

#BW CSV data import RTIOXA 43 x 5 days in a row  ----

BW_data_43 <- read_csv("../data/BW.csv") %>% 
  filter(COHORT == 9) %>% # 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
      ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
      ID %in% c(8095, 8100, 8103) ~ "Vehicle"),
    day_rel = DATE - first(DATE)) %>% 
  filter(as.Date(DATE) >= as.Date("2025-04-16"))

n_distinct(BW_data_43$ID) #9 animals in total, 3 veh, 3 rtioxa 43 from medchem and 3 rtioxa 43 from yanan


# Step 1. Prepare baseline vs endpoint using bw_rel
BW_compare <- BW_data_43 %>%
  filter(COMMENTS %in% c("DAY_1_INJECTIONS", "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION")) %>%
  mutate(COMMENTS = recode(COMMENTS,
                           "DAY_1_INJECTIONS" = "baseline",
                           "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION" = "endpoint")) %>%
  select(ID, DRUG, COMMENTS, bw_rel) %>%
  pivot_wider(names_from = COMMENTS, values_from = bw_rel) %>%
  filter(!is.na(endpoint)) %>%
  mutate(DRUG = factor(DRUG, levels = c("Vehicle", "RTI_43_Y", "RTI_43_M")))

# Step 2. Remove groups with <2 animals
valid_groups <- BW_compare %>%
  group_by(DRUG) %>%
  filter(n() >= 2) %>%
  ungroup()

# plot 7 RTI43 PLOT BW% GAIN ----

# Step 3. Plot with raw data, mean ± SEM bars, and statistics
ggplot(valid_groups, aes(x = DRUG, y = endpoint, fill = DRUG)) +
  # Mean ± SEM bars
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  # Individual points
  geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
  # T-test statistics
  stat_compare_means(
    comparisons = list(c("Vehicle", "RTI_43_Y"),
                       c("Vehicle", "RTI_43_M")),
    method = "t.test",
    label = "p.signif",
    label.y = c(max(valid_groups$endpoint) + 5,
                max(valid_groups$endpoint) + 10),
    tip.length = 0.02
  ) +
  labs(
    x = "",
    y = "% of BW gain from baseline)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray70", "#66C2A5", "darkgreen"))

# --- Plot 6: RTIOXA-47 % BW change (Panel A) ----
BW_delta <- BW_data %>%
  filter(STATUS %in% c("BW maintenance", "BW regain")) %>%
  filter(DRUG %in% c("vehicle", "RTIOXA_47")) %>%
  select(ID, DRUG, STATUS, BW, STRAIN, SEX) %>%
  pivot_wider(names_from = STATUS, values_from = BW) %>%
  filter(!is.na(`BW maintenance`) & !is.na(`BW regain`)) %>%
  group_by(ID, SEX, STRAIN) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    bw_rel = 100 * (`BW regain` - `BW maintenance`) / `BW maintenance`
  )

plot_A <- ggplot(BW_delta, aes(x = DRUG, y = bw_rel, fill = DRUG)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  geom_point(position = position_jitter(width = 0.1), size = 2, color = "black") +
  stat_compare_means(
    comparisons = list(c("vehicle", "RTIOXA_47")),
    method = "t.test",
    label = "p.signif",
    label.y = max(BW_delta$bw_rel, na.rm = TRUE) + 2
  ) +
  labs(
    x = "",
    y = "% BW change from BW maintenance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 18)
  ) +
  scale_fill_manual(values = c("vehicle" = "white", "RTIOXA_47" = "#8DA0CB")) +
  facet_grid(~ STRAIN * SEX) +
  coord_cartesian(ylim = c(0, 40))  # <<--- unified Y axis


# --- Plot 7: RTIOXA-43 % BW gain (Panel B) ----
BW_data_43 <- read_csv("../data/BW.csv") %>%
  filter(COHORT == 9) %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTIOXA_43_donated",
      ID %in% c(8097, 8098, 8101) ~ "RTIOXA_43_Medchem",
      ID %in% c(8095, 8100, 8103) ~ "vehicle"
    ),
    day_rel = DATE - first(DATE),
    STRAIN = "C57BL6/J"
  ) %>%
  filter(as.Date(DATE) >= as.Date("2025-04-16"))

BW_compare <- BW_data_43 %>%
  filter(COMMENTS %in% c("DAY_1_INJECTIONS", "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION")) %>%
  mutate(COMMENTS = recode(COMMENTS,
                           "DAY_1_INJECTIONS" = "baseline",
                           "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION" = "endpoint")) %>%
  select(ID, DRUG, COMMENTS, bw_rel, STRAIN, SEX) %>%
  pivot_wider(names_from = COMMENTS, values_from = bw_rel) %>%
  filter(!is.na(endpoint)) %>%
  mutate(DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_donated", "RTIOXA_43_Medchem")))

valid_groups <- BW_compare %>%
  group_by(DRUG) %>%
  filter(n() >= 2) %>%
  ungroup()

plot_B <- ggplot(valid_groups, aes(x = DRUG, y = endpoint, fill = DRUG)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  geom_point(position = position_jitter(width = 0.1), size = 2, color = "black") +
  stat_compare_means(
    comparisons = list(
      c("vehicle", "RTIOXA_43_donated"),
      c("vehicle", "RTIOXA_43_Medchem")
    ),
    method = "t.test",
    label = "p.signif",
    label.y = c(max(valid_groups$endpoint) + 5,
                max(valid_groups$endpoint) + 10),
    tip.length = 0.02
  ) +
  labs(
    x = "",
    y = "% BW gain from baseline"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 18)
  ) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_Medchem" = "darkgreen"
  )) +
  facet_grid(~ STRAIN * SEX) +
  coord_cartesian(ylim = c(0, 40))  # <<--- unified Y axis


# --- Combine plots A and B ----
library(patchwork)

(plot_B | plot_A) + plot_annotation(tag_levels = "A")

