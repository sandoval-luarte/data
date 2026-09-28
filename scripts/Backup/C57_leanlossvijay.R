# This script aims to explore changes lean mass in middle age C57 after different stages of feeding:
#from peak of obesity to acute body weight loss

#libraries
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

echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 2) %>% # Just C57 males and females
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"))%>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP,SEX, DIET_FORMULA) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      Date == as.Date("2025-03-07") ~ "peak obesity",
      Date == as.Date("2025-04-21") ~ "BW loss",
      Date == as.Date("2025-06-05") ~ "BW maintenance",
      Date %in% as.Date(c("2025-09-11", "2025-09-10","2025-09-05","2025-09-04",
                          "2025-09-02","2025-09-01","2025-08-28","2025-08-27")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>% 
  filter( STATUS %in% c( "BW loss", "peak obesity")) %>% 
  filter(DIET_FORMULA=="D12451i")

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c( "peak obesity", "BW loss")))


# DATA FOR PLOT


plot_data <- echoMRI_data %>%
  mutate(
    PlotGroup = case_when(
      STATUS == "peak obesity" ~ "all",
      STATUS == "BW loss" ~ GROUP
    )
  ) %>%
  filter(!is.na(PlotGroup))


# Make sure the order is:
# all -> restricted -> ad lib
plot_data$PlotGroup <- factor(
  plot_data$PlotGroup,
  levels = c("all", "restricted", "ad lib")
)

# SUMMARY DATA


summary_data <- plot_data %>%
  group_by(SEX, STATUS, PlotGroup) %>%
  summarise(
    mean = mean(Lean, na.rm = TRUE),
    SE = sd(Lean, na.rm = TRUE) / sqrt(sum(!is.na(Lean))),
    .groups = "drop"
  )


bar_data <- summary_data %>%
  mutate(
    xmin = case_when(
      STATUS == "peak obesity" ~ 0.75,
      STATUS == "BW loss" & PlotGroup == "ad lib" ~ 1.75,
      STATUS == "BW loss" & PlotGroup == "restricted" ~ 1.75
    ),
    
    xmax = case_when(
      STATUS == "peak obesity" ~ 1.25,
      STATUS == "BW loss" ~ 2.25
    )
  )


# PLOT 1----

plot <- ggplot() +
# INDIVIDUAL ID TRAJECTORIES

geom_line(
  data = plot_data,
  aes(
    x = ifelse(STATUS == "peak obesity", 1, 2),
    y = Lean,
    group = ID
  ),
  color = "black",
  linewidth = 0.7,
  alpha = 0.9
) +
geom_col(
  data = bar_data %>%
    filter(STATUS == "peak obesity"),
  aes(
    x = 1,
    y = mean,
    fill = PlotGroup
  ),
  width = 0.55,
  color = "#8C6F00",
  linewidth = 0.8
) +
  geom_col(
    data = bar_data %>%
      filter(
        STATUS == "BW loss",
        PlotGroup == "ad lib"
      ),
    aes(
      x = 2,
      y = mean,
      fill = PlotGroup
    ),
    width = 0.55,
    color = "#317EC2FF",
    linewidth = 0.8
  ) +
  geom_col(
    data = bar_data %>%
      filter(
        STATUS == "BW loss",
        PlotGroup == "restricted"
      ),
    aes(
      x = 2,
      y = mean,
      fill = PlotGroup
    ),
    width = 0.55,
    color = "#317EC2FF",
    linewidth = 0.8
  ) +
geom_errorbar(
  data = bar_data %>%
    filter(STATUS == "peak obesity"),
  aes(
    x = 1,
    ymin = mean - SE,
    ymax = mean + SE
  ),
  width = 0.15,
  linewidth = 1
) +
geom_errorbar(
  data = bar_data %>%
    filter(STATUS == "BW loss"),
  aes(
    x = 2,
    ymin = mean - SE,
    ymax = mean + SE
  ),
  width = 0.15,
  linewidth = 1
) +
geom_point(
  data = plot_data,
  aes(
    x = ifelse(STATUS == "peak obesity", 1, 2),
    y = Lean
  ),
  shape = 21,
  fill = "white",
  color = "#8AA7D8",
  size = 2.5,
  stroke = 0.9,
  position = position_jitter(
    width = 0.04,
    height = 0
  )
) +
scale_x_continuous(
  breaks = c(1, 2),
  labels = c(
    "peak obesity",
    "BW loss"
  ),
  limits = c(0.5, 2.5)
) +
scale_fill_manual(
  values = c(
    "all" = "#F5F1D5",
    "restricted" = "#BFE3E0",
    "ad lib" = "#8B001F"
  ),
  labels = c(
    "all" = "all",
    "restricted" = "restricted",
    "ad lib" = "ad libitum"
  )
) +
facet_wrap(~SEX) +

labs(
  x = "STATUS",
  y = "Lean mass (g)",
  fill = "Group"
) +
theme_classic() +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = "Helvetica",
      size = 14,
      face = "plain"
    ),
    
    axis.text = element_text(
      family = "Helvetica",
      size = 13,
      color = "black"
    ),
    
    axis.title = element_text(
      family = "Helvetica",
      size = 14,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.title = element_text(
      family = "Helvetica",
      size = 13
    ),
    
    legend.text = element_text(
      family = "Helvetica",
      size = 12
    ),
    
    panel.grid = element_blank(),
    
    panel.border = element_blank(),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    panel.spacing = unit(1.2, "lines")
  )


plot

# CALCULATE LEAN MASS AS % OF PEAK OBESITY


plot_data <- echoMRI_data %>%
  group_by(ID) %>%
  mutate(
    peak_lean = Lean[STATUS == "peak obesity"][1],
    Lean_percent_peak = (Lean / peak_lean) * 100
  ) %>%
  ungroup()

# DATA FOR PLOT


plot_data <- plot_data %>%
  mutate(
    PlotGroup = case_when(
      STATUS == "peak obesity" ~ "all",
      STATUS == "BW loss" ~ GROUP
    )
  ) %>%
  filter(!is.na(PlotGroup))


plot_data$PlotGroup <- factor(
  plot_data$PlotGroup,
  levels = c("all", "restricted", "ad lib")
)

# SUMMARY DATA


summary_data <- plot_data %>%
  group_by(SEX, STATUS, PlotGroup) %>%
  summarise(
    mean = mean(Lean_percent_peak, na.rm = TRUE),
    SE = sd(Lean_percent_peak, na.rm = TRUE) /
      sqrt(sum(!is.na(Lean_percent_peak))),
    .groups = "drop"
  )

# PLOT


plot_perc <- ggplot() +
  
  # ----------------------------------------------------------
# INDIVIDUAL ID TRAJECTORIES
# ----------------------------------------------------------

geom_line(
  data = plot_data,
  aes(
    x = ifelse(STATUS == "peak obesity", 1, 2),
    y = Lean_percent_peak,
    group = ID
  ),
  color = "black",
  linewidth = 0.7,
  alpha = 0.9
) +
  
  # ----------------------------------------------------------
# PEAK OBESITY — ALL
# ----------------------------------------------------------

geom_col(
  data = summary_data %>%
    filter(
      STATUS == "peak obesity",
      PlotGroup == "all"
    ),
  aes(
    x = 1,
    y = mean,
    fill = PlotGroup
  ),
  width = 0.55,
  color = "#8C6F00",
  linewidth = 0.8
) +
  
  # ----------------------------------------------------------
# BW LOSS — AD LIBITUM
# ----------------------------------------------------------

geom_col(
  data = summary_data %>%
    filter(
      STATUS == "BW loss",
      PlotGroup == "ad lib"
    ),
  aes(
    x = 2,
    y = mean,
    fill = PlotGroup
  ),
  width = 0.55,
  color = "#317EC2FF",
  linewidth = 0.8
) +
  
  # ----------------------------------------------------------
# BW LOSS — RESTRICTED
# ----------------------------------------------------------

geom_col(
  data = summary_data %>%
    filter(
      STATUS == "BW loss",
      PlotGroup == "restricted"
    ),
  aes(
    x = 2,
    y = mean,
    fill = PlotGroup
  ),
  width = 0.55,
  color = "#317EC2FF",
  linewidth = 0.8
) +
  
  # ----------------------------------------------------------
# SEM — PEAK OBESITY
# ----------------------------------------------------------

geom_errorbar(
  data = summary_data %>%
    filter(STATUS == "peak obesity"),
  aes(
    x = 1,
    ymin = mean - SE,
    ymax = mean + SE
  ),
  width = 0.15,
  linewidth = 1
) +
  
  # ----------------------------------------------------------
# SEM — BW LOSS
# ----------------------------------------------------------

geom_errorbar(
  data = summary_data %>%
    filter(STATUS == "BW loss"),
  aes(
    x = 2,
    ymin = mean - SE,
    ymax = mean + SE
  ),
  width = 0.15,
  linewidth = 1
) +
  
  # ----------------------------------------------------------
# INDIVIDUAL ANIMALS
# ----------------------------------------------------------

geom_point(
  data = plot_data,
  aes(
    x = ifelse(STATUS == "peak obesity", 1, 2),
    y = Lean_percent_peak
  ),
  shape = 21,
  fill = "white",
  color = "#8AA7D8",
  size = 2.5,
  stroke = 0.9,
  position = position_jitter(
    width = 0.04,
    height = 0
  )
) +
  
  # ----------------------------------------------------------
# X AXIS
# ----------------------------------------------------------

scale_x_continuous(
  breaks = c(1, 2),
  labels = c(
    "peak obesity",
    "BW loss"
  ),
  limits = c(0.5, 2.5)
) +
  
  # ----------------------------------------------------------
# Y AXIS
# ----------------------------------------------------------

scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0, 0.05))
) +
  
  # ----------------------------------------------------------
# COLORS
# ----------------------------------------------------------

scale_fill_manual(
  values = c(
    "all" = "#F5F1D5",
    "restricted" = "#BFE3E0",
    "ad lib" = "#8B001F"
  ),
  labels = c(
    "all" = "all",
    "restricted" = "restricted",
    "ad lib" = "ad libitum"
  )
) +
  
  # ----------------------------------------------------------
# FACETS
# ----------------------------------------------------------

facet_wrap(~SEX) +
  
  # ----------------------------------------------------------
# LABELS
# ----------------------------------------------------------

labs(
  x = "STATUS",
  y = "Lean mass (% of peak obesity)",
  fill = "Group"
) +
  
  # ----------------------------------------------------------
# THEME
# ----------------------------------------------------------

theme_classic() +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = "Arial",
      size = 14,
      face = "plain"
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 13,
      color = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 14,
      color = "black"
    ),
    
    axis.text.x = element_text(
      family = "Arial",
      size = 13,
      angle = 45,
      hjust = 1
    ),
    
    legend.title = element_text(
      family = "Arial",
      size = 13
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 12
    ),
    
    panel.grid = element_blank(),
    
    panel.border = element_blank(),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    panel.spacing = unit(1.2, "lines")
  )

plot_perc

plot <- plot +
  labs(tag = "A")

plot_perc <- plot_perc +
  labs(tag = "B")

combined_plot <- (plot | plot_perc) 

combined_plot


stats_data <- echoMRI_data %>%
  mutate(
    STATUS = factor(
      STATUS,
      levels = c("peak obesity", "BW loss")
    ),
    GROUP = factor(
      GROUP,
      levels = c("ad lib", "restricted")
    ),
    SEX = factor(SEX),
    ID = factor(ID)
  )
model <- lmer(
  Lean ~ STATUS * GROUP + (1 | ID),
  data = stats_data
)

anova(model)

emm <- emmeans(
  model,
  ~ STATUS * GROUP
)

emm
pairs(
  emmeans(model, ~ GROUP | STATUS)
)

pairs(
  emmeans(model, ~ STATUS | GROUP)
)

model_sex <- lmer(
  Lean ~ STATUS * GROUP * SEX + (1 | ID),
  data = stats_data
)

anova(model_sex, ddf = "Kenward-Roger")

library(lmerTest)
class(model_sex)
anova(model_sex)

library(lme4)
library(lmerTest)

# Refit AFTER loading lmerTest
model_sex <- lmer(
  Lean ~ STATUS * GROUP * SEX + (1 | ID),
  data = stats_data
)

# Check
class(model_sex)

# ANOVA
anova(model_sex)


# ============================================================
# THIRD PLOT: SEX COLLAPSED
# ============================================================

plot_sex_collapsed <- plot_data %>%
  mutate(
    PlotGroup = case_when(
      STATUS == "peak obesity" ~ "Peak obesity",
      STATUS == "BW loss" & GROUP == "restricted" ~ "Restricted",
      STATUS == "BW loss" & GROUP == "ad lib" ~ "Ad lib"
    )
  ) %>%
  filter(!is.na(PlotGroup)) %>%
  mutate(
    PlotGroup = factor(
      PlotGroup,
      levels = c("Peak obesity", "Restricted", "Ad lib")
    )
  )

# Summary statistics
summary_sex_collapsed <- plot_sex_collapsed %>%
  group_by(PlotGroup) %>%
  summarise(
    mean_lean = mean(Lean_percent_peak, na.rm = TRUE),
    sem_lean = sd(Lean_percent_peak, na.rm = TRUE) /
      sqrt(sum(!is.na(Lean_percent_peak))),
    .groups = "drop"
  )

# Plot
plot_sex_collapsed <- ggplot() +
  
  # Bars
  geom_col(
    data = summary_sex_collapsed,
    aes(x = PlotGroup, y = mean_lean),
    fill = "white",
    color = "black",
    width = 0.65
  ) +
  
  # SEM
  geom_errorbar(
    data = summary_sex_collapsed,
    aes(
      x = PlotGroup,
      ymin = mean_lean - sem_lean,
      ymax = mean_lean + sem_lean
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual points
  geom_jitter(
    data = plot_sex_collapsed,
    aes(x = PlotGroup, y = Lean_percent_peak),
    width = 0.08,
    size = 1.5,
    color = "black",
    alpha = 0.8
  ) +
  
  # 100% reference line
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  labs(
    x = NULL,
    y = "Lean mass (% of peak obesity)"
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.text = element_text(
      family = "Arial",
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      family = "Arial",
      size = 13,
      color = "black"
    ),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

plot_sex_collapsed 

plot <- plot + labs(tag = "A")
plot_perc <- plot_perc + labs(tag = "B")
plot_sex_collapsed <- plot_sex_collapsed + labs(tag = "C")

combined_plot <- plot + plot_perc + plot_sex_collapsed

combined_plot

library(dplyr)

lean_change <- stats_data %>%
  group_by(ID) %>%
  mutate(
    peak_lean = Lean[STATUS == "peak obesity"],
    Lean_percent_change = ((Lean - peak_lean) / peak_lean) * 100
  ) %>%
  ungroup()
lean_change_bwloss <- lean_change %>%
  filter(STATUS == "BW loss")
model_percent <- lm(
  Lean_percent_change ~ GROUP,
  data = lean_change_bwloss
)

anova(model_percent)
summary(model_percent)
library(emmeans)

emmeans(model_percent, ~ GROUP)
pairs(emmeans(model_percent, ~ GROUP))

model_percent_sex <- lm(
  Lean_percent_change ~ GROUP * SEX,
  data = lean_change_bwloss
)

anova(model_percent_sex)
emmeans(model_percent_sex, ~ GROUP | SEX)
pairs(emmeans(model_percent_sex, ~ GROUP | SEX))
