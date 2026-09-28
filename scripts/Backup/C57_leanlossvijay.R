# This script aims to explore changes lean mass in middle age C57 from peak of obesity to acute body weight loss

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

stats_data <- echoMRI_data %>%
  filter(
    STATUS %in% c("peak obesity", "BW loss"),
    DIET_FORMULA == "D12451i"
  ) %>%
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

lean_change <- stats_data %>%
  group_by(ID) %>%
  mutate(
    peak_lean = Lean[STATUS == "peak obesity"][1],
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

# Panel C: Lean mass % change from peak obesity, sex collapsed

# Calculate % change relative to each animal's peak-obesity lean mass
lean_change <- stats_data %>%
  group_by(ID) %>%
  mutate(
    peak_lean = Lean[STATUS == "peak obesity"][1],
    Lean_percent_change = ((Lean - peak_lean) / peak_lean) * 100
  ) %>%
  ungroup()

# Keep only BW loss
lean_change_bwloss <- lean_change %>%
  filter(STATUS == "BW loss") %>%
  mutate(
    GROUP = factor(
      GROUP,
      levels = c("restricted", "ad lib")
    )
  )

# Summary statistics
plot_data_C <- lean_change_bwloss %>%
  group_by(GROUP) %>%
  summarise(
    mean = mean(Lean_percent_change, na.rm = TRUE),
    SEM = sd(Lean_percent_change, na.rm = TRUE) /
      sqrt(sum(!is.na(Lean_percent_change))),
    .groups = "drop"
  )

# Colors
group_colors <- c(
  "ad lib" = "#6BAED6",      
  "restricted" = "#7FBF7B"  
)

# Panel C
pC <- ggplot(plot_data_C,
             aes(x = GROUP, y = mean, fill = GROUP)) +
  
  # Bars
  geom_col(
    width = 0.65,
    color = "black"
  ) +
  
  # Error bars
  geom_errorbar(
    aes(
      ymin = mean - SEM,
      ymax = mean + SEM
    ),
    width = 0.15,
    linewidth = 0.6,
    color = "black"
  ) +
  
  # Individual animals
  geom_jitter(
    data = lean_change_bwloss,
    aes(
      x = GROUP,
      y = Lean_percent_change
    ),
    width = 0.08,
    size = 1.5,
    color = "black",
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  
  # Peak obesity reference (0%)
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  
  # Significance bar for restricted group
  annotate(
    "segment",
    x = 1, xend = 1,
    y = 0, yend = -11
  ) +
  
  annotate(
    "segment",
    x = 0.92, xend = 1.08,
    y = 0, yend = 0
  ) +
  
  annotate(
    "segment",
    x = 0.92, xend = 1.08,
    y = -11, yend = -11
  ) +
  
  annotate(
    "text",
    x = 1.15,
    y = -5.5,
    label = "p < 0.0001",
    angle = 90,
    size = 3.5,
    family = "Arial"
  ) +
  
  labs(
    x = NULL,
    y = "Lean mass change from peak obesity (%)"
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_x_discrete(
    labels = c(
      "restricted" = "Restricted",
      "ad lib" = "Ad libitum"
    )
  ) +
  
  theme_classic() +
  theme(
    text = element_text(
      family = "Arial"
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(
      color = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    ),
    legend.position = "none"
  )

pC

