# This script aims to explore changes fat mass in middle age NZO after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

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
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      Date == as.Date("2025-02-20") ~ "peak obesity",
      Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      Date == as.Date("2025-05-27") ~ "BW maintenance",
      Date %in% as.Date(c("2025-07-22", "2025-07-21","2025-07-17","2025-07-16",
                          "2025-07-14","2025-07-09","2025-07-08")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS)) %>%  # <-- optional
  filter(!(ID == 3726 & Date == as.Date("2025-04-28")))  #repeated

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("baseline", "peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))
#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14),
  
  # remove background grid lines only
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  # keep axis lines
  axis.line = element_line(color = "black")
)

plot <- echoMRI_data %>%
  ggplot(aes(x = STATUS, y = Fat, fill = GROUP)) +
  
  # mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    color = "black", width = 0.7, alpha = 0.7
  ) +
  
  # error bars (mean ± SE)
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3
  ) +
  
  # individual data points
  geom_point(
    aes(color = DRUG),
    position = position_dodge(width = 0.8),
    alpha = 0.7, size = 2
  ) +
  
  scaleFill +
  theme_minimal() +
  labs(y = "Fat in grams", fill = "Group", color = "Drug") +
  format.plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot

fat_plotdata <-echoMRI_data  %>%
  mutate(
    PlotGroup = case_when(
      STATUS %in% c("baseline", "peak obesity") ~ "all",          # collapse all
      STATUS %in% c("BW loss", "BW maintenance") ~ GROUP,         # separate by GROUP
      STATUS == "BW regain" ~ paste(GROUP, DRUG, sep = "_")       # GROUP × DRUG
    )
  )

# Define custom colors
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#E67E22",              # orange
  "restricted" = "#3498DB",          # sky blue
  "ad lib_vehicle" = "#E67E22",      # darker orange
  "ad lib_RTIOXA_47" = "#F39C12",    # lighter orange
  "restricted_vehicle" = "#3498DB",  # darker blue
  "restricted_RTIOXA_47" = "#5DADE2" # lighter blue
)

plot <- fat_plotdata %>%
  ggplot(aes(x = STATUS, y = Fat, fill = PlotGroup)) +
  
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7, alpha = 0.7) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  
  geom_point(aes(color = PlotGroup),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  theme_minimal() +
  labs(y = "Fat (g)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  format.plot

plot

# Fit model
model <- lmer(Fat ~ STATUS * GROUP * DRUG + (1|ID), data = echoMRI_data)
summary(model)

# Save emmeans results
emmeans_results <- emmeans(model, pairwise ~ STATUS * GROUP * DRUG, adjust = "tukey")
emmeans_results
#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm <- as.data.frame(emmeans_results$emmeans)   # estimated means
df_pairs <- as.data.frame(emmeans_results$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm, n = Inf)
print(df_pairs, n = Inf)
View(df_pairs)   # opens spreadsheet-style viewer in RStudio
