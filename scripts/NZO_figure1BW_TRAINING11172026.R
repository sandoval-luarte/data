# This script aims to explore changes in body weight in middle age NZO after different stages of feeding:
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

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% # Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    ),
    day_rel = DATE - first(DATE)
  ) %>%
  mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      DATE == as.Date("2025-02-24") ~ "peak obesity",
      day_rel == 161 ~ "BW loss",
      COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance", 
      COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(STATUS)) 


# Make STATUS an ordered factor
BW_data <- BW_data %>%
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

plot <- BW_data %>%
  ggplot(aes(x = STATUS, y = BW, fill = GROUP)) +
  
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
  labs(y = "Body Weight (BW) in grams", fill = "Group", color = "Drug") +
  format.plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot

BW_plotdata <- BW_data %>%
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
  "ad lib_vehicle" = "red",      # darker orange
  "ad lib_RTIOXA_47" = "green",    # lighter orange
  "restricted_vehicle" = "pink",  # darker blue
  "restricted_RTIOXA_47" = "brown" # lighter blue
)

plot <- BW_plotdata %>%
  ggplot(aes(x = STATUS, y = BW, fill = PlotGroup)) +
  
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
  labs(y = "Body Weight (g)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  format.plot

plot

# Fit model
#fits a linear mixed-effects model to evaluate how body weight (BW) changes 
#across different experimental stages (STATUS), between experimental groups (GROUP), 
#and depending on treatment (DRUG)

model <- lmer(BW ~ STATUS * GROUP * DRUG + (1|ID), data = BW_data) # 1|ID is important here because you measure body weight repeatedly in the same animals, 
#observations from the same animal are not independent, the random intercept allows each animal to have its own baseline body weight.
summary(model)

#The asterisks tell R to include:
#1. Main effects: STATUS, GROUP and DRUG
#2. Two-way interactions: STATUS × GROUP, STATUS × DRUG, GROUP × DRUG
#3. Three-way interaction: STATUS × GROUP × DRUG

#In other words, you are testing whether body weight changes across experimental stages
#and whether those changes differ depending on the experimental group and drug treatment.

#Your model is:model <- lmer(BW ~ STATUS * GROUP * DRUG + (1|ID), data = BW_data), You have: 110 BW observations, 22 individual animals,
#Five STATUS categories (one reference category plus four displayed in the output), Two GROUP categories, Two DRUG categories.
#The model considers all combinations of STATUS, GROUP, and DRUG.

#The random effect section Tdescribes the variability in your data.
#Component	    Variance	    Meaning: 
#ID (intercept)	9.997	        Differences in average BW between animals.
#Residual    	5.996	        Variation in BW not explained by the model, including within-animal variation.

#For example, some animals might consistently weigh more than others, even when they are at the same experimental stage.
#The model accounts for these differences rather than treating all 110 observations as independent.
#You have 22 animals contributing to those 110 observations.

#3. Understanding fixed effects, This is the most important section of your output.
#Before interpreting individual coefficients, we need to understand something called the reference category.
#R uses one category from each categorical variable as the reference against which the other categories are compared.
#Based on your output, your reference categories appear to be: STATUS: Baseline (the omitted category),  GROUP: Unrestricted (the omitted category) and DRUG: RTI-47 (the omitted category).
#This is important because all the coefficients in your output are interpreted relative to this reference condition, unless they represent interactions.

# Save emmeans results
emmeans_results <- emmeans(model, pairwise ~ STATUS * GROUP, adjust = "tukey")

#to evaluate baseline ad lib - baseline restricted   p=1
# to evaluate peak obesity ad lib - peak obesity restricted p=0.99

# Convert to data frame
df_emm <- as.data.frame(emmeans_results$emmeans)   # estimated means
df_pairs <- as.data.frame(emmeans_results$contrasts)  # pairwise comparisons

# Print all rows
print(df_emm, n = Inf)
print(df_pairs, n = Inf)
# Keep only significant contrasts
df_sig <- df_pairs %>%
  filter(p.value <= 0.05)

# View the results
print(df_sig, n = Inf)
View(df_sig)

