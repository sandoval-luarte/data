# This script aims to explore changes in body weight in middle age C57 after different stages of feeding:
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
  filter(COHORT == 2) %>% # Just C57 males and females
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ),
    day_rel = DATE - first(DATE)
  ) %>%mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      DATE == as.Date("2025-06-04") ~ "BW maintenance", 
      COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      day_rel == 234 ~ "BW loss",
     day_rel == 190 ~ "peak obesity",
      TRUE ~ NA_character_
    ) ) %>%
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
  "ad lib_vehicle" = "#E67E22",      # darker orange
  "ad lib_RTIOXA_47" = "#F39C12",    # lighter orange
  "restricted_vehicle" = "#3498DB",  # darker blue
  "restricted_RTIOXA_47" = "#5DADE2" # lighter blue
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
  facet_wrap(~SEX*DIET_FORMULA)+
  format.plot

plot

# Fit model
model <- lmer(BW ~ SEX * STATUS * GROUP * DRUG * DIET_FORMULA + (1|ID), data = BW_data)
summary(model)

# Save emmeans results
emmeans_results <- emmeans(model, pairwise ~ SEX *STATUS * GROUP * DRUG * DIET_FORMULA, adjust = "tukey")

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


