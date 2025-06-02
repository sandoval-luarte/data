# We aim to test if BW is statistically different after 100 days ( ~12 weeks) with LFD in NZO female mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(ggpubr)

#format plot

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

##bw discrete ####
##body weight baseline and peak of obesity####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID== 3715) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(day_rel = DATE - first(DATE)) %>% 
  filter(day_rel < 100) %>% 
  filter(day_rel %in% c(0, 98))

scale_x_discrete(labels = function(x) {
  x <- as.character(x)
  x[x == "0"] <- "baseline"
  x[x == "98"] <- "peak of obesity"
  return(x)
})

# Reshape to wide format for easier comparison
BW_wide <- BW_data %>%
  select(ID, day_rel, BW) %>%
  pivot_wider(names_from = day_rel, values_from = BW)

# Run paired t-test
t.test(BW_wide$`0`, BW_wide$`98`, paired = TRUE)

plot_2 <- BW_data %>%
  ggplot(aes(x = as.factor(day_rel), y = BW)) +
  geom_point(aes(group = ID), alpha = 0.1) +
  geom_line(aes(group = ID), alpha = 0.1) +
  stat_summary(fun.data = mean_se, 
               geom = "ribbon", fill = "grey70", color = NA, 
               aes(group = 1), alpha = 0.7) +
  stat_summary(fun = mean, 
               geom = "line", color = "black", 
               aes(group = 1)) +
  labs(
    x = "Weeks",
    y = "Body weight (g)"
  ) +
  format.plot +
  scale_x_discrete(labels = function(x) {
    x <- as.character(x)
    x[x == "0"] <- "0"
    x[x == "98"] <- "12"
    return(x)
  })
plot_2


# Create data frame for annotation
pvalue_df <- data.frame(
  group1 = "0",
  group2 = "98",  # or "12" if you've relabeled already
  y.position = 52,  # Adjust based on your plot height
  p = 0.0001,       # Or the actual p-value you got
  label = "***"
)

plot_2 <- plot_2 +
  stat_pvalue_manual(
    pvalue_df,
    label = "label",
    tip.length = 0.01
  )

