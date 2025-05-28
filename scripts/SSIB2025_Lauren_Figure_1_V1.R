# We aim to test if BW is statistically different after 100 days ( ~14 weeks) with LFD in NZO female mice

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(ggpubr)

#format plot
scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))

format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),          
        panel.spacing.y = unit(1.5, "lines"),  
        legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 13),
        axis.title = element_text(family = "Helvetica", size = 14))

##bw over time - timecourse####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID== 3715) %>% #died during study
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE)) %>% 
  mutate(relative_weeks = as.integer(day_rel / 7) + 1) %>% 
  filter(relative_weeks <=13)

n_distinct(BW_data$ID) #here we know there is 22 animals

# Custom function for mean ± SEM
mean_se <- function(x) {
  se <- sd(x) / sqrt(length(x))
  return(c(y = mean(x), ymin = mean(x) - se, ymax = mean(x) + se))
}

plot_1 <- BW_data %>%
  ggplot(aes(x = relative_weeks, y = BW)) +
  geom_point(aes(group = ID), alpha = 0.1) +             # individual points
  geom_line(aes(group = ID), alpha = 0.1) +              # individual lines
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
  scale_y_continuous(limits = c(25, 55)) +
  scale_x_continuous(breaks = seq(3, 12, by = 3)) 
plot_1

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

plot_2 <- plot_2 +
  annotate("text", 
           x = 1.5,  # midpoint between x = "0" and "98"
           y = max(BW_data$BW, na.rm = TRUE) + 2,  # a bit above the max BW
           label = "t(22) = -17.64, p < 0.001\nΔ = 14.2 g [95% CI: 12.5–15.9]", 
           size = 4, hjust = 0.5)

plot_2


# Create an arranged plot
ggarrange(plot_1, plot_2, 
                  nrow = 1, 
                  align = "hv", 
                  widths = c(0.65, 0.35),  # Increase width of panel A, decrease B
                  labels = c("A", "B"))
          
