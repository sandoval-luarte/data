# This script aims to explore changes in body weight, food intake,
# body composition and locomotion in middle age NZO female mice after 
# different stages of feeding:1: Basal, 2: peak obesity, 3: Acute body 
#weight loss 4: BW maintenance
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()
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

#NZO####
##body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 2 & COHORT < 6) %>% #Just NZO females
  filter(!ID %in% c(3712, 3715)) %>% #died during study %>% 
  filter(DATE <= "2025-02-24") %>% # day 2 of restriction diet
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
         body_lag = (lag(BW) - BW),
         GROUP = case_when(
           ID %in% c(3706, 3707, 3709, 3711, 3712, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
           ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted")) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE)) %>% 
  filter(day_rel <=100)

n_distinct(BW_data$ID) #here we know there is 22 animals

# Custom function for mean ± SEM
mean_se <- function(x) {
  se <- sd(x) / sqrt(length(x))
  return(c(y = mean(x), ymin = mean(x) - se, ymax = mean(x) + se))
}

plot_1 <- BW_data %>%
  ggplot(aes(x = day_rel, y = bw_rel)) +
  geom_point(aes(group = ID), alpha = 0.1) +             # individual points
  geom_line(aes(group = ID), alpha = 0.1) +              # individual lines
  stat_summary(fun.data = mean_se, 
               geom = "ribbon", fill = "grey70", color = NA, 
               aes(group = 1), alpha = 0.7) +
  stat_summary(fun = mean, 
               geom = "line", color = "black", 
               aes(group = 1)) +
  labs(
    x = "Days",
    y = "%Body weight gain"
  ) +
  format.plot +
 # scale_y_continuous(limits = c(25, 60)) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) 
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
    x = "Days",
    y = "Body weight (g)"
  ) +
  format.plot +
  scale_x_discrete(labels = function(x) {
    x <- as.character(x)
    x[x == "0"] <- "0"
    x[x == "98"] <- "100"
    return(x)
  })
plot_2



# Create an arranged plot
ggarrange(plot_1, plot_2, 
                  nrow = 1, 
                  align = "hv", 
                  widths = c(0.65, 0.35),  # Increase width of panel A, decrease B
                  labels = c("A", "B"))
          
