#This script aims to evaluate changes in BW, FI and body comp of animals after 5 days of IP injections of RTIOXA 43 or RTIOXA 47
#libraries
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) 
library(tidyr)  # to use drop-na()


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


#body weight####
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(9, 12)) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    GROUP = case_when(
      ID %in% c(8100, 8103, 8095) ~ "VEHICLE",
      ID %in% c(8079, 8076, 8074) ~ "VEHICLE_2",
      ID %in% c(8078, 8077, 8075) ~ "RTI_47_10mgkg",
      ID %in% c(8099, 8096, 8102) ~ "RTI_43_10mgkg_Y",
      ID %in% c(8098, 8101, 8097) ~ "RTI_43_10mgkg_MED"
    )
  ) %>% 
  group_by(ID) %>% 
  mutate(day_rel = DATE - first(DATE))

n_distinct(BW_data$ID) #here we know there is 22 animals


plot <- BW_data %>% 
  ggplot(aes(day_rel, bw_rel, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP, scales = "free_x") +  # free x-axis scales
  geom_smooth() +
  geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + 
  labs(
    x = "Date",
    y = "BW (grams)"
  )  + 
  format.plot

deposit_summary <- read_csv("../data/WHITE_DEPOSIT.csv") %>% 
  group_by(DRUG) 

deposit_summary <-deposit_summary %>% 
  summarise(
    n = n(),
    mean = mean(DEPOSIT_AMOUNT_MG, na.rm = TRUE),
    sd = sd(DEPOSIT_AMOUNT_MG, na.rm = TRUE)
  ) %>% 
  filter(!is.na(mean))
  

# Plot
ggplot(deposit_summary, aes(x = DRUG, y = mean)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(
    x = "Drug Group",
    y = "Mean Deposit (mg)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  format.plot

