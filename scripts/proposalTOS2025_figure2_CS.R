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
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  filter(SEX == "F") %>% # Just females from both strains
  group_by(ID) %>%
  arrange(Date) %>%
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726)) %>% # ad lib NZO
  filter(!ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                    7882, 7883)) %>% # ad lib C57
  filter(!ID %in% c(3723, 3724, 3725)) %>% # CAGE 5 ISSUES NZO AND 3724 CAGE 6 ISSUES
  filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
  filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index,STRAIN,DIET_FORMULA) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      n_measurement == 1 ~ "baseline",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-03-07") ~ "peak obesity",
      STRAIN == "NZO/HlLtJ" & Date == as.Date("2025-02-20") ~ "peak obesity",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-03-31") ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & Date %in% as.Date(c("2025-04-28", "2025-05-05","2025-05-05","2025-05-06")) ~ "BW loss",
      STRAIN == "C57BL6/J" & Date == as.Date("2025-06-05") ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & Date == as.Date("2025-05-27") ~ "BW maintenance",
      STRAIN == "C57BL6/J" & Date %in% as.Date(c("2025-08-28", #7878 last day of echo
                                                 "2025-09-01", #7872 last day of echo
                                                 "2025-09-02", #7874 last day of echo
                                                 "2025-09-05")#7877 last day of echo
                                               ) ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & Date %in% as.Date(c("2025-07-08", #3708 last day of echo
                                                  "2025-07-09",#3710 last day of echo
                                                  "2025-07-14",#3714 last day of echo
                                                  "2025-07-17", #3720, 3721, 3722 last day of echo
                                                  "2025-07-22"#3727, 3728, 3729 last day of echo
                                                  ) 
                                                ) ~ "BW regain",
      
      TRUE ~ NA_character_
    )
  )
n_distinct(echoMRI_data$ID) #9 NZO and 4 C57 so 13 ID in total

highlight_data <- echoMRI_data %>%
  filter(STATUS == "peak obesity")

highlight_data_2 <- echoMRI_data %>%
  filter(STATUS == "BW loss")

highlight_data_3 <- echoMRI_data %>%
  filter(STATUS == "BW maintenance")

highlight_data_4 <- echoMRI_data %>%
  filter(STATUS == "BW regain")

#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14))


plot <- echoMRI_data %>%
  ggplot(aes(day_rel, Fat, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ID) +
  #geom_smooth()+
  #geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + # ID label
  labs(
    x = "Days relative to baseline",
    y = "Fat mass (grams)") +
  geom_vline(data = highlight_data, 
             aes(xintercept = day_rel), 
             linetype = "dashed", color = "red", alpha = 0.7)+ #red means peak obesity
  geom_vline(data = highlight_data_2, 
             aes(xintercept = day_rel), 
             linetype = "dashed", color = "green", alpha = 0.7)+ # green means end of BW loss
  geom_vline(data = highlight_data_3, 
             aes(xintercept = day_rel), 
             linetype = "dashed", color = "blue", alpha = 0.7)+ #blue means end of BW maintenance
  geom_vline(data = highlight_data_4, 
             aes(xintercept = day_rel), 
             linetype = "dashed", color = "orange", alpha = 0.7)+ #orange means end of BW regain
  format.plot

plot


