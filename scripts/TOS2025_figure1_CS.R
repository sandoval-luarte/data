#changes body weight in c57 and NZO mice in baseline, peak obesity, bw loss, 
#BW maintenance and BW regain
#Figure 1 poster TOS CS

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(tidyverse)
library(lme4)  # For mixed-effects models
library(ggpubr)  
library(lmerTest)
library(emmeans)

BW_data <- read_csv("../data/BW.csv") %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  filter(SEX == "F") %>% # Just females from both strains
  filter(!ID %in% c(3712, 3715)) %>% # died during study
  filter(!ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726)) %>% # ad lib NZO
  filter(!ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                    7882, 7883)) %>% # ad lib C57
  filter(!ID %in% c(3723, 3724, 3725)) %>% # CAGE 5 ISSUES NZO AND 3724 CAGE 6 ISSUES 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    body_lag = (lag(BW) - BW),
    day_rel = DATE - first(DATE)) %>% 
mutate(
    STATUS = case_when(
      day_rel == 0 ~ "baseline", 
      STRAIN == "NZO/HlLtJ" & DATE == as.Date("2025-02-21") ~ "peak obesity",
      STRAIN == "C57BL6/J" & DATE == as.Date("2025-02-26") ~ "peak obesity",
      STRAIN == "C57BL6/J" & day_rel == 234 ~ "BW loss",
      STRAIN == "NZO/HlLtJ" & day_rel == 161 ~ "BW loss",
      STRAIN == "C57BL6/J" & day_rel == 300 ~ "BW maintenance", 
      STRAIN == "NZO/HlLtJ" & day_rel == 196 ~ "BW maintenance",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    ) ) 

n_distinct(BW_data$ID) #9 NZO in total and 4 C57 in total

highlight_data <- BW_data %>%
  filter(STATUS == "peak obesity")

highlight_data_2 <- BW_data %>%
  filter(STATUS == "BW loss")

highlight_data_3 <- BW_data %>%
  filter(STATUS == "BW maintenance")

highlight_data_4 <- BW_data %>%
  filter(STATUS == "BW regain")

#format plot

scaleFill <- scale_fill_manual(values = c("#C03830FF", "#317EC2FF"))


format.plot <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),          
  panel.spacing.y = unit(1.5, "lines"),  
  axis.text = element_text(family = "Helvetica", size = 13),
  axis.title = element_text(family = "Helvetica", size = 14))


plot <- BW_data %>%
  ggplot(aes(day_rel, BW, group = ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~STRAIN) +
  geom_smooth()+
  #geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.6) + # ID label
  labs(
    x = "Days relative to baseline",
    y = "BW (grams)") +
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

# get the day of peak obesity for each mouse
peak_days <- BW_data %>%
  filter(STATUS == "peak obesity") %>%
  select(ID, peak_day = day_rel)

# join back to original data and keep baseline → peak obesity
mdl_data_bw <- BW_data %>%
  left_join(peak_days, by = "ID") %>%
  filter(day_rel >= 0 & day_rel <= peak_day)   # keep only baseline → peak obesity

mdl_data_bw %>%
  group_by(ID) %>%
  summarise(min_day = min(day_rel), max_day = max(day_rel), .groups = "drop")

mdl_data_bw <- mdl_data_bw %>% 
  group_by(ID) %>% 
  mutate(BW_rel=((BW - min(BW))/min(BW)*100)) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  select(
    ID, BW, time, BW_rel, DATE
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols="BW") %>% 
  select(
    ID, time, name, value
  ) 

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD for NZO and HFD for C57)
lmer_bw <- lmerTest::lmer(
  data = mdl_data_bw,
  value ~ time + (1+time|ID)
)
summary(lmer_bw)
view(coef(lmer_bw))

coef_bw_peakobesity <- coef(lmer_bw)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_bw= `(Intercept)`)

write_csv(coef_bw_peakobesity, "coef_bw_peakobesity.csv")  # Save the dataframe as an csv file


