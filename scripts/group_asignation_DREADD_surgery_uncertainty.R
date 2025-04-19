#This script aim to explore adiposity index (fat/lean mass) in orexin-cre and wt mice
#we aim split the orexin cre males into four groups:
#n=4 Males orexin - cre with inhibitory DREADD - uncertainty
#n=4 Males WT with control DREADD - no uncertainty
#n=4 Males orexin - cre with control DREADD - uncertainty
#n=4 Males orexin - cre with control DREAD - no uncertainty

#libraries####
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggrepel)
#Data import####
echomri_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echomri_data <- echomri_data %>% 
  filter(COHORT ==11) %>% 
  select(ID, adiposity_index,SEX,STRAIN) %>% #good so here we checked there is 27 animals in total (cohort 11)
  filter(SEX=="M") %>% 
  filter(ID != c(298, 315)) %>%  # died in surgery
  mutate(GROUP = case_when(
    ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_NO_UNCERT",
    ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_NO_UNCERT",
    ID %in% c(307, 316, 319) ~ "OREXIN_CRE_CONTROL_UNCERT"
    )) %>% 
  drop_na()
# Summarize the data####
summary_data <- echomri_data %>%
  group_by(GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

#plot####
plot <- ggplot() +
  geom_point(data = echomri_data, aes(x = GROUP, y = adiposity_index), alpha = 0.7) +  # Optional: raw data points
  geom_point(data = summary_data, aes(x = GROUP, y = mean_adiposity), 
             color = "red", size = 3) +
  geom_errorbar(data = summary_data, aes(x = GROUP, 
                                         ymin = mean_adiposity - se_adiposity, 
                                         ymax = mean_adiposity + se_adiposity),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("Adiposity Index") +
  xlab("Virus Group")
plot
# Check if adiposity index differs between groups using ANOVA####
model <- lm(adiposity_index ~ GROUP, data = echomri_data)
anova_result <- anova(model)
print(anova_result) #no diff among groups 

#postoperative period####
bw <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import
bw <- bw %>% 
  filter(COHORT ==11) %>% 
  select(ID,BW, DATE,STRAIN,COMMENTS) %>% #good so here we checked there is 27 animals in total (cohort 11)
  filter(ID != c(298, 315)) %>%  # died in surgery
  mutate(GROUP = case_when(
    ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_NO_UNCERT",
    ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_NO_UNCERT",
    ID %in% c(307, 316, 319) ~ "OREXIN_CRE_CONTROL_UNCERT"
  )) %>% 
  filter(!is.na(GROUP))

#plot####
# Create a surgery_events dataframe
surgery_events <- bw %>%
  filter(COMMENTS %in% c("BILATERAL_CONTROL_DREADD_INJECTION_SURGERY", 
                         "BILATERAL_INH_DREADD_INJECTION_SURGERY"))

# Then, update your plot with an additional geom_point layer to mark those dates
bw %>%
  ggplot(aes(x = DATE, y = BW)) + 
  geom_point(alpha = 0.2) +
  geom_line(aes(group = ID)) +
  geom_text(aes(label = ID), size = 2) +
  geom_point(data = surgery_events, aes(x = DATE, y = BW), 
             color = "red", size = 3, shape = 21, fill = "yellow") +  # Highlight the day of the surgery
  labs(
    x = "Date",
    y = "body weight (g)"
  ) +
  facet_wrap(~STRAIN) +
  theme_classic()

