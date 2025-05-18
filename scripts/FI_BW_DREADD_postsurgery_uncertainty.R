#This script aim to explore BW changes after virus injection in orexin-cre and wt mice
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr) #to use drop_na()
library(lme4)
#BW####
BW_data <-read_csv("~/Documents/GitHub/data/data/BW.csv") #data import
  
BW_data <- BW_data %>% 
    filter(COHORT == 11) %>% 
    filter(SEX == "M") %>% #We just want to follow the animals with surgery 
    filter(!ID %in% c(298, 315)) %>% #died in surgery
    filter(!ID %in% c(308, 317, 324)) %>% #C57BL6J not used in surgery
    mutate(GROUP = case_when(
      ID %in% c(297, 313,318, 320) ~ "OREXIN_CRE_DREADD",
      ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL",
      ID %in% c(307, 316, 319,321,323,325) ~ "OREXIN_CRE_CONTROL"
     )) 
  # %>% 
  # ungroup() %>% 
  # mutate(GROUP = case_when(
  #   ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_UNCERT",
  #   ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
  #   ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_NO_UNCERT",
  #   ID %in% c(307,316, 319) ~ "OREXIN_CRE_CONTROL_NO_UNCERT"
  # )) %>% 
  # filter(ID %in% c(297, 305, 306, 307, 314, 316, 321,325)) #only calorimetry guys

BW_data_1 <- BW_data %>% 
      group_by(ID) %>% 
    arrange(DATE) %>% 
     mutate(bw_rel = 100 * (BW - first(BW)) / first(BW),
            body_lag = lag(BW) - BW) %>% 
    group_by(ID) %>% 
    mutate(day_rel = DATE - first(DATE))
  

n_distinct(BW_data$ID) #here we know there is 14 animals
  
# Subset rows where COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE"
highlight_data <- BW_data %>%
  filter(COMMENTS == "FIRST_DAY_JUST_FED3_BASELINE")

plot <- BW_data_1 %>% 
  ggplot(aes(DATE, BW, group = ID, color = SEX)) +
  geom_point() +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = BW_data_1,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE) +
  # Highlight using vertical lines
  geom_vline(data = highlight_data, 
             aes(xintercept = as.numeric(DATE)), 
             linetype = "dashed", color = "red", alpha = 0.7)
plot

#echoMRI data####
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data <- echoMRI_data %>% 
  filter(COHORT == 11) %>% 
  filter(SEX == "M") %>% 
  filter(!ID %in% c(298, 315)) %>% #died in surgery
  filter(!ID %in% c(308, 317, 324)) %>% #C57BL6J not used in surgery
  mutate(GROUP = case_when(
    ID %in% c(297, 313,318, 320) ~ "OREXIN_CRE_DREADD",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL",
    ID %in% c(307, 316, 319,321,323,325) ~ "OREXIN_CRE_CONTROL"
  )) %>% 
  select(ID,Date,Fat,Lean,Weight,n_measurement,adiposity_index,GROUP) 

plot_echo <- echoMRI_data %>% 
  ggplot(aes(Date, adiposity_index, group = ID)) +
  geom_line() +
  facet_wrap(~GROUP) +
  geom_text(data = echoMRI_data,
            aes(label = ID), 
            hjust = -0.2, vjust = 0.5, 
            size = 3, show.legend = FALSE)
plot_echo

# body comp analysis
echoMRI_data <-read_csv("~/Documents/GitHub/data/data/echomri.csv") #data import

echoMRI_data_analysis <- echoMRI_data %>% 
  filter(COHORT == 11) %>% 
  filter(SEX == "M") %>% 
  filter(!ID %in% c(298, 315)) %>% #died in surgery
  filter(!ID %in% c(308, 317, 324)) %>% #C57BL6J not used in surgery
  mutate(GROUP = case_when(
    ID %in% c(321, 323, 325) ~ "OREXIN_CRE_CONTROL_UNCERT",
    ID %in% c(297, 313, 318, 320) ~ "OREXIN_CRE_DREADD_UNCERT",
    ID %in% c(305, 306, 314, 322) ~ "WT_CONTROL_NO_UNCERT",
    ID %in% c(307,316, 319) ~ "OREXIN_CRE_CONTROL_NO_UNCERT"
  )) %>%
  filter(n_measurement == 2)  #only cares what happen FEDs experience 

# Summarize the data####
summary_data <- echoMRI_data_analysis %>%
  group_by(GROUP) %>%
  summarise(
    mean_adiposity = mean(adiposity_index, na.rm = TRUE),
    se_adiposity = sd(adiposity_index, na.rm = TRUE) / sqrt(n()),
    mean_Fat = mean(Fat, na.rm = TRUE),
    se_Fat = sd(Fat, na.rm = TRUE) / sqrt(n()),
    mean_Lean = mean(Lean, na.rm = TRUE),
    se_Lean = sd(Lean, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) 

#plot####
plot <- ggplot() +
  geom_point(data = echoMRI_data_analysis, aes(x = GROUP, y = Lean), alpha = 0.7) +
  geom_text(data = echoMRI_data_analysis, 
            aes(x = GROUP, y = Lean, label = ID), 
            vjust = -0.5, size = 3) +  # Adjust vjust/size as needed
  geom_point(data = summary_data, aes(x = GROUP, y = mean_Lean), 
             color = "red", size = 3) +
  geom_errorbar(data = summary_data, aes(x = GROUP, 
                                         ymin = mean_Lean - se_Lean, 
                                         ymax = mean_Lean + se_Lean),
                width = 0.2, color = "red") +
  theme_minimal() +
  ylab("Lean (g)") +
  xlab("Virus Group")
plot

# Check if adiposity index differs between groups using ANOVA####
model <- lm(Lean ~ GROUP, data = echoMRI_data_analysis)
anova_result <- anova(model)
print(anova_result) #no diff among groups 

#Now this model is just for lean mass but I replaced those variables for "Fat" and "Adiposity_index" and
#also NO DIFFERENCES AMONG GROUPS
#Osea a todos los animales las cirugias y los FEDs les pegan por igual en relacion al body comp. 

