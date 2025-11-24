#COHORT 6 24h LOCOMOTION
#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels
library(lme4)
library(stringr)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#locomotion####
# build the summarized dataset
sable_loc_data <- sable_dwn %>% 
  filter(COHORT == 6) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
         SEX = case_when(
           startsWith(as.character(ID), "1") ~ "M",
           startsWith(as.character(ID), "2") ~ "F",
           TRUE ~ NA_character_
         )) %>% 
  mutate(DRUG= case_when(
    ID==1001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==1001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==1004 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==1004 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1005 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1005 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1006 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1006 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==1006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==1007 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==1007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1007 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==2001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==2002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==2002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
    ID==2003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==2003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2003 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==2004 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2004 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2004 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==2005 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2005 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_donated",
    ID==2005 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2006 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==2006 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==2006 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_medchem",
    ID==2007 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==2007 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==2007 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
  )) %>% 
  filter(grepl("AllMeters_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,SEX) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX,DRUG) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # ✅ Keep only the first 4 hours of the dark period
  filter(lights == "off", zt_time < 9) %>% 
  
  
  # calculate LOCOMOTION for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, DRUG, lights,SEX) %>% 
  mutate(loc_act = value - lag(value)) %>%
  filter(loc_act >= 0) %>%
  summarise(total_act = sum(loc_act), .groups = "drop") %>%
  
  # keep both complete days
  filter(!ID %in% c(1003,2001), is_complete_day == 1, complete_days %in% c(1,2)) %>% #ID 1003 DIED DURING EXPERIMENT AND ID 2005 WAS IN CAGE 5
  
  # average across the 2 days per ID × DRUG × lights
  group_by(ID, DRUG,lights,SEX) %>% 
  summarise(total_act = mean(total_act), .groups = "drop") %>%
  mutate(
    DRUG = factor(DRUG, 
                  levels = c("vehicle", 
                             "RTIOXA_43_medchem", 
                             "RTIOXA_43_donated")))


sable_loc_data %>%
  group_by(DRUG,SEX) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

# ---- 1. Ensure factors are correct ----
sable_loc_data <- sable_loc_data %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_medchem", "RTIOXA_43_donated")),
    SEX = factor(SEX)#,
  #  lights = factor(lights)
  )

# ---- Plot mean ± SEM bars with individual points ----
ggplot(sable_loc_data, aes(x = DRUG, y = total_act, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.1)) +
  geom_text(aes(label = ID),
            position = position_jitter(width = 0.1, height = 0),
            vjust = -0.8, size = 3, color = "black") +
  facet_grid( ~ SEX) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",
    "RTIOXA_43_medchem" = "darkgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "locomotion (m/period)", x = "")



# --- Ensure factors are correct ---
sable_loc_data <- sable_loc_data %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_medchem", "RTIOXA_43_donated")),
    SEX = factor(SEX, levels = c("F","M")),
    lights = factor(lights, levels = c("off","on"))
  )

# --- Function to fit LMM, get emmeans, pairwise vs vehicle ---
get_emm_contrasts <- function(data_sex) {
  
  # Fit model
  model <- lmer(total_act ~ DRUG + (1 | ID), data = data_sex)
  
  # Estimated marginal means
  emm <- emmeans(model, ~ DRUG)
  
  # Pairwise comparisons vs vehicle
  posthoc <- pairs(emm, adjust = "bonferroni")
  posthoc_df <- as.data.frame(posthoc) %>%
    mutate(
      group1 = sub(" - .*", "", contrast),
      group2 = sub(".*- ", "", contrast),
      label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  list(emm = emm, posthoc = posthoc_df)
}

# --- Split by sex ---
sable_loc_F <- subset(sable_loc_data, SEX == "F")
sable_loc_M <- subset(sable_loc_data, SEX == "M")

res_F <- get_emm_contrasts(sable_loc_F)
res_M <- get_emm_contrasts(sable_loc_M)

# --- Combine results for plotting ---
# Add SEX column to posthoc
res_F$posthoc$SEX <- "F"
res_M$posthoc$SEX <- "M"
posthoc_all <- bind_rows(res_F$posthoc, res_M$posthoc)

# --- Plot mean ± SEM bars with individual points and significance ---
ggplot(sable_loc_data, aes(x = DRUG, y = total_act, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.1)) +
  facet_grid(SEX ~ lights) +
  scale_fill_manual(values = c(
    "vehicle" = "white",
    "RTIOXA_43_donated" = "#66C2A5",    # light green
    "RTIOXA_43_medchem" = "darkgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    y = "Locomotion (m/period)",
    x = ""
  ) +
  # Add significance labels above bars
  geom_text(
    data = posthoc_all %>% filter(group1 == "vehicle"), 
    aes(
      x = group2, 
      y = max(sable_loc_data$total_act) * 1.05, 
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

