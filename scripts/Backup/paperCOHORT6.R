#COHORT 6 TEE
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

#TEE####
# build the summarized dataset
sable_TEE_data <- sable_dwn %>% 
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
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, DRUG,SEX) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days,SEX) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, DRUG, lights,SEX) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(1003,2001), is_complete_day == 1, complete_days %in% c(1,2)) %>% #ID 1003 DIED DURING EXPERIMENT AND ID 2005 WAS IN CAGE 5

  # average across the 2 days per ID × DRUG × lights
  group_by(ID, DRUG,lights,SEX) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>% 
  mutate(
    DRUG = factor(DRUG, 
                   levels = c("vehicle", 
                              "RTIOXA_43_medchem", 
                              "RTIOXA_43_donated")))


sable_TEE_data %>%
  group_by(DRUG,SEX) %>%
  summarise(n_ID = n_distinct(ID)) #this is good

# ---- 1. Ensure factors are correct ----
sable_TEE_data <- sable_TEE_data %>%
  mutate(
    ID = factor(ID),
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_43_medchem", "RTIOXA_43_donated")),
    SEX = factor(SEX),
    lights = factor(lights)
  )

# ---- 2. Fit the linear mixed model ----
# Each animal is a random intercept (repeated measures)
model <- lmer(tee ~ DRUG * SEX * lights + (1 | ID), data = sable_TEE_data)
anova(model, type = 3)  # Type III ANOVA table

# ---- 3. Compute emmeans and pairwise comparisons ----
emm <- emmeans(model, ~ DRUG | SEX * lights)          # marginal means per facet


# Convert emmeans object to data frame
emm_df <- as.data.frame(emm)

# The columns in emm_df include: DRUG, SEX, lights, emmean, SE, df, lower.CL, upper.CL
# Now do pairwise contrasts
posthoc_df <- pairs(emm, adjust = "bonferroni") %>% as.data.frame()

# Add the facet variables (SEX and lights) by matching the 'emm' object
# The 'emmGrid' labels are stored in rownames, but easier is to use the
# original emm_df to merge by DRUG within each SEX*lights

# First, create a helper key in both
posthoc_df <- posthoc_df %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast)
  )

# ---- Plot mean ± SEM bars with individual points ----
ggplot(sable_TEE_data, aes(x = DRUG, y = tee, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(aes(group = ID), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.1)) +
  facet_grid(~ lights*SEX) +
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
    y = "TEE (kcal/24h)",
    x = ""
  )


ggplot(sable_TEE_data, aes(x = DRUG, y = tee, fill = DRUG)) +
  
  # Bars: mean ± SEM
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  
  # Lines connecting the same mouse
  geom_line(aes(group = ID), color = "gray50", alpha = 0.7) +
  
  # Points for each mouse
  geom_point(aes(group = ID), alpha = 0.7, size = 2,
             position = position_jitter(width = 0.1)) +
  
  facet_grid(~ lights*SEX) +
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
    y = "TEE (kcal/24h)",
    x = ""
  )



