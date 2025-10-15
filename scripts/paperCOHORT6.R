#COHORT 6
#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(ggrepel) # optional, but better for labels

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#TEE####
# build the summarized dataset
sable_TEE_data <- sable_dwn %>% 
  filter(COHORT == 6) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(DRUG= case_when(
    ID==1001 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_medchem",
    ID==1001 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "vehicle",
    ID==1001 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1002 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "vehicle",
    ID==1002 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1002 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3") ~ "RTIOXA_43_donated",
    ID==1003 & sable_idx %in% c("SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6") ~ "RTIOXA_43_medchem",
    ID==1003 & sable_idx %in% c("SABLE_DAY_7","SABLE_DAY_8","SABLE_DAY_9") ~ "vehicle",
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
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day *and lights period*
  group_by(ID, complete_days, is_complete_day, SABLE, lights) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% #3709, 3717, 3718, 3723, 3725 has cage5 issues and 3724 cage 6 is was not registered correctly
  
  # average across the 2 days per ID × SABLE × lights
  group_by(ID, SABLE,lights) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  )

sable_TEE_data <- sable_TEE_data %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )


ggplot(sable_TEE_data, aes(x = SABLE, y = tee, color = GROUP, group = ID)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_text(aes(label = ID), size = 2.5, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", aes(group = GROUP), size = 1.2) +
  facet_wrap(~lights) +
  labs(y = "TEE (kcal/day)", color = "Group") +
  theme_minimal()


#--without lights division

sable_TEE_data <- sable_dwn %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%   
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                     "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                     "SABLE_DAY_7") ~ "baseline",
    sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
    sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain"
  )) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% 
  ungroup() %>% 
  group_by(ID, SABLE) %>% 
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr!=lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time==0 & is_zt_init == 1,1,0))
  ) %>% 
  ungroup() %>% 
  group_by(ID, complete_days) %>% 
  mutate(is_complete_day = if_else(min(zt_time)==0 & max(zt_time)==23, 1, 0)) %>% 
  ungroup() %>% 
  
  # calculate TEE for each day (no lights split!)
  group_by(ID, complete_days, is_complete_day, SABLE) %>% 
  summarise(tee = sum(value)*(1/60), .groups="drop") %>% 
  
  # keep both complete days
  filter(!ID %in% c(3715,3712), is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!ID %in% c(3709, 3717, 3718, 3723, 3724, 3725)) %>% 
  
  # average across the 2 days per ID × SABLE
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729) ~ "RTIOXA_47"
    )
  ) %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain"))
  )
ggplot(sable_TEE_data, aes(x = SABLE, y = tee, color = DRUG, group = ID)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 2, alpha = 0.5) +
  geom_text(aes(label = ID), size = 2.5, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", aes(group = DRUG), size = 1.2) +
  facet_wrap(~GROUP) +  # only group, no lights
  labs(y = "TEE (kcal/day)", color = "Drug") +
  theme_minimal()

#--colapsing data

# Make a plotting dataset for TEE
TEE_plotdata <- sable_TEE_data %>%
  mutate(
    PlotGroup = case_when(
      SABLE %in% c("baseline", "peak obesity") ~ "all",          # collapse all
      SABLE %in% c("BW loss", "BW maintenance") ~ GROUP,         # separate by GROUP
      SABLE == "BW regain" ~ paste(GROUP, DRUG, sep = "_")       # GROUP × DRUG
    )
  )

# Define custom colors (same as BW for consistency)
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#E67E22",              # orange
  "restricted" = "#3498DB",          # sky blue
  "ad lib_vehicle" = "#E67E22",      # darker orange
  "ad lib_RTIOXA_47" = "#F39C12",    # lighter orange
  "restricted_vehicle" = "#3498DB",  # darker blue
  "restricted_RTIOXA_47" = "#5DADE2" # lighter blue
)

# Plot
TEE_plot <- TEE_plotdata %>%
  ggplot(aes(x = SABLE, y = tee, fill = PlotGroup)) +
  
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7, alpha = 0.7) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  
  geom_point(aes(color = PlotGroup),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  theme_minimal() +
  labs(y = "TEE (kcal/day)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

TEE_plot
# Fit model for TEE
TEE_model <- lmer(tee ~ SABLE * GROUP * DRUG + (1|ID), data = sable_TEE_data)
summary(TEE_model)

# Estimated marginal means + pairwise contrasts
TEE_emm <- emmeans(TEE_model, pairwise ~ SABLE * GROUP * DRUG, adjust = "tukey")

# Convert to data frames
df_emm_TEE <- as.data.frame(TEE_emm$emmeans)     # estimated means
df_pairs_TEE <- as.data.frame(TEE_emm$contrasts) # pairwise comparisons

# Print results
print(df_emm_TEE, n = Inf)
print(df_pairs_TEE, n = Inf)

# Keep only significant contrasts
df_sig_TEE <- df_pairs_TEE %>%
  filter(p.value <= 0.05)

# View significant results
print(df_sig_TEE, n = Inf)
View(df_sig_TEE)
