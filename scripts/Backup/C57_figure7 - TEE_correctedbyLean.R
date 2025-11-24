#we aim to evaluate changes in 24h TEE in middle age C57 after different stages of feeding:
#1 from baseline to peak obesity,
#2:from peak of obesity to acute body weight loss
#3 from acute body weight loss to body weight maintenance
#4 from body weight maintenance to body weight gain after RTIOXA-47 injections

#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)
library(ggpubr)

#functions####
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr-20, hr+4))
}

sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds") 

#TEE####
sable_TEE_data <- sable_dwn %>% # Load the data
  filter(COHORT ==2) %>% # JustC57
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  mutate(SABLE= case_when(
    sable_idx %in% c("SABLE_DAY_1",
                     "SABLE_DAY_2",
                     "SABLE_DAY_3",
                     "SABLE_DAY_4") ~ "peak obesity",
    sable_idx %in% c("SABLE_DAY_5",
                     "SABLE_DAY_6",
                     "SABLE_DAY_7",
                     "SABLE_DAY_8") ~ "BW loss",
    sable_idx %in% c("SABLE_DAY_9",
                     "SABLE_DAY_10",
                     "SABLE_DAY_11",
                     "SABLE_DAY_12") ~ "BW maintenance", 
    sable_idx %in% c("SABLE_DAY_13",
                     "SABLE_DAY_14",
                     "SABLE_DAY_15",
                     "SABLE_DAY_16") ~ "BW regain"
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
  filter(is_complete_day == 1, complete_days %in% c(1,2)) %>% 
  filter(!(ID %in% c(7866, 7874, 7877, 7879, 7864, 7881))) %>%  #cage 5 in at least one SABLE stage
  filter(!(ID %in% c(7865, 7875, 7882 ))) %>%  #cage 6 issues in SABLE stage BW Mainten or regain
  
  # average across the 2 days per ID × SABLE
  group_by(ID, SABLE) %>% 
  summarise(tee = mean(tee), .groups = "drop") %>% 
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )) %>%
  mutate(
    SABLE = factor(SABLE, 
                   levels = c("baseline", 
                              "peak obesity", 
                              "BW loss", 
                              "BW maintenance", 
                              "BW regain")))


echoMRI_data <- read_csv("~/Documents/GitHub/data/data/echomri.csv") %>%
  filter(COHORT == 2) %>% # Just C57 males and females
  mutate(ID = as.factor(ID)) %>% 
  group_by(ID) %>%
  arrange(Date) %>%
  mutate(
    GROUP = case_when(
      ID %in% c(7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ))%>%
  select(ID, Date, Fat, Lean, Weight, n_measurement, adiposity_index, GROUP, DRUG,SEX, DIET_FORMULA) %>%
  mutate(
    day_rel = Date - first(Date),
    STATUS = case_when(
      Date == as.Date("2025-03-07") ~ "peak obesity",
      Date == as.Date("2025-04-21") ~ "BW loss",
      Date == as.Date("2025-06-05") ~ "BW maintenance",
      Date %in% as.Date(c("2025-09-11", "2025-09-10","2025-09-05","2025-09-04",
                          "2025-09-02","2025-09-01","2025-08-28","2025-08-27")) ~ "BW regain",
      TRUE ~ NA_character_
    )) %>% 
  filter(!is.na(STATUS))

# Make STATUS an ordered factor
echoMRI_data <- echoMRI_data %>%
  mutate(STATUS = factor(STATUS, 
                         levels = c("peak obesity", "BW loss", 
                                    "BW maintenance", "BW regain")))

# Rename STATUS to SABLE for merging
echoMRI_data <- echoMRI_data %>%
  rename(SABLE = STATUS)

# Left join lean mass info into TEE dataset
sable_TEE_adj <- sable_TEE_data %>%
  left_join(
    echoMRI_data %>% select(ID, SABLE, Lean,SEX,DIET_FORMULA),
    by = c("ID", "SABLE")
  ) 

#statistical model
model_TEE_lean <- lmer(tee ~ SABLE * GROUP *SEX + Lean + (1 | ID), data = sable_TEE_adj)
summary(model_TEE_lean)

sable_TEE_adj %>%
  group_by(SEX) %>%
  summarise(n_ID = n_distinct(ID))


emm_TEE <- emmeans(model_TEE_lean, ~ SABLE * GROUP*SEX, cov.reduce = mean)
emm_TEE_df <- as.data.frame(emm_TEE)

# Pairwise contrasts within each GROUP
contrasts_by_group <- contrast(emm_TEE, method = "pairwise", by = "GROUP")

# Convert to a data frame
contrasts_df <- as.data.frame(contrasts_by_group)


ggplot(emm_TEE_df, aes(x = SABLE, y = emmean, color = GROUP, group = GROUP)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  geom_line(position = position_dodge(0.2), linewidth = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.1, position = position_dodge(0.2)) +
  facet_wrap(~SEX) +
  theme_minimal(base_size = 14) +
  labs(y = "TEE (adjusted for Lean mass)", x = "SABLE phase",
       color = "Group",
       title = "TEE across SABLE phases (adjusted for Lean mass)") +
  theme(legend.position = "top")

sable_TEE_adj %>%
  group_by(SEX, SABLE) %>%
  summarise(mean_lean = mean(Lean), .groups = "drop")

