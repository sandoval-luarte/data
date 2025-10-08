# --- Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(ggpubr)

# --- Function ---
zt_time <- function(hr){
  return(if_else(hr >= 20 & hr <= 23, hr - 20, hr + 4))
}

# --- Load data ---
sable_dwn <- readRDS(file = "../data/sable_downsampled_data.rds")

# --- Locomotion data prep ---
sable_locomotion_data <- sable_dwn %>% 
  filter(COHORT > 1 & COHORT < 6) %>% # Just NZO AND C57
  mutate(SABLE = case_when(
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3",
                                             "SABLE_DAY_4","SABLE_DAY_5","SABLE_DAY_6",
                                             "SABLE_DAY_7") ~ "baseline",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_8","SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11") ~ "peak obesity",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_12","SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15") ~ "BW loss", 
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_16","SABLE_DAY_17","SABLE_DAY_18","SABLE_DAY_19") ~ "BW maintenance",
    STRAIN == "NZO/HlLtJ" & sable_idx %in% c("SABLE_DAY_20","SABLE_DAY_21","SABLE_DAY_22","SABLE_DAY_23") ~ "BW regain",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_1","SABLE_DAY_2","SABLE_DAY_3","SABLE_DAY_4") ~ "peak obesity",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_5","SABLE_DAY_6","SABLE_DAY_7","SABLE_DAY_8") ~ "BW loss",
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_9","SABLE_DAY_10","SABLE_DAY_11","SABLE_DAY_12") ~ "BW maintenance", 
    STRAIN == "C57BL6/J" & sable_idx %in% c("SABLE_DAY_13","SABLE_DAY_14","SABLE_DAY_15","SABLE_DAY_16") ~ "BW regain"
  )           
) %>%          
  filter(!is.na(SABLE)) %>%
  
  filter(grepl("AllMeters_*", parameter)) %>%
  group_by(ID, SABLE,STRAIN,SEX,DIET_FORMULA) %>%
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>%
  ungroup() %>%
  group_by(ID, complete_days,STRAIN,SEX,DIET_FORMULA) %>%
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>%
  ungroup() %>%
  
  # compute locomotor activity incrementally
  group_by(ID, complete_days, SABLE,STRAIN, SEX,DIET_FORMULA) %>%
  mutate(loc_act = value - lag(value)) %>%
  filter(loc_act >= 0) %>%
  summarise(total_act = sum(loc_act), .groups = "drop") %>%
  
  # exclude problematic IDs (same as TEE)
  filter(!ID %in% c(3712, 3715)) %>% # NZO mice which died during study
  filter(complete_days %in% c(1, 2)) %>%
  
  # average across 2 complete days
  group_by(ID, SABLE,STRAIN,SEX,DIET_FORMULA) %>%
  summarise(total_act = mean(total_act), .groups = "drop") %>%
  
  # reattach GROUP and DRUG
  mutate(
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 7876, 7879, 7880, 7881,
                7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728, 7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871, 7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729, 7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    )
  ) %>%
  mutate(
    SABLE = factor(SABLE,
                   levels = c("baseline", "peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))
  ) %>% 
    filter(!is.na(SABLE)) %>% 
  filter(SEX=="F") 

write_csv(sable_locomotion_data, "../data/sable_locomotion_data.csv") # Save as CSV


# --- Collapse into PlotGroup for graph consistency ---
loc_plotdata <- sable_locomotion_data %>%
  mutate(
    PlotGroup = case_when(
      SABLE %in% c("baseline", "peak obesity") ~ "all",
      SABLE %in% c("BW loss", "BW maintenance") ~ GROUP,
      SABLE == "BW regain" ~ paste(GROUP, DRUG, sep = "_")
    )
  )

# --- Custom color palette (same as TEE) ---
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#E67E22",
  "restricted" = "#3498DB",
  "ad lib_vehicle" = "#E67E22",
  "ad lib_RTIOXA_47" = "#F39C12",
  "restricted_vehicle" = "#3498DB",
  "restricted_RTIOXA_47" = "#5DADE2"
)

# --- Plot Locomotion ---
loc_plot <- loc_plotdata %>%
  ggplot(aes(x = SABLE, y = total_act, fill = PlotGroup)) +
  
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
  labs(y = "Locomotion (m/day)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~STRAIN)

loc_plot

# --- Mixed model for locomotion ---
loc_model <- lmer(total_act ~ STRAIN* SABLE * GROUP * DRUG  + (1|ID), data = sable_locomotion_data)
summary(loc_model)

# --- Estimated marginal means ---
loc_emm <- emmeans(loc_model, pairwise ~ STRAIN* SABLE * GROUP * DRUG, adjust = "tukey")

# --- Convert to dataframes for results ---
df_emm_loc <- as.data.frame(loc_emm$emmeans)
df_pairs_loc <- as.data.frame(loc_emm$contrasts)

# --- Keep significant results only ---
df_sig_loc <- df_pairs_loc %>%
  filter(p.value <= 0.05)

print(df_sig_loc, n = Inf)
View(df_sig_loc)

