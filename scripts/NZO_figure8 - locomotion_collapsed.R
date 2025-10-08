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
  filter(COHORT %in% c(3, 4, 5)) %>%
  mutate(
    lights = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on"),
    SABLE = case_when(
      sable_idx %in% paste0("SABLE_DAY_", 1:7) ~ "baseline",
      sable_idx %in% paste0("SABLE_DAY_", 8:11) ~ "peak obesity",
      sable_idx %in% paste0("SABLE_DAY_", 12:15) ~ "BW loss",
      sable_idx %in% paste0("SABLE_DAY_", 16:19) ~ "BW maintenance",
      sable_idx %in% paste0("SABLE_DAY_", 20:23) ~ "BW regain"
    )
  ) %>%
  filter(grepl("AllMeters_*", parameter)) %>%
  group_by(ID, SABLE) %>%
  mutate(
    zt_time = zt_time(hr),
    is_zt_init = replace_na(as.numeric(hr != lag(hr)), 0),
    complete_days = cumsum(if_else(zt_time == 0 & is_zt_init == 1, 1, 0))
  ) %>%
  ungroup() %>%
  group_by(ID, complete_days) %>%
  mutate(is_complete_day = if_else(min(zt_time) == 0 & max(zt_time) == 23, 1, 0)) %>%
  ungroup() %>%
  
  # compute locomotor activity incrementally
  group_by(ID, complete_days, SABLE, lights) %>%
  mutate(loc_act = value - lag(value)) %>%
  filter(loc_act >= 0) %>%
  summarise(total_act = sum(loc_act), .groups = "drop") %>%
  
  # exclude problematic IDs (same as TEE)
  filter(!ID %in% c(3715,3712,3709,3717,3718,3723,3724,3725)) %>%
  filter(complete_days %in% c(1, 2)) %>%
  
  # average across 2 complete days
  group_by(ID, SABLE, lights) %>%
  summarise(total_act = mean(total_act), .groups = "drop") %>%
  
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
                   levels = c("baseline", "peak obesity", "BW loss", 
                              "BW maintenance", "BW regain"))
  )

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
  facet_wrap(~lights) 

loc_plot

# --- Mixed model for locomotion ---
loc_model <- lmer(total_act ~ SABLE * GROUP * DRUG *lights + (1|ID), data = sable_locomotion_data)
summary(loc_model)

# --- Estimated marginal means ---
loc_emm <- emmeans(loc_model, pairwise ~ SABLE * GROUP * DRUG*lights, adjust = "tukey")

# --- Convert to dataframes for results ---
df_emm_loc <- as.data.frame(loc_emm$emmeans)
df_pairs_loc <- as.data.frame(loc_emm$contrasts)

# --- Keep significant results only ---
df_sig_loc <- df_pairs_loc %>%
  filter(p.value <= 0.05)

print(df_sig_loc, n = Inf)
View(df_sig_loc)


# --- Plot Locomotion split by lights ---
loc_plot_lights <- loc_plotdata %>%
  ggplot(aes(x = SABLE, y = total_act, color = GROUP, group = ID)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", aes(group = GROUP), size = 1.2) +
  facet_wrap(~lights) +
  labs(y = "Locomotion (m/day)", color = "Group") +
  theme_minimal()

loc_plot_lights
