# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(patchwork)

# ---- Read BW data ----
BW_data <- read_csv("../data/BW.csv") %>% 
  filter(!ID %in% c(3712, 3715)) %>%  # remove animals that died
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(
    # Define GROUP
    GROUP = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3717, 3716, 3719, 3718, 3726,
                7860, 7862, 7864, 7867, 7868, 7869, 7870, 7871, 7873, 7875, 
                7876, 7879, 7880, 7881, 7882, 7883) ~ "ad lib",
      ID %in% c(3708, 3714, 3720, 3721, 3710, 3722, 3723, 3724, 3725, 3727, 3728, 3729,
                7861, 7863, 7865, 7866, 7872, 7874, 7877, 7878) ~ "restricted"
    ),
    # Define DRUG
    DRUG = case_when(
      ID %in% c(3706, 3707, 3709, 3711, 3713, 3714, 3720, 3724, 3725, 3727, 3728,
                7861, 7863, 7864, 7878, 7867, 7872, 7875, 7876, 7869, 7870, 7871,
                7868, 7880, 7881, 7882, 7883) ~ "vehicle",
      ID %in% c(3708, 3710, 3716, 3717, 3718, 3719, 3721, 3722, 3723, 3726, 3729,
                7862, 7865, 7873, 7874, 7877, 7866, 7879, 7860) ~ "RTIOXA_47"
    ),
    # Define STATUS
    STATUS = case_when(
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_1_INJECTIONS" ~ "BW maintenance",
      STRAIN == "C57BL6/J" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      STRAIN == "NZO/HlLtJ" & COMMENTS == "DAY_4_SABLE_AND_SAC" ~ "BW regain",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(STATUS)) %>%
  ungroup()

# ---- 5-week RTIOXA-47 delta ----
BW_delta <- BW_data %>%
  filter(STATUS %in% c("BW maintenance", "BW regain")) %>%
  filter(DRUG %in% c("vehicle", "RTIOXA_47")) %>%
  select(ID, DRUG, STATUS, BW, STRAIN, SEX, GROUP) %>%
  pivot_wider(names_from = STATUS, values_from = BW) %>%
  filter(!is.na(`BW maintenance`) & !is.na(`BW regain`)) %>%
  group_by(ID, SEX, STRAIN, GROUP, DRUG) %>%
  mutate(
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47")),
    bw_rel = 100 * (`BW regain` - `BW maintenance`) / `BW maintenance`
  ) %>%
  ungroup()

# ---- Pairwise p-values for 5-week plot ----
pairwise_5W <- BW_delta %>%
  group_by(STRAIN, SEX, GROUP) %>%
  summarise(
    p_val = t.test(bw_rel[DRUG=="vehicle"], bw_rel[DRUG=="RTIOXA_47"])$p.value,
    .groups = "drop"
  )

# ---- 5-week plot (split by GROUP) ----
plot_BW_47_5W <- ggplot(BW_delta, aes(x = DRUG, y = bw_rel, fill = DRUG)) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_point(position = position_jitter(width = 0.15), size = 2, alpha = 0.7) +
  facet_grid(GROUP ~ STRAIN + SEX) +
  scale_fill_manual(values = c("vehicle" = "gray70", "RTIOXA_47" = "#8DA0CB")) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "% BW change from BW maintenance") +
  geom_text(
    data = pairwise_5W,
    aes(x = 1.5, y = max(BW_delta$bw_rel) + 3, label = paste0("p = ", signif(p_val,3))),
    inherit.aes = FALSE
  ) +
  theme(legend.position = "none")

# ---- 5-day RTIOXA-43 ----
BW_data_43 <- read_csv("../data/BW.csv") %>% 
  filter(COHORT == 9) %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    bw_rel = 100 * (BW - first(BW)) / first(BW),
    DRUG = case_when(
      ID %in% c(8096, 8099, 8102) ~ "RTI_43_Y",
      ID %in% c(8097, 8098, 8101) ~ "RTI_43_M",
      ID %in% c(8095, 8100, 8103) ~ "Vehicle"
    ),
    day_rel = DATE - first(DATE)
  ) %>%
  filter(as.Date(DATE) >= as.Date("2025-04-16"))

BW_compare_43 <- BW_data_43 %>%
  filter(COMMENTS %in% c("DAY_1_INJECTIONS", "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION")) %>%
  mutate(COMMENTS = recode(COMMENTS,
                           "DAY_1_INJECTIONS" = "baseline",
                           "SAC_AND_WHITE_DEPOSIT_QUANTIFICATION" = "endpoint")) %>%
  select(ID, DRUG, COMMENTS, bw_rel) %>%
  pivot_wider(names_from = COMMENTS, values_from = bw_rel) %>%
  filter(!is.na(endpoint)) %>%
  mutate(DRUG = factor(DRUG, levels = c("Vehicle", "RTI_43_Y", "RTI_43_M")))

# Compute p-values
pairwise_43 <- BW_compare_43 %>%
  group_by(DRUG) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(n >= 2) # keep groups >=2

t_vehicle_Y <- t.test(BW_compare_43$endpoint[BW_compare_43$DRUG=="Vehicle"],
                      BW_compare_43$endpoint[BW_compare_43$DRUG=="RTI_43_Y"])$p.value
t_vehicle_M <- t.test(BW_compare_43$endpoint[BW_compare_43$DRUG=="Vehicle"],
                      BW_compare_43$endpoint[BW_compare_43$DRUG=="RTI_43_M"])$p.value

plot_BW_43_5D <- ggplot(BW_compare_43, aes(x = DRUG, y = endpoint, fill = DRUG)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
  geom_text(aes(x = 1.5, y = max(endpoint) + 5, label = paste0("p = ", signif(t_vehicle_Y,3))),
            inherit.aes = FALSE) +
  geom_text(aes(x = 2.5, y = max(endpoint) + 10, label = paste0("p = ", signif(t_vehicle_M,3))),
            inherit.aes = FALSE) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "", y = "% BW gain from baseline") +
  scale_fill_manual(values = c("gray70", "#66C2A5", "darkgreen"))

# ---- 5-day RTIOXA-47 ----
BW_data_47 <- read_csv("../data/BW.csv") %>%
  filter(COHORT == 12) %>%
  group_by(ID) %>%
  arrange(DATE) %>%
  mutate(
    DRUG = case_when(
      ID %in% c(8075, 8077, 8078) ~ "RTIOXA_47",
      ID %in% c(8074, 8076, 8079) ~ "vehicle"
    ),
    STRAIN = "C57BL6/J"
  ) 

BW_compare_47 <- BW_data_47 %>%
  filter(COMMENTS %in% c("INJECTION_DAY_1", "ECHO_MRI_AND_SAC")) %>%
  pivot_wider(id_cols = c(ID, DRUG, STRAIN, SEX),
              names_from = COMMENTS, values_from = BW) %>%
  mutate(
    bw_rel = 100 * (ECHO_MRI_AND_SAC - INJECTION_DAY_1)/INJECTION_DAY_1,
    DRUG = factor(DRUG, levels = c("vehicle", "RTIOXA_47"))
  )

p_47 <- t.test(BW_compare_47$bw_rel[BW_compare_47$DRUG=="vehicle"],
               BW_compare_47$bw_rel[BW_compare_47$DRUG=="RTIOXA_47"])$p.value

plot_BW_47_5D <- ggplot(BW_compare_47, aes(x = DRUG, y = bw_rel, fill = DRUG)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
  geom_text(aes(x = 1.5, y = max(bw_rel)+3, label = paste0("p = ", signif(p_47,3))),
            inherit.aes = FALSE) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "", y = "% BW change from baseline") +
  scale_fill_manual(values = c("gray70", "#8DA0CB"))

# ---- Combine plots ----
combined_BW <- (plot_BW_43_5D | plot_BW_47_5D) / plot_BW_47_5W +
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = "A")

combined_BW
