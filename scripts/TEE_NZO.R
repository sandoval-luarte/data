#we aim to evaluate the basal SPA of NZO female prior to LFD

#Libraries
library(dplyr) #to open a RDS and use pipe

sable_hr_data <- readRDS(file = "../data/sable/sable_hr_data.rds") %>% # Load the data
  filter(COHORT > 2 & COHORT < 6) %>%  #we only want NZO mice
  mutate(lights  = if_else(hr %in% c(20,21,22,23,0,1,2,3,4,5), "off", "on")) %>% 
  filter(grepl("kcal_hr_*", parameter)) %>% # just to see TEE in kcal first
group_by(ID,sable_idx,lights) %>% 
  summarise(mean_hr_tee = mean(value))

# Create the plot
sable_hr_plot <- sable_hr_data %>%
  ggplot(aes(x = sable_idx, y = mean_hr_tee, color = as.factor(ID))) +  # or just value to see the hourly value
  geom_line(aes(group = ID)) +              # Connect lines across days for each ID
  geom_point(size = 3, alpha = 0.8) +                 # Add individual points
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) +
  facet_wrap(~lights) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    size = 0.5,
    shape = 21,
    color = "black",
    fill = "red",
    position = position_dodge(width = 0.2),  # Mean and SEM as big points,
    aes(group = sable_idx)
  ) +
  geom_text(aes(label = ID), size = 3, vjust = -1, hjust = 1) + 
  labs(
    title = "TEE",
    x = "Days in sable",
    y = "cumsum TEE (kcal/h)"
  ) +
  theme_minimal()                               # Separate by light condition
sable_hr_plot

mdl_ee <- lme4::lmer(
  data = sable_hr_data,
  (mean_hr_tee) ~ lights + (1|ID)
)
summary(mdl_ee)

emmeans::emmeans(
  mdl_ee,
  pairwise ~ lights,
  type = "response"
)

sable_hr_plot

