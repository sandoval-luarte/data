# libs ----
pacman::p_load(
    tidyverse,
    ggplot2,
    ggpubr,
    furrr,
    lme4,
    robustlmm,
    coxme,
    ggsurvfit,
    tidycmprsk
)

# import data ----
# sable analysis ----

setwd(this.path::here())
sable_data <- readRDS("./RTIOXA_43/sable_downsampled_data.rds")
sable_data_injections <- readRDS("./RTIOXA_43/before_after_analysis.rds")

## all meters ----

# 2004 has weird jump in all meters

allmeters <- sable_data_injections %>% 
    filter(parameter=="AllMeters",
           event_flag=="after",
           ID != 2004) %>% 
    group_by(ID, drug) %>% 
    mutate(
        injection_value = corrected_value[abs(time_from_injection) == min(abs(time_from_injection))][1],
        corrected_value_rel = (corrected_value - injection_value),
        corrected_value_perc = (corrected_value - injection_value)/injection_value,
        hr_factor = trunc(time_from_injection/60) %>% as.factor(),
        time_scaled = scale(time_from_injection)
    ) %>% 
  mutate(time_from_injection = as.numeric(gsub(" mins", "", time_from_injection))) 
  
  
allmeters

allmeters_mdl <- lme4::lmer(
    data = allmeters,
    corrected_value_rel ~ drug * SEX * hr_factor + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
summary(allmeters_mdl)

allmeters_mdl_2 <- lmerTest::lmer(
  data = allmeters,
  log(corrected_value_rel+1) ~ drug * SEX * time_from_injection + (1|ID),
  control = lmerControl(optimizer = "bobyqa")
)
summary(allmeters_mdl_2)


#use this values to describe the statistics in report 1 v1 - discrete
allmeters_emm <- emmeans::emmeans( 
    allmeters_mdl,
    pairwise ~ drug|SEX*hr_factor,
    type = "response",
    lmer.df = "satterthwaite"
)
allmeters_emm

allmeters_emtrend_2 <- emmeans::emtrends( 
    allmeters_mdl_2,
     pairwise ~ drug|SEX,
    var= "time_from_injection"
)
allmeters_emtrend_2

#use this values to describe the statistics in report 1 v1
allmeters_emm_2 <- emmeans::emmeans( 
  allmeters_mdl_2,
  pairwise ~ drug|SEX*time_from_injection,
  at= list(time_from_injection=c(1440)),
 regrid = "response"
)
allmeters_emm_2

allmeters_emm_p <- allmeters_emm$emmeans %>% 
    broom::tidy(conf.int=TRUE) %>% 
    mutate(
        hr_factor = as.numeric(hr_factor)*60
    )
allmeters_emm_p

allmeters_emm_contrasts <- allmeters_emm$contrasts %>% 
    broom::tidy(conf.int=TRUE) %>% 
    mutate(
        hr_factor = as.numeric(hr_factor)*60
    )
allmeters_emm_contrasts

pl <- allmeters %>%
  ggplot(aes(
    Time.from.injection, corrected_value_rel,
    color = drug
  )) +
  geom_line(aes(group = interaction(ID, drug)), alpha = 0.25) +
  geom_line(data = allmeters_emm, aes(Time.from.injection, estimate)) +
  geom_pointrange(data = allmeters_emm, aes(Time.from.injection, estimate,
                                            ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~SEX) +
  ggpubr::theme_pubr() +
  ylab(" SPA (meters)") +
  xlab("Minutes from injection") +
  scale_color_manual(values = c("veh" = "indianred1","RTI_43_M" = "springgreen3", "RTI_43_Y" = "cornflowerblue")) 
pl


p3 <- allmeters %>%
  ungroup() %>%
  group_by(time_from_injection, drug, SEX) %>%
  filter(time_from_injection %in% c(60, 240, 1440)) %>%
  ungroup() %>%
  select(parameter, value, ID, SEX, corrected_value, time_from_injection, drug, hr_factor, corrected_value_rel)

# Reorder the drug factor levels
p3$drug <- factor(p3$drug, levels = c("veh", "RTI_43_M", "RTI_43_Y"))

# Create the boxplot without connected IDs
ggplot(p3, aes(x = drug, y = corrected_value_rel, fill = drug)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(aes(color = drug), position = position_dodge(width = 0.75), size = 1.5) +
  labs(
    x = "Time from injection (minutes)",
    y = "SPA (meters)"
  ) +
  theme_pubr() +
  facet_wrap(~SEX*as.factor(time_from_injection),scales="free")+
  theme(
    strip.background = element_rect(fill = "white", color = NA),  # White background, no border
    strip.text = element_text(color = "black")  # Ensure text remains visible
  )

p4 <- allmeters %>%
  ungroup() %>%
  group_by(time_from_injection, drug, SEX) %>%
  filter(time_from_injection ==1380) %>%
  ungroup() %>%
  select(parameter, value, ID, SEX, corrected_value, time_from_injection, drug, hr_factor, corrected_value_rel)

# Reorder the drug factor levels
p4$drug <- factor(p4$drug, levels = c("veh", "RTI_43_M", "RTI_43_Y"))

# Create the boxplot without connected IDs
ggplot(p4, aes(x = drug, y = corrected_value_rel, fill = drug)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(aes(color = drug), position = position_dodge(width = 0.75), size = 1.5) +
  labs(
    x = "Time from injection (minutes)",
    y = "SPA (meters)"
  ) +
  theme_pubr() +
  facet_wrap(~SEX*as.factor(time_from_injection),scales="free")+
  theme(
    strip.background = element_rect(fill = "white", color = NA),  # White background, no border
    strip.text = element_text(color = "black")  # Ensure text remains visible
  )

