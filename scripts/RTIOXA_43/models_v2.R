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
sable_data <- readRDS("../data/sable/sable_downsampled_data.rds")
sable_data_injections <- readRDS("../data/sable/before_after_analysis.rds")

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

allmeters_emm <- emmeans::emmeans(
    allmeters_mdl,
    pairwise ~ drug|SEX*hr_factor,
    type = "response",
    lmer.df = "satterthwaite"
)
view(allmeters_emm)

allmeters_emm_p <- allmeters_emm$emmeans %>% 
    broom::tidy(conf.int=TRUE) %>% 
    mutate(
        hr_factor = as.numeric(hr_factor)*60
    )
view(allmeters_emm_p)

allmeters_emm_contrasts <- allmeters_emm$contrasts %>% 
    broom::tidy(conf.int=TRUE) %>% 
    mutate(
        hr_factor = as.numeric(hr_factor)*60
    )
view(allmeters_emm_contrasts)


p1<- allmeters%>%
ggplot(aes(
        time_from_injection, corrected_value_rel,
        color = drug
    )) +
    geom_line(aes(group = interaction(ID, drug)), alpha=0.25) +
    geom_line(data = allmeters_emm_p,
                aes(hr_factor, estimate)) +
    geom_pointrange(data = allmeters_emm_p,
                aes(hr_factor, estimate,
                    ymin=conf.low, ymax=conf.high)) +
    facet_wrap(~SEX) +
    ggpubr::theme_pubr() +
    ylab("SPA (Meters after injection)")+
  xlab("Time from injection (minutes)") 
  
p1

p2 <- allmeters_emm_contrasts %>% 
  filter(hr_factor==60) %>% 
    ggplot(aes(
        contrast, estimate,
        ymin=conf.low, ymax=conf.high
    )) +
    geom_hline(yintercept = 0) +
    geom_pointrange() +
    facet_wrap(~hr_factor*SEX) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") +
    scale_x_discrete(labels=c("RTI (M) - RTI (Y)",
                              "RTI (M) - Veh",
                              "RTI (Y) - Veh"))
p2


p3 <- allmeters %>%
  ungroup() %>%
  group_by(time_from_injection, drug, SEX) %>%
  filter(time_from_injection %in% c(60, 240, 1440)) %>%
  ungroup() %>%
  select(parameter, value, ID, SEX, corrected_value, time_from_injection, drug, hr_factor, corrected_value_rel)

# Reorder the drug factor levels
p3$drug <- factor(p3$drug, levels = c("veh", "RTI_43_M", "RTI_43_Y"))

# Create the boxplot without connected IDs
ggplot(p3, aes(x = as.factor(time_from_injection), y = corrected_value_rel, fill = drug)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(aes(color = drug), position = position_dodge(width = 0.75), size = 1.5) +
  labs(
    x = "Time from injection (min)",
    y = "SPA (meters)"
  ) +
  theme_pubr() +
  facet_wrap(~SEX)


## time course analysis ----

rel_sable_injections <- sable_data_injections %>% 
  group_by(ID, drug) %>% 
  mutate(
    parameter = `str_remove(parameter, "_[0-9]+")`,
    t = as.numeric(time_from_injection),
    injection_hour = hr[which.min(abs(t))],
    fix_after_injection = hr - injection_hour,
    fix_before_after = if_else(fix_after_injection>=0, "after", "before")
  ) %>% 
  ungroup() %>% 
  group_by(ID, drug, parameter) %>% 
  mutate(rel_corrected_value = corrected_value - corrected_value[which.min(abs(t))]) %>% 
  ungroup()


time_course_p1 <- rel_sable_injections %>% 
  filter(SEX == "F", parameter == "AllMeters", fix_before_after == "after") %>%
  ggplot(aes(
    fix_after_injection, rel_corrected_value, color = drug
  )) +
  geom_point(aes(group=interaction(ID, drug))) +
  geom_line(aes(group=interaction(ID, drug))) +
  facet_wrap(~ID, scales = "free")
time_course_p1

mdl_allmeters <- lmerTest::lmer(
  data = rel_sable_injections %>%
    filter(parameter=="AllMeters", SEX=="F", fix_before_after == "after"),
  rel_corrected_value ~ drug * fix_after_injection + (1|ID)
)
summary(mdl_allmeters)

emmeans::emtrends(
  mdl_allmeters,
  pairwise ~ drug,
  var = "fix_after_injection"
)

rel_sable_injections %>% 
  filter(fix_before_after == "after") %>% 
  group_by(ID, drug, parameter) %>% 
  slice_max(order_by = corrected_value, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(ID, drug, parameter, corrected_value) %>% 
  pivot_wider(., names_from = c("parameter"), values_from = "corrected_value") %>%
  ggplot(aes(
    log(AllMeters), log(kcal_hr)
  )) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~drug)

time_course_p2 <- rel_sable_injections %>% 
  ggplot(aes(
    time_from_injection, corrected_value, color = drug
  )) +
  geom_smooth(
    aes(group = drug),
    method = "loess",
    se = FALSE
  ) +
  facet_wrap(~`str_remove(parameter, "_[0-9]+")`*SEX, scales = "free")
time_course_p2
