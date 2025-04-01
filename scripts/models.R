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


#format plot
format.plot <- theme_pubr() +
    theme(strip.background = element_blank(), 
          #    strip.text = element_blank(),
          plot.margin = unit(rep(0.2,4), "cm"), 
          #legend.position = "none",
          axis.text.x = element_text(family = "Helvetica", size = 14),
          axis.text.y = element_text(family = "Helvetica", size = 14),
          axis.title.y.left =element_text(family = "Helvetica", size = 16),
          axis.title.x.bottom =element_text(family = "Helvetica", size = 16))

# import data ----
## remember to set the path to script location
META <- read_csv("../data/META.csv")
BW_RAW <- read_csv("../data/BW.csv")
FI_RAW <- read_csv("../data/FI.csv")
ECHOMRI_RAW <- read_csv("../data/echomri.csv")

# food intake analysis ----
FI <- FI_RAW %>% 
    group_by(ID) %>% 
    arrange(DATE, .by_group = TRUE) %>% 
    mutate(
        REL_DATE = as.numeric(as.factor(DATE))
    ) %>% 
    group_by(ID) %>% 
    arrange(.by_group = TRUE) %>% 
    drop_na(corrected_intake_kcal) %>% 
    mutate(
        cumulative_intake_kcal = cumsum(corrected_intake_kcal)
    )

# food intake plot ----

FI %>% 
    filter(COHORT > 1) %>% 
    ggplot(aes(
        REL_DATE, cumulative_intake_kcal, group = ID, color = DIET_FORMULA
    )) +
    geom_line() +
    facet_wrap(~SEX*STRAIN)+
    format.plot

mdl_cumintake_slope <- lmerTest::lmer(
    data = FI %>% filter(SEX == "F", DIET_FORMULA == "D12450Ki_research_diets"),
    cumulative_intake_kcal ~ REL_DATE * STRAIN + (1 | ID)
)
summary(mdl_cumintake_slope)

emm_cumintake_slope <- emmeans::emtrends(
    mdl_cumintake_slope,
    pairwise ~ STRAIN,
    var = "REL_DATE"
)
emm_cumintake_slope


# food intake and body weight gain correlation ----

BW_RAW %>% 
    group_by(ID) %>% 
    summarise(gain_bw = tail(BW, 1) - head(BW, 1),
              time_elapsed = as.numeric(tail(DATE, 1) - head(DATE, 1)),
              bw_gain_per_day = gain_bw / time_elapsed)

mdl_weight_gain <- BW_RAW %>% 
    group_by(ID) %>% 
    mutate(date_diff = cumsum(replace_na(as.numeric(DATE - lag(DATE)), 0))) %>% 
    group_split() %>% 
    map(
        ., function(X){
            mdl <- lm(data=X, BW ~ date_diff)
            return(broom::tidy(mdl) %>% slice(2) %>% mutate(ID = X$ID[1]))
        }
    ) %>% 
    bind_rows()
mdl_weight_gain

mdl_fi <- FI_RAW %>% 
    filter(DIET != "CHOW") %>% 
    drop_na(corrected_intake_kcal) %>% 
    group_by(ID) %>% 
    summarise(daily_intake_kcal = mean(corrected_intake_kcal))
mdl_fi

fi_bw_corr <- left_join(mdl_fi, mdl_weight_gain, by = "ID") %>% 
    rename(daily_weight_gain = estimate) %>% 
    left_join(., META, by = "ID")

fi_bw_corr %>% 
    ggplot(aes(daily_intake_kcal, daily_weight_gain)) +
    geom_point(aes(color = STRAIN)) +
    geom_smooth(aes(group = STRAIN, color = STRAIN), method = "lm") +
    facet_wrap(~SEX) + 
    format.plot

# echomri ----

adiposity_index_gain <- echomri_data %>% 
    filter(COHORT > 1) %>% #For now we want to evaluate what happen with these cohorts
    ungroup() %>% 
    group_by(ID, SEX, STRAIN, DIET_CODE) %>% 
    arrange(Date, .by_group = TRUE) %>% 
    summarise(
        adi_gain = tail(adiposity_index, 1) - head(adiposity_index, 1),
        time_elapsed = as.numeric(tail(Date, 1) - head(Date, 1)),
        corrected_gain = (adi_gain / time_elapsed) * 100
    ) %>% 
    filter(time_elapsed != 0)
adiposity_index_gain

# plot adiposity index ----

echomri_data <- ECHOMRI_RAW%>% 
    filter(COHORT > 1) 

echomri_data %>% 
    ggplot(aes(
        n_measurement, adiposity_index,
        color = interaction(STRAIN, SEX)
    )) +
    geom_point(alpha = 0.25) +
    geom_line(aes(group = ID), alpha = 0.25) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        aes(group = interaction(STRAIN, SEX))
    ) +
    facet_wrap(~DIET_CODE)+
    format.plot

adiposity_index_gain %>% 
    ggplot(aes(
        interaction(STRAIN, SEX), corrected_gain
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    facet_wrap(~DIET_CODE)+
    format.plot

adiposity_index_gain %>% 
    group_by(SEX, STRAIN, DIET_CODE) %>% 
    summarise(count = n())


adiposity_index_gain %>% 
    group_by(SEX, STRAIN, DIET_CODE) %>% 
    summarise(
        m = mean(corrected_gain),
        se = sd(corrected_gain)/sqrt(n())
    )

# bodyweight over time (derivative) ----

bw_derivative <- BW_RAW %>% 
    group_by(ID) %>% 
    mutate(
        DELTA_TIME = replace_na(as.numeric(DATE - lag(DATE)), 0),
        DELTA_BW = replace_na(as.numeric(BW - lag(BW)), 0),
        BW_TIME_DERIVATIVE = (DELTA_BW / DELTA_TIME),
        REL_DATE = cumsum(replace_na(as.numeric(DATE - lag(DATE)), 0)),
        REL_DATE_BIN = as.factor(cut(REL_DATE, breaks = 3, labels = FALSE)),
        tmp = log(BW_TIME_DERIVATIVE)
    )
bw_derivative

bw_deriv_nzo <- bw_derivative %>% 
    filter(COHORT %in% c(2, 3, 4, 5)) %>% 
    drop_na(BW_TIME_DERIVATIVE)
bw_deriv_nzo

last_n_days_nzo_deriv <- bw_deriv_nzo %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    arrange(DATE, .by_group = TRUE) %>% 
    mutate(time_from_last_measurement = abs(as.numeric(DATE - tail(DATE, 1)))) %>% 
    filter(time_from_last_measurement <= 30)
last_n_days_nzo_deriv

last_n_slopes <- last_n_days_nzo_deriv %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    group_split() %>% 
    map(., function(X){
        mdl <- lm(data = X, BW ~ REL_DATE)
        out <- broom::tidy(mdl) %>% 
            mutate(ID = X$ID[1])
    }) %>% 
    bind_rows() %>% 
    filter(term == "REL_DATE") %>% 
    select(-term)
last_n_slopes
    

baseline_derivative <- bw_derivative %>% filter(REL_DATE < 50) %>% pull(DELTA_BW)

baseline_mean <- mean(baseline_derivative)
baseline_sd <- sd(baseline_derivative)
baseline_mean+(2*baseline_sd)
baseline_mean-(2*baseline_sd)

bw_derivative %>% 
    filter(COHORT %in% c(3, 4, 5)) %>% 
    ggplot(aes(x = REL_DATE, y = BW)) +
    geom_line(aes(group = ID)) +
    geom_smooth(method = "lm", aes(group = ID), se = FALSE) +
    facet_wrap(~ID, scales = "free")

last_n_days_nzo_deriv %>% 
    filter(COHORT %in% c(2,3, 4, 5)) %>% 
    left_join(., last_n_slopes, by = "ID") %>%
    mutate(gain_last_n = if_else(estimate > 0, "Gained", "Lost")) %>% 
    ggplot(aes(x = REL_DATE, y = BW)) +
    geom_line(aes(group = ID)) +
    geom_point(aes(group = ID, color = gain_last_n), size = 2) +
    geom_smooth(method = "lm", aes(group = ID, color = gain_last_n), se = FALSE) +
    facet_wrap(~ID, scales = "free") +
    scale_color_manual(values = c("Gained" = "darkgreen", "Lost" = "red"))


bw_derivative %>% 
    filter(COHORT %in% c(3, 4, 5)) %>% 
    ggplot(aes(REL_DATE, BW_TIME_DERIVATIVE)) +
    geom_point(aes(group = ID)) +
    geom_line(aes(group = ID)) +
    geom_hline(yintercept = 0) +
    facet_wrap(~ID, scales = "free_x")

bw_deriv_nzo %>% 
    left_join(., last_n_slopes, by = "ID") %>%
    mutate(gain_last_n = if_else(estimate > 0, "Gained", "Lost")) %>% 
    ggplot(aes(REL_DATE, BW_TIME_DERIVATIVE)) +
    geom_point(aes(group = ID, color = gain_last_n)) +
    geom_line(aes(group = ID)) +
    geom_hline(yintercept = 0) +
    facet_wrap(~ID, scales = "free_x") +
    scale_color_manual(values = c("green", "red"))


bw_derivative %>% 
    filter(REL_DATE < 50) %>% 
    ggplot(aes(BW_TIME_DERIVATIVE)) +
    geom_histogram()

# sable analysis ----

setwd(this.path::here())
sable_data <- readRDS("../data/sable/sable_downsampled_data.rds")
sable_data_injections <- readRDS("../data/sable/before_after_analysis.rds")

## test ----

test <- sable_data_injections %>% 
    filter(
        parameter %in% c("RQ", "AllMeters"),
        event_flag == "after"
    ) %>% 
    group_by(ID, drug, parameter) %>% 
    mutate(
        vel = (corrected_value-lag(corrected_value)),
        acc = (vel-lag(vel)),
        jerk = (acc-lag(acc)),
        scaled_inst_change = scale(vel)
    ) %>% 
    ungroup() %>% 
    group_by(ID, parameter) %>% 
    mutate(
        quantile = percent_rank(jerk)
    ) %>% 
    ungroup()

t <- test %>% filter(ID==2006,drug=="veh") %>% 
    mutate(
        fft_ = abs(fft(corrected_value)))

t %>% 
    filter(time_from_injection<100) %>% 
    ggplot(aes(
        time_from_injection, fft_
    )) +
    geom_point()


test %>% 
    filter(parameter=="AllMeters", time_from_injection>240) %>% 
    group_by(ID, SEX, drug) %>% 
    mutate(jerk_sum=cumsum(abs(replace_na(acc,0)))) %>% 
    ggplot(aes(
        time_from_injection, jerk_sum, group=interaction(ID,drug),color=drug
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~SEX*ID,scale="free")

grp_speed <- test %>% 
    group_by(drug) %>% 
    group_split() %>% 
    map_dfr(., function(X){
        X %>% 
            select(ID, SEX, drug, time_from_injection, parameter, inst_change) %>% 
            pivot_wider(names_from = "parameter",
                        values_from = "inst_change"
                        )
    })

grp_speed %>% 
    filter(RQ>0,RQ<1,
           drug!="RTI_43_Y") %>% 
    ggplot(aes(
        log(AllMeters+1), log(RQ)
    )) +
    geom_point() +
    facet_wrap(~drug*SEX) +
    geom_smooth(method="lm")

test %>% 
    filter(parameter=="kcal_hr") %>% 
    mutate(
        cnt = as.numeric(if_else(quantile>0.9,1,0)) %>% 
            replace_na(., 0)
    ) %>%
    ungroup() %>% 
    group_by(ID, drug) %>% 
    mutate(
        cum_cnt = cumsum(cnt)
    ) %>% 
    filter(time_from_injection<=240) %>% 
    ggplot(aes(
        time_from_injection, cum_cnt, color=drug
    )) +
    geom_line() +
    facet_wrap(~ID*SEX,scale="free")


## bodyweight ----

bw_sable <- sable_data_injections %>% 
    filter(
        parameter == "BodyMass",
        event_flag == "after"
    ) %>% 
    filter(corrected_value > 0) %>% 
    ungroup() %>% 
    mutate(
        hr_factor = trunc(time_from_injection/60) %>% as.factor(),
        time = as.numeric(time_from_injection)
    )

bw_mdl <- lme4::lmer(
    data = bw_sable,
    corrected_value ~ drug * SEX * time + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
summary(bw_mdl)

bw_emm <- emmeans::emtrends(
    bw_mdl,
    pairwise ~ drug | SEX,
    var = "time",
    type = "response"
)
bw_emm

bw_sable %>% 
    group_by(ID, drug) %>% 
    mutate(
        delta_bw = corrected_value - corrected_value[1]
    ) %>% 
    ggplot(aes(
        time, delta_bw,
        group=interaction(ID, drug), color=drug
    )) +
    geom_line(alpha=0.25) +
    geom_smooth(method="lm", se=FALSE, aes(group=drug)) +
    facet_wrap(~SEX)


## total energy expenditure ----

tee_correction <- sable_data_injections %>% 
    filter(
        parameter %in% c("kcal_hr", "BodyMass"),
        event_flag == "after",
        drug == "veh"
    ) %>% 
    ungroup() %>% 
    group_by(ID, drug) %>% 
    select(ID, drug, DateTime, parameter, corrected_value, time_from_injection) %>% 
    pivot_wider(names_from = "parameter", values_from = "corrected_value") %>% 
    filter(kcal_hr>0.2, BodyMass>0, ID!=2001)
tee_correction

tee_corr_mdl <- lme4::lmer(
    data = tee_correction,
    kcal_hr ~ BodyMass + time_from_injection + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
summary(tee_corr_mdl)

bw <- sable_data_injections %>% 
    filter(parameter=="BodyMass",
           event_flag=="after") %>% 
    pull(corrected_value)
intake <- sable_data_injections %>% 
    filter(parameter=="FoodA",
           event_flag=="after") %>% 
    pull(corrected_value)
tee <- sable_data_injections %>% 
    filter(
        parameter == "kcal_hr",
        event_flag == "after"
    ) %>% 
    mutate(BodyMass = bw,
           intake_gr = intake,
           kcal_hr = corrected_value) %>% 
    ungroup() %>% 
    group_by(ID, drug) %>% 
    filter(corrected_value>0) %>% 
    mutate(
        cum_tee = cumsum(corrected_value/60),
        cum_intake = abs(intake_gr-lag(intake_gr)) %>% replace_na(.,0) %>% cumsum(),
        energy_balance = (cum_intake)-(cum_tee),
        hr_factor = trunc(time_from_injection/60) %>% as.factor(),
        time = as.numeric(time_from_injection)
    ) %>% 
    ungroup() %>% 
    mutate(
        tee_mdl = predict(tee_corr_mdl, newdata=., re.form=NA),
        tee_resid = kcal_hr - tee_mdl,
        tee_resid_cumsum = cumsum(tee_resid)
    )
tee

tee %>% 
    ggplot(aes(
        time, energy_balance,group=interaction(ID,drug),color=drug
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~SEX)

total_tee <- tee %>% 
    ungroup() %>% 
    group_by(ID, drug) %>% 
    summarise(
        kcal = max(cum_tee)
    )


total_intake <- sable_data_injections %>% 
    filter(
        parameter == "FoodA",
        event_flag == "after"
    ) %>% 
    ungroup() %>% 
    group_by(ID, drug, SEX) %>% 
    mutate(
        intake = max(corrected_value) - corrected_value,
        time = as.numeric(time_from_injection),
        hr_factor = as.factor(time/60)
    )

intake_mdl <- lme4::lmer(
    data = total_intake,
    intake ~ drug * SEX * time + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)

total_intake %>% 
    ungroup() %>% 
    mutate(preds = predict(intake_mdl, re.form=NA)) %>% 
    ggplot(aes(
        time, intake,
        group=interaction(ID, drug), color=drug
    )) +
    geom_line(alpha=0.25) +
    geom_line(aes(time, preds), linewidth=2) +
    facet_wrap(~SEX)

intake_emm <- emmeans::emmeans(
    intake_mdl,
    pairwise ~ drug | SEX,
    at = list(time=c(1440)),
    type="response"
)
intake_emm

intake_emm <- emmeans::emtrends(
    intake_mdl,
    pairwise ~ drug | SEX,
    var = "time",
    type="response"
)
intake_emm


total_intake %>%
    ggplot(aes(
        time, intake,
        group=interaction(ID,drug), color=drug
    )) +
    geom_line() +
    facet_wrap(~SEX)

intake_kcal <- total_tee %>% 
    left_join(., total_intake, by = c("ID", "drug"))

intake_kcal %>% 
    ggplot(aes(
        drug, intake-(kcal/60), color=drug
    )) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE,aes(group=drug)) +
    facet_wrap(~SEX)

tee %>% 
    ggplot(aes(
        corrected_value, group=drug, color=drug
    )) +
    geom_density() +
    facet_wrap(~SEX)

tee %>% 
    filter(time < 100) %>% 
    ggplot(aes(
        time, corrected_value,
        group=interaction(ID,drug),color=drug
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~ID*SEX, scale="free_y")

tee_mdl <- lmer(
    data = tee %>% filter(!(ID %in% c(2001, 2006))),
    energy_balance ~ drug * SEX * hr_factor + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)


tee_emm <- emmeans::emmeans(
    tee_mdl,
    pairwise ~ drug | SEX * hr_factor,
    type = "response",
    lmer.df = "satterthwaite",
    lmerTest.limit=13044
)
tee_emm

tee %>% 
    filter(hr_factor == "4") %>%
    mutate(
        drug = factor(as.factor(drug), levels=c("veh","RTI_43_M","RTI_43_Y"))
    ) %>% 
    ungroup() %>% 
    group_by(ID, drug, SEX) %>% 
    summarise(
        energy_balance = max(energy_balance)
    ) %>% 
    ggplot(aes(
        drug, energy_balance
    )) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_line(aes(group=ID)) +
    facet_wrap(~SEX) +
    ggtitle("4 hrs energy balance")

tee_emm_p <- tee_emm$emmeans %>% 
    broom::tidy(conf.int=TRUE)
tee_emm_p

p1 <- tee %>% 
    ggplot(aes(
        time, energy_balance, color = drug
    )) +
    geom_line(alpha=0.25, linewidth=3, aes(group = interaction(ID, drug))) +
    facet_wrap(~SEX*ID,scale="free") 
p1

## microstructure ----

microstructure <- sable_data_injections %>% 
    filter(parameter=="FoodA",
           event_flag=="after",
           ID != 2004) %>% 
    group_by(ID, drug) %>% 
    mutate(
        is_bout = if_else(corrected_value<lag(corrected_value),"bout", "stable") %>% 
            replace_na(., "stable"),
        bout_tag = data.table::rleid(is_bout),
        hr_factor = trunc(time_from_injection/60) %>% as.factor(),
    ) %>% 
    filter(is_bout == "bout") %>% 
    ungroup() %>% 
    group_by(ID, drug, bout_tag, SEX) %>% 
    summarise(
        init = min(time_from_injection),
        end = max(time_from_injection),
        centroid = (init+end)/2,
        bout_length_sec = as.numeric(max(DateTime) - min(DateTime)),
        intake_gr = max(corrected_value)-min(corrected_value)
    ) %>%
    ungroup() %>% 
    group_by(ID, drug) %>% 
    mutate(
        cum_bouts = row_number()
    )
microstructure


microstructure_cnt_mdl <- lme4::lmer(
    data = microstructure,
    cum_bouts ~ drug * SEX * centroid + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
microstructure_cnt_mdl

microstructure_cnt_emm <- emmeans::emmeans(
    microstructure_cnt_mdl,
    pairwise ~ drug | SEX*centroid,
    at = list(centroid = seq(0, 1500, 60)),
    type = "response"
)
microstructure_cnt_emm

microstructure_cnt_p <- microstructure_cnt_emm$emmeans %>% 
    broom::tidy(conf.int = TRUE)
microstructure_cnt_p


p1 <- microstructure %>% 
    ggplot(aes(
        centroid, cum_bouts,
        group = interaction(ID, drug), color = drug
    )) +
    geom_line(alpha=0.25) +
    geom_pointrange(
        data = microstructure_cnt_p,
        aes(centroid, estimate,ymin=conf.low,ymax=conf.high,
            group=drug,fill=drug)
    ) +
    geom_pointrange(
        data = microstructure_cnt_p,
        aes(centroid, estimate,ymin=conf.low,ymax=conf.high,
            group=drug,fill=drug)
    ) +
    facet_wrap(~SEX) +
    ggpubr::theme_pubr() +
  labs(x = "Time after injection (min)", y = "Cumulative number of bouts")

p1

microstructure_len <- lme4::lmer(
    data = microstructure,
    bout_length_sec ~ drug * SEX * centroid + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)

microstructure_len_emm <- emmeans::emmeans(
    microstructure_len,
    pairwise ~ drug | centroid * SEX,
    at = list(centroid = seq(0, 1500, 60)),
    type = "response"
)
microstructure_len_emm

microstructure_in <- lme4::lmer(
    data = microstructure %>%
      ungroup() %>% 
      group_by(ID, drug) %>% 
      mutate(cum_intake_gr=cumsum(intake_gr),
             hr_factor=trunc(centroid/60) %>% as.factor()),
    cum_intake_gr ~ drug * SEX * centroid + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)

microstructure %>%
  filter(drug != "RTI_43_Y") %>% 
  ungroup() %>% 
  group_by(ID, drug) %>% 
  mutate(cum_intake_gr=cumsum(intake_gr),
         hr_factor=trunc(centroid/60) %>% as.factor()) %>% 
  slice_max(order_by = cum_intake_gr, n=1) %>% 
  ggplot(aes(
    drug, cum_intake_gr
  )) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_line(aes(group=ID))+
  facet_wrap(~SEX)+
  labs(x = "Drug", y = "Cumulative intake in grams after 24h of injection")


microstructure_in_emm <- emmeans::emmeans(
  microstructure_in,
  pairwise ~ drug | SEX * hr_factor,
  type = "response"
)
microstructure_in_emm

microstructure_in_emm <- emmeans::emmeans(
    microstructure_in,
    pairwise ~ drug | centroid * SEX,
    at = list(centroid = seq(0, 1500, 60)),
    type = "response"
)
microstructure_in_emm

microstructure_in_emm$emmeans %>% 
  broom::tidy(conf.int=TRUE) %>% 
  ggplot(aes(
    centroid, estimate,
    group = drug, color = drug
  ))+
  geom_line() +
  facet_wrap(~SEX)


p1 <- microstructure %>% 
    ungroup() %>% 
    group_by(ID, hr_factor, drug) %>% 
    summarise(
        bout_length_sec = mean(bout_length_sec)
    ) %>% 
    ggplot(aes(
        drug, bout_length_sec,
        group = interaction(ID, drug), color = drug
    )) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        aes(group=drug)
    ) +
    geom_point(alpha=0.25) +
    facet_wrap(~hr_factor, scales="free_y")
p1

p2 <- microstructure %>% 
    ungroup() %>% 
    group_by(ID, hr_factor, drug) %>% 
    summarise(
        intake_gr = mean(intake_gr)
    ) %>% 
    ggplot(aes(
        drug, intake_gr,
        group = interaction(ID, drug), color = drug
    )) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        aes(group=drug)
    ) +
    geom_point(alpha=0.25) +
    facet_wrap(~hr_factor, scales="free_y")
p2

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
    )
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
allmeters_emm

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
view(allmeters_emm_contrasts)

p1 <- allmeters %>% 
  filter(hr_factor==360) %>% 
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
    ylab("Meters after injection") +
    xlab("Minutes from injection")
p1

p2 <- allmeters_emm_contrasts %>% 
    filter(hr_factor==1440) %>% 
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

## before after injection ----

before_after_analysis <- readRDS("../data/sable/before_after_analysis.rds")

rti_oxa43 <- before_after_analysis %>% 
  select(ID, hr, datetime, parameter, corrected_value, event_flag, drug) %>% 
  mutate(parameter2 = str_remove(parameter, "_[0-9]+")) %>% 
  filter(parameter == "RQ", ID != 1003) %>% 
  group_by(ID, parameter, drug, event_flag) %>% 
  mutate(
    time = as.numeric(as.factor(datetime)),
    RQ = corrected_value,
    time_cont = if_else(event_flag == "after", time + 12, time)
  )
  

### RQ analysis ----


#### RQ plots ----

rti_oxa43 %>% 
  filter(ID == 1004, drug == "RTI_43_Y", parameter == "RQ_6") %>% 
  ggplot(aes(event_flag, corrected_value)) +
  geom_point()

# compute the slope in the before and after injection use
# this to feed the model



## plot ----
mdl_RQ_data %>% 
    filter(ID != 1003) %>% 
    ggplot(aes(
        time, RQ - RQ_baseline, color = event_flag, group = interaction(as.factor(ID), event_flag)
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~drug)

## mdl RQ ----
mdl_RQ <- lmer(
    data = mdl_RQ_data,
    RQ ~ drug * time + RQ_baseline + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
summary(mdl_RQ)

RQ_emm <- emmeans::emmeans(
    mdl_RQ,
    pairwise ~ drug + RQ_baseline | time,
    type = "response"
)
RQ_emm

RQ_trend <- emmeans::emtrends(
    mdl_RQ,
    ~drug : RQ_baseline,
    var = "time",
    infer = TRUE
)
RQ_trend 
    

mdl_data_before_after <- before_after_analysis %>% 
    select(ID, hr, datetime, parameter, corrected_value, event_flag, drug) %>% 
    mutate(parameter = str_remove(parameter, "_[0-9]+")) %>% 
    filter(parameter == "AllMeters") %>% 
    group_by(ID, parameter, drug, event_flag) %>% 
    mutate(
        time = as.numeric(as.factor(datetime)),
        meters = corrected_value  - head(corrected_value, 1)
    )
mdl_data_before_after

mdl_data_before_after %>% 
    group_by(ID, drug) %>% 
    filter(event_flag == "before") %>% 
    mutate(delta_move = meters) %>% 
    filter(time == 12) %>% 
    ggplot(aes(
        ID, delta_move, color = as.factor(ID)
    )) +
    stat_summary(fun.data = "mean_se", geom = "pointrange") +
    geom_point()

mdl_data_before_after %>% 
    filter(event_flag == "after") %>% 
    ggplot(aes(
        time, meters
    )) +
    geom_point() +
    geom_line(aes(group = ID)) +
    geom_boxplot(
        data = mdl_data_before_after %>% 
            group_by(ID) %>% 
            filter(time == 12, event_flag == "after"),
        outlier.shape = NA,
        aes(color = drug)
    ) +
    facet_wrap(~drug, scales = "free_x") +
    scale_x_continuous(breaks = seq(0, 12, 1)) +
    ggpubr::theme_classic2() +
    ylab("Meters after injection") +
    theme(legend.position = "none")

analysis_per_hour <- mdl_data_before_after %>% 
    filter(time > 1, event_flag == "after") %>% 
    group_by(time) %>% 
    group_split() %>% 
    map(., function(X){
        mdl <- lmer(
            data = X,
            meters ~ drug + (1|ID),
            control = lmerControl(optimizer = "bobyqa"))
        emm <- emmeans::emmeans(
            mdl,
            pairwise ~ drug,
            type = "response"
        )
        pairwise_analysis <- broom::tidy(emm$contrasts, conf.int = TRUE) %>% 
            mutate(
                time = X$time[1]
            )
        marginals <- broom::tidy(emm$emmeans, conf.int = TRUE) %>% 
            mutate(
                time = X$time[1]
            )
        return(pairwise_analysis)
    })


mdl1 <- lmer(
    data = mdl_data_before_after %>% 
        filter(event_flag == "after"),
    meters ~ drug + (1|ID),
    control = lmerControl(optimizer = "bobyqa")
)
summary(mdl1)

mdl1_emm <- emmeans::emmeans(
    mdl1, 
    pairwise ~ drug,
    type = "response"
)
pairs(mdl1_emm, adjust = "none")

mdl1_emm$emmeans %>% 
    broom::tidy(conf.int = TRUE) %>% 
    ggplot(aes(
        drug, estimate, ymin = conf.low, ymax = conf.high
    )) +
    geom_pointrange() +
    ggpubr::theme_classic2()

bind_rows(analysis_per_hour) %>% 
    ggplot(aes(
        time, estimate, color = contrast,
        ymin = estimate - std.error, ymax = estimate + std.error
    )) +
    geom_hline(yintercept = 0) +
    geom_pointrange() +
    geom_line() +
    facet_wrap(~contrast) +
    scale_x_continuous(breaks = seq(2, 12, 1))

bind_rows(analysis_per_hour) %>% 
    ggplot(aes(
        time, estimate, color = drug,
        ymin = estimate - std.error, ymax = estimate + std.error
    )) +
    geom_pointrange() +
    geom_line() +
    scale_x_continuous(breaks = seq(2, 12, 1))

