pacman::p_load(
    tidyverse,
    ggplot2,
    ggpubr,
    furrr
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

sable_data <- readRDS("../data/sable/sable_hr_data.rds")


## before after injection ----

before_after_analysis <- readRDS("../data/sable/before_after_analysis.rds")

## RQ analysis ----
mdl_RQ_data_ <- before_after_analysis %>% 
    select(ID, hr, datetime, parameter, corrected_value, event_flag, drug, COHORT) %>% 
    mutate(parameter = str_remove(parameter, "_[0-9]+")) %>% 
    filter(parameter == "RQ", COHORT == 6, ID != 1003) %>% 
    group_by(ID, parameter, drug, event_flag) %>% 
    mutate(
        time = as.numeric(as.factor(datetime)),
        RQ = corrected_value
    ) 

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

