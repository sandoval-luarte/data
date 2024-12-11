pacman::p_load(
    tidyverse,
    ggplot2
)

library(ggpubr)

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



