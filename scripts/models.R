pacman::p_load(
    tidyverse,
    ggplot2
)

# import data ----
META <- read_csv("../data/META.csv")
BW_RAW <- read_csv("../data/BW.csv")
FI_RAW <- read_csv("../data/FI.csv")

# food intake analysis ----
FI <- FI_RAW %>% 
    drop_na() %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    mutate(
        REL_DATE = as.numeric(as.factor(DATE))
    ) %>% 
    group_by(ID) %>% 
    arrange(.by_group = TRUE) %>% 
    mutate(
        cumulative_intake_kcal = cumsum(corrected_intake_kcal)
    )

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
    facet_wrap(~SEX)