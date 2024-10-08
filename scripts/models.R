pacman::p_load(
    tidyverse,
    ggplot2
)

# import data ----
META <- read_csv("data/eegap/META.csv") %>% select(-LOC) %>% unique()
EMRI_RAW <- read_csv("data/eegap/EMRI.csv")
BW_RAW <- read_csv("data/eegap/BW.csv")
FI_RAW <- read_csv("data/eegap/FI.csv")
DELTA_BW <- BW_RAW %>% 
    mutate(num_date = as.numeric(as.factor(DATE))) %>% 
    group_by(ID) %>% 
    group_split() %>% 
    map_dfr(
        ., function(X){
            mdl <- lm(BW ~ num_date, data = X) %>% 
                broom::tidy() %>% 
                mutate(ID = X$ID[1])
        }
    ) %>% 
    mutate(
        term = if_else(term=="num_date", "BW_ROC", "INITIAL_BW")
    ) %>% 
    select(term, estimate, ID) %>% 
    pivot_wider(
        names_from = "term",
        values_from = "estimate"
    )
DELTA_FI <- FI_RAW %>% 
    mutate(num_date = as.numeric(as.factor(DATE))) %>% 
    group_by(ID) %>% 
    group_split() %>% 
    map_dfr(
        ., function(X){
            mdl <- lm(KCAL ~ num_date, data = X) %>% 
                broom::tidy() %>% 
                mutate(ID = X$ID[1])
        }
    ) %>% 
    mutate(
        term = if_else(term=="num_date", "FI_ROC", "INITIAL_FI")
    ) %>% 
    select(term, estimate, ID) %>% 
    pivot_wider(
        names_from = "term",
        values_from = "estimate"
    )
EMRI <- EMRI_RAW %>% 
    pivot_wider(
        names_from = "PARAM",
        values_from = "VALUE"
    ) %>% 
    mutate(
        ADIPOSITY_INDEX = Fat / (Fat + Lean + FreeWater),
        FAT_DELTA = (Fat - Fat[DATE=="2024-07-02"])/Fat[DATE=="2024-07-02"],
        LEAN_DELTA = (Lean - Lean[DATE=="2024-07-02"])/Lean[DATE=="2024-07-02"]
    ) %>% 
    drop_na() %>% 
    pivot_longer(
        cols = 4:10,
        names_to = "PARAM",
        values_to = "VALUE"
    ) %>% 
    left_join(META, by = c("ID")) %>%
    left_join(DELTA_BW, by = c("ID")) %>% 
    left_join(DELTA_FI, by = c("ID")) %>% 
    mutate(
        ID = as.factor(ID),
        PARAM = as.factor(PARAM),
        DIET = as.factor(DIET),
        SEX = as.factor(SEX)
    ) %>% 
    filter(DATE == "2024-08-26")

# all PARAM ~ DIET

mdl_data <- EMRI %>% 
    group_by(PARAM) %>% 
    group_split()

 # mdl 0 ~ DIET -----
mdls_0 <- mdl_data %>% 
    map_dfr(
        ., function(PARAM){
            mdl <- lm(
                data = PARAM,
                VALUE ~ DIET
            )
            emm <- emmeans::emmeans(
                mdl,
                pairwise ~ DIET,
                type = "response"
            )$contrast %>% 
                as_tibble() %>% 
                mutate(
                    PARAM = PARAM$PARAM[1]
                )
            return(emm)
        }
    )
mdls_0

 # mdl 1 ~ DIET + Weight -----
mdls_1 <- mdl_data %>% 
    map_dfr(
        ., function(PARAM){
            mdl <- lm(
                data = PARAM,
                VALUE ~ DIET + Weight
            )
            emm <- emmeans::emmeans(
                mdl,
                pairwise ~ DIET + Weight,
                type = "response"
            )$contrast %>% 
                as_tibble() %>% 
                mutate(
                    PARAM = PARAM$PARAM[1]
                )
            return(emm)
        }
    )
mdls_1

 # mdl 2 ~ DIET + BW_ROC -----
mdls_2 <- mdl_data %>% 
    map_dfr(
        ., function(PARAM){
            mdl <- lm(
                data = PARAM,
                VALUE ~ DIET + BW_ROC
            )
            emm <- emmeans::emmeans(
                mdl,
                pairwise ~ DIET + BW_ROC,
                type = "response"
            )$contrast %>% 
                as_tibble() %>% 
                mutate(
                    PARAM = PARAM$PARAM[1]
                )
            return(emm)
        }
    )
mdls_2

# mdl 3 ~ DIET + FI_ROC ----
mdls_3 <- mdl_data %>% 
    map_dfr(
        ., function(PARAM){
            mdl <- lm(
                data = PARAM,
                VALUE ~ DIET + FI_ROC
            )
            emm <- emmeans::emmeans(
                mdl,
                pairwise ~ DIET + FI_ROC,
                type = "response"
            )$contrast %>% 
                as_tibble() %>% 
                mutate(
                    PARAM = PARAM$PARAM[1]
                )
            return(emm)
        }
    )
mdls_3


# BW_ROC HFD vs LFD ----
mdls_bw_roc <- lm(
    data = EMRI %>% filter(PARAM=="ADIPOSITY_INDEX"),
    VALUE ~ DIET * BW_ROC
)
summary(mdls_bw_roc)

mdls_bw_roc_emtrend <- emmeans::emtrends(
    mdls_bw_roc,
    pairwise ~ DIET,
    var = "BW_ROC"
)
mdls_bw_roc_emtrend

# plots ----
EMRI %>% 
    filter(PARAM == "ADIPOSITY_INDEX") %>% 
    ggplot(aes(
        DIET, VALUE
    )) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange"
    ) +
    geom_point()

EMRI %>% 
    filter(PARAM == "ADIPOSITY_INDEX") %>% 
    ggplot(aes(
        Weight, VALUE, color = DIET
    )) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)

EMRI %>% 
    filter(PARAM == "ADIPOSITY_INDEX") %>% 
    ggplot(aes(
        BW_ROC, VALUE, color = DIET
    )) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)

EMRI %>% 
    filter(PARAM == "ADIPOSITY_INDEX") %>% 
    ggplot(aes(
        FI_ROC, VALUE, color = DIET
    )) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)

BW_RAW %>% 
    left_join(META, by = c("ID")) %>% 
    ggplot(aes(
        DATE, BW, group = ID, color = DIET
    )) +
    geom_line()

FI_RAW %>% 
    drop_na() %>% 
    left_join(META, by = c("ID")) %>% 
    ggplot(aes(
        DATE, (KCAL), group = ID, color = DIET
    )) +
    geom_line()
