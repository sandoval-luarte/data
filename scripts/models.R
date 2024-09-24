pacman::p_load(
    tidyverse,
    ggplot2
)

# import data ----
META <- read_csv("data/eegap/META.csv") %>% select(-LOC) %>% unique()
EMRI_RAW <- read_csv("data/eegap/EMRI.csv")
EMRI <- EMRI_RAW %>% 
    pivot_wider(
        names_from = "PARAM",
        values_from = "VALUE"
    ) %>% 
    mutate(
        ADIPOSITY_INDEX = Fat / Lean,
        FAT_DELTA = (Fat - Fat[DATE=="2024-07-02"]),
        LEAN_DELTA = (Lean - Lean[DATE=="2024-07-02"]),
        ADIPOSITY_INDEX_DELTA = FAT_DELTA / LEAN_DELTA
    ) %>% 
    drop_na() %>% 
    pivot_longer(
        cols = 4:11,
        names_to = "PARAM",
        values_to = "VALUE"
    ) %>% 
    left_join(META, by = c("ID")) %>%
    mutate(
        ID = as.factor(ID),
        PARAM = as.factor(PARAM),
        DIET = as.factor(DIET),
        SEX = as.factor(SEX)
    )

# all PARAM ~ DIET

mdl_data <- EMRI %>% 
    group_by(PARAM) %>% 
    group_split()

mdls <- mdl_data %>% 
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

