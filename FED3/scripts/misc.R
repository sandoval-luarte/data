pacman::p_load(
    tidyverse,
    ggplot2
)

setwd(this.path::here())

file_path <- list.files(
    path = "~/repos_sync/data/FED3/data",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "*.CSV"
)
file_path

dat <- file_path %>% 
    map_dfr(., function(X){
        dat <- read_csv(X, show_col_type = FALSE, col_types = cols(.default = col_character()))
        COLS <- c(
                "MM:DD:YYYY hh:mm:ss",
                "Library_Version",
                "Session_type",
                "Device_Number",
                "Battery_Voltage",
                "Motor_Turns",
                "FR",
                "Event",
                "Active_Poke",
                "Left_Poke_Count",
                "Right_Poke_Count",
                "Failed_Retrieval_Count",
                "Pellet_Count",
                "Block_Pellet_Count",
                "Retrieval_Time",
                "InterPelletInterval",
                "Poke_Time",
                "Delay"
            )
        col_match <- sum(colnames(dat) %in% COLS)
        if (col_match == 18){
            return(dat)
        }
        else if (col_match == 17){
            return(dat %>% mutate(Delay = "1"))
        }
        else {
            return(tibble::tibble(!!!COLS, .rows = 0, .name_repair = ~COLS))
        }
    })

dat_s1 <- dat %>% 
    rename(ymd_hms = `MM:DD:YYYY hh:mm:ss`) %>% 
    mutate(
        ymd_hms = lubridate::mdy_hms(ymd_hms),
        Pellet_Count = as.numeric(Pellet_Count)
    ) %>% 
    filter(lubridate::year(ymd_hms) == 2025)

case_1 <- dat_s1 %>% 
    filter(Session_type %in% c("certainty", "uncertainty")) %>% 
    group_by(Device_Number) %>% 
    mutate(
        Poke_Time = replace_na(as.numeric(Poke_Time), 0),
        Delay = if_else(Delay == "nan", NA, Delay) %>% as.numeric(),
        Retrieval_Time = if_else(Retrieval_Time == "nan", "-1", Retrieval_Time),
        Delay = if_else(Retrieval_Time=="-1", NA_real_, Delay)
    ) %>% 
    group_by(Device_Number, Session_type) %>% 
    arrange(ymd_hms, .by_group = TRUE) %>% 
    mutate(
        Delay_group = rev(cumsum(rev(!is.na(Delay))))
    ) %>% 
    fill(Delay, .direction = c("up"))

fs <- case_1 %>% 
    group_by(Device_Number, Delay_group, Delay, Session_type,
             date = lubridate::date(ymd_hms)) %>% 
    summarise(
        fs = n(),
        fs_time = sum(Poke_Time),
        rel_delay = log(abs(Delay - 15)+1),
        ret_time = max(as.numeric(Retrieval_Time))
    ) %>% 
    ungroup() %>% 
    group_by(Device_Number) %>% 
    mutate(
        date = as.numeric(as.factor(date))
    ) %>% 
    ungroup() %>% 
    group_by(Device_Number, Delay_group, Session_type, date, Delay, rel_delay) %>% 
    summarise(
        fs_poke = n(),
        fs_time = max(fs_time),
        ret_time = max(ret_time)
    )

mdl <- lme4::glmer.nb(
    data = fs %>% filter(date!=4, fs_poke>0),
    fs_poke ~ Session_type*date*rel_delay + (1|Device_Number)
)
summary(mdl)

fs %>% 
    filter() %>% 
    ggplot(aes(
        rel_delay, fs_poke, group = Device_Number, color = Session_type
    )) +
    geom_point() +
    geom_vline(xintercept = log(15)) +
    geom_smooth(aes(group=interaction(Session_type, date))) +
    scale_y_continuous(transform = "log") +
    scale_x_continuous(transform = "log") +
    facet_wrap(~date)

fs %>% 
    ggplot(aes(Delay)) +
    geom_histogram() +
    facet_wrap(~Session_type)

emmeans::emmeans(
    mdl,
    pairwise ~ Session_type*date|rel_delay,
    at = list(rel_delay = seq(0, 300, 50)),
    type = "link"
) %>%
    {.$emmeans} %>% 
    broom::tidy() %>% 
    ggplot(aes(rel_delay, estimate, group=Session_type, color=Session_type)) +
    geom_line()

fs %>% 
    ungroup() %>% 
    group_by(Device_Number, Session_type, date) %>% 
    summarise(
        total_fs = mean(fs),
        total_fs_time = mean(fs_time),
        pellets = length(unique(Delay_group)),
        rel_total_fs = total_fs / pellets
    ) %>% 
    ggplot(aes(
        Session_type, total_fs
    )) +
    geom_point() +
    facet_wrap(~date)

fs %>% 
    ggplot(aes(
        rel_delay+1, fs, color = Device_Number
    )) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_continuous(transform = "log") +
    scale_x_continuous(transform = "log") +
    facet_wrap(~Device_Number * Session_type, scales = "free")

