pacman::p_load(
    ggplot2,
    tidyverse
)

# run this line to select .csv file
dat <- read_csv(file.choose()) %>% 
    mutate(
        DateTime = lubridate::mdy_hms(DateTime)
    )

# bodyweight check ----
bw <- dat %>% 
    select(DateTime, contains("BodyMass_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

# mean mass
bw %>% 
    group_by(name) %>% 
    summarise(
        mean_weight = mean(value)
    )

# plot
bw %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~name)

# water check ----
water <- dat %>% 
    select(DateTime, contains("Water_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

# plot
water %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~name, scales="free")

# locomotor check ----
locomotor <- dat %>% 
    select(DateTime, contains("AllMeters_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

locomotor %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_line() +
    facet_wrap(~name)

# tee check ----
tee <- dat %>% 
    select(DateTime, contains("kcal_hr_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

tee %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_line() +
    facet_wrap(~name)

# O2 check ----
VO2 <- dat %>% 
    select(DateTime, contains("VO2_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

VO2 %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_point() +
    facet_wrap(~name)

# co2 check ----
VCO2 <- dat %>% 
    select(DateTime, contains("VCO2_")) %>% 
    pivot_longer(cols=-DateTime) %>% 
    filter(
        lubridate::second(DateTime) == 0
    )

VCO2 %>% 
    ggplot(aes(
        DateTime, value
    )) +
    geom_line() +
    facet_wrap(~name)

# weir equation ----
wdat <- dat %>% 
    filter(
        lubridate::second(DateTime) == 0
    ) %>% 
    select(DateTime, contains("VO2_"), contains("VCO2_")) %>% 
    pivot_longer(-DateTime) %>% 
    mutate(ID = str_extract(pattern = "_[0-9]+", name),
           name = str_extract(pattern = "[[:alnum:]]+_", name),
           value = if_else(ID=="_5" & name=="VO2_", value+2, value)) %>% 
    ungroup() %>% 
    group_by(ID, DateTime) %>% 
    mutate(
        weir_eq = (1.44 * (3.94 * value[name=="VO2_"] 
                          + 1.11 * value[name=="VCO2_"]))/24
    )

wdat %>% 
    ggplot(aes(
        DateTime, weir_eq
    )) +
    geom_point() +
    geom_line() +
    facet_wrap(~ID)
