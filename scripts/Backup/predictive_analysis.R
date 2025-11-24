pacman::p_load(
    tidyverse,
    ggplot2,
    mmand
)

# set current source file location to here
setwd(this.path::here())

# get single animal data set
data <- read_csv("../data/sable/202412310738-extd_m.csv")

# I want to get feeding bouts

# first take food hopper vector
food_hopper <- data %>% 
    select(
        DateTime, FoodA_1
    ) %>% 
    mutate(
        DateTime = lubridate::mdy_hms(DateTime)
    )
food_hopper

food_hopper %>% 
    ggplot(aes(
        DateTime, FoodA_1
    )) +
    geom_line()

# do morphological operation 
# opening gets rid of upwards peaks
# closing gets rid of downwards peaks
# I smoothed out the curve to get rid of small jitters
# finally I correct all measurements that are not strictly going down

food_hopper_opening <- food_hopper %>% 
    mutate(
        opening = mmand::gaussianSmooth(FoodA_1, 61) %>% 
            mmand::opening(., rep(1, 121)) %>% 
            mmand::closing(., rep(1, 601)),
        smoothed_curve = if_else(opening>cummin(opening), NA_real_, opening)
    ) %>% 
    fill(smoothed_curve) %>% 
    mutate(
        bouts = if_else(smoothed_curve<lag(smoothed_curve), "bout", "stable")
    )
food_hopper_opening

# get only bouts data
bouts_data <- food_hopper_opening %>% 
    drop_na() %>% 
    mutate(bout_tag = data.table::rleid(bouts)) %>% 
    filter(bouts == "bout") %>% 
    mutate(
        bout_tag = as.numeric(as.factor(bout_tag))
    )
bouts_data

# grab only the timestamps

bout_stamps <- bouts_data %>% 
    group_by(bout_tag) %>% 
    summarise(
        init = min(DateTime),
        end = max(DateTime),
        centroid = ((end - init)/2) + init,
        bout_tag = unique(bout_tag)
    )
bout_stamps

# RQ data
rq <- data %>% 
    select(
        DateTime, RQ_1
    ) %>% 
    mutate(
        DateTime = lubridate::mdy_hms(DateTime)
    )

# lets prepare another metric
kcal <- data %>% 
    select(
        DateTime, kcal_hr_1
    ) %>% 
    mutate(
        DateTime = lubridate::mdy_hms(DateTime)
    )
kcal

kcal %>% 
    filter(kcal_hr_1 < 1) %>% 
    ggplot(aes(
        DateTime, kcal_hr_1
    )) +
    geom_line() +
    geom_vline(xintercept = bout_stamps$centroid, color = "red")

rq %>% 
    filter(RQ_1 < 1) %>% 
    ggplot(aes(
        DateTime, (RQ_1)
    )) +
    geom_point() +
    geom_line()

# map all the centroids and get a time window of 1 minute

time_win <- bout_stamps %>% 
    group_by(row_number()) %>% 
    group_split %>% 
    map_dfr(., function(X){
        centroid <- X$centroid
        bout_dur <- as.numeric(X$end - X$init)
        init <- centroid - lubridate::minutes(5)
        end <- centroid + lubridate::minutes(5)
        bout_tag <- X$bout_tag[1]
        data_out <- kcal %>% 
            filter(
                DateTime >= init,
                DateTime <= end
            ) %>% 
            mutate(bout_tag = bout_tag,
                   centroid = centroid,
                   z_rq = scale(kcal_hr_1),
                   bout_len = bout_dur)
        return(data_out)
    })

test_win <- time_win %>% 
    group_by(bout_tag) %>% 
    mutate(
        rel_time = DateTime - centroid,
        bf = if_else(rel_time<0, "before", "after")
    )
test_win

mdl <- lm(
    data = test_win,
    kcal_hr_1 ~ bf * bout_len
)
summary(mdl)

emmeans::emmeans(
    mdl,
    pairwise ~ bf | bout_len,
    at = list(bout_len = c(seq(0, 100, 1))),
    type = "response"
)$contrasts %>% 
    broom::tidy() %>% 
    ggplot(aes(
        bout_len, estimate, ymin = estimate-std.error, ymax = estimate+std.error
    )) +
    geom_line() +
    geom_errorbar()

time_win %>% 
    group_by(bout_tag) %>% 
    mutate(
        rel_time = DateTime - centroid
    ) %>% 
    ggplot(aes(
        rel_time, z_rq, group = bout_tag
    )) +
    geom_line() +
    geom_smooth(method = "loess", color = "red", aes(group=1))


# get some idea of how long bouts typically are

bouts_data %>% 
    group_by(bout_tag) %>% 
    summarise(
        len = n()
    ) %>%
    ggplot(aes(len)) +
    geom_histogram()



bouts_data %>% 
    ggplot(aes(
        DateTime, smoothed_curve, group = 1
    )) +
    geom_line(aes(color = bout_tag == 41))


food_hopper_opening %>% 
    ggplot(aes(
        DateTime, smoothed_curve, color = bouts, group = 1
    )) +
    geom_line()


