pacman::p_load(
    tidyverse,
    ggplot2,
    clipr,
    bigsnpr,
    distributions3
)

# compute the entropy of a beta distribution
differential_entropy <- function(alpha, beta){
        tryCatch(-integrate(function(x) { 
            dbeta(x, alpha, beta) * log(dbeta(x, alpha, beta)) 
        }, lower = 0, upper = 1)$value, error=function(err) NA)
}

# sample from a set of delays to compute the total amount of time to get 150 pellets
# given 15 sec, 2250 seconds per day / 0.83 hours is what im aiming for
get_time <- function(delay_vector, probability_vector){
    sample(x = delay_vector, size = 1000, prob = probability_vector)
}

# grid for the beta distribution parameters
beta_params <- expand_grid(
    shape1 = seq(0, 1, 0.01),
    shape2 = seq(0, 1, 0.01)
)
beta_params

# compute entropy and time for all the distributions
beta_entropy <- beta_params %>% 
    group_split(row_number()) %>% 
    map_dfr(., function(X){
        expected_delay <- (X$shape1 / (X$shape1+X$shape2))*150*150
        expected_entropy <- differential_entropy(X$shape1, X$shape2)
        return(tibble(
            alpha = X$shape1,
            beta = X$shape2,
            expected_delay = expected_delay,
            expected_entropy = expected_entropy
        ))
    })

# get those with an expected delay of 3000
sub_beta_entropy <- beta_entropy %>% 
    filter(expected_delay<3000, expected_delay>2000) %>% 
    slice_max(order_by = expected_entropy, n=1)
sub_beta_entropy

rbeta(shape1 = sub_beta_entropy$alpha, shape2 = sub_beta_entropy$beta,
      n=1000) %>%  scales::rescale(., to = c(0,150)) %>% 
    trunc()

table((rbeta(shape1 = sub_beta_entropy$alpha, shape2 = sub_beta_entropy$beta,
           n=150)) %>% 
    scales::rescale(., to = c(0,150)) %>% 
        cut(., breaks=seq(0,150,15)))

