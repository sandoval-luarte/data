pacman::p_load(
    tidyverse,
    ggplot2,
    lmtest,
    caret
)

# change the directory to source file location
setwd(this.path::here())

# main data ----
dat <- read_csv("parameter_compilate.csv") %>% 
    mutate(
        ID = as.factor(ID)
    )

# all single variable models ----

# set data to long format
dat_long <- dat %>% 
    pivot_longer(
        cols = -c(ID, tee_delta)
    )
dat_long

single_var_mdls <- dat_long %>% 
    group_by(name) %>% 
    group_split() %>% 
    map_dfr(., function(dat){
        mdl <- lm(data=dat %>%
                      mutate(value=(value)), tee_delta ~ value)
        mdl_tbl <- mdl %>% 
            broom::tidy(conf.int = TRUE) %>% 
            mutate(var = dat$name[1]) %>% 
            filter(term=="value")
        return(mdl_tbl)
    })

p1 <- single_var_mdls %>% 
    ggplot(aes(
        var, estimate,
        ymin = conf.low, ymax = conf.high
    )) +
    geom_hline(yintercept = 0) +
    geom_pointrange() +
    theme(axis.text.x = element_text(angle = 90))
p1

# models ----

dat_scaled <- dat %>% 
    mutate(across(starting_bw:tee_delta, scale))

# starting body composition model ----
mdl_1 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ starting_bw * starting_lean * starting_fat
)
summary(mdl_1)

# development of body composition model ----
mdl_2 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ time_bw * time_lean * time_fat
)
summary(mdl_2)

# body composition start + development ----
mdl_3 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ time_bw * time_lean * time_fat +
        starting_bw + starting_lean + starting_fat
)
summary(mdl_3)

# intake starting ----
mdl_4 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ starting_fi 
)
summary(mdl_4)

# intake development ----
mdl_5 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ time_fi 
)
summary(mdl_5)

# intake start + development ----
mdl_6 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ time_fi + starting_fi
)
summary(mdl_6)

# null model ----
mdl_0 <- lm(
    data = dat_scaled %>% select(-ID),
    tee_delta ~ 1
)
summary(mdl_0)

# cv ----

ctrl <- trainControl(method = "repeatedcv", number=5, repeats=100)

# starting values
mdl_0_cv <- train(time_bw ~ null,
                  data = dat_scaled %>% mutate(null=rnorm(23,0,1)), method="lm",
                  trControl=ctrl)
mdl_1_cv <- train(time_bw ~ starting_bw + starting_fi + starting_ai +
                      starting_ai + starting_fat + starting_lean,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
# time course values
mdl_2_cv <- train(time_bw ~ time_fi + time_ai + time_fat + time_lean +
                      tee_delta,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
# full model
mdl_3_cv <- train(time_bw ~ time_fi + time_ai + time_fat + time_lean +
                      tee_delta + starting_bw + starting_fi + starting_ai +
                      starting_fat + starting_lean,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
# individual models
mdl_4_cv <- train(time_bw ~ starting_bw,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_5_cv <- train(time_bw ~ starting_fi,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_6_cv <- train(time_bw ~ starting_ai,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_7_cv <- train(time_bw ~ starting_fat,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_8_cv <- train(time_bw ~ starting_lean,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_9_cv <- train(time_bw ~ time_fi,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_10_cv <- train(time_bw ~ time_ai,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_11_cv <- train(time_bw ~ time_fat,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_12_cv <- train(time_bw ~ time_lean,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)
mdl_13_cv <- train(time_bw ~ tee_delta,
                  data = dat_scaled, method="lm",
                  trControl=ctrl)


results <- resamples(list(null_model=mdl_0_cv,
                          starting_model=mdl_1_cv,
                          timecourse_model=mdl_2_cv,
                          full_model=mdl_3_cv,
                          starting_bw=mdl_4_cv,
                          starting_fi=mdl_5_cv,
                          starting_ai=mdl_6_cv,
                          starting_fat=mdl_7_cv,
                          starting_lean=mdl_8_cv,
                          time_fi=mdl_9_cv,
                          time_ai=mdl_10_cv,
                          time_fat=mdl_11_cv,
                          time_lean=mdl_12_cv,
                          tee_delta=mdl_13_cv
                          ))
dotplot(results)
# compare best model
timecourse_model <- mdl_2_cv$finalModel
null_model <- mdl_0_cv$finalModel
tee_model <- mdl_13_cv$finalModel

anova(null_model, timecourse_model) #null vs all variables
anova(null_model, tee_model)  #null vs tee
anova(tee_model, timecourse_model)


# analysis of best model
coeffs <- summary(timecourse_model) %>% 
    broom::tidy(conf.int=TRUE)
coeffs

coeffs %>% 
    ggplot(aes(
        term, estimate,
        ymin = conf.low, ymax = conf.high
    )) +
    geom_hline(yintercept = 0) +
    geom_pointrange()

best_mdl_emm <- emmeans::emmeans(
    timecourse_model,
    ~ time_fi + time_ai + time_fat + time_lean + tee_delta,
    type = "response"
)
best_mdl_emm
