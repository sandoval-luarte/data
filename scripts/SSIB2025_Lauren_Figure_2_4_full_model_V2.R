pacman::p_load(
    tidyverse,
    ggplot2,
    lmtest,
    caret,
    ggpubr,
    latex2exp,
    patchwork
)

# change the directory to source file location
setwd(this.path::here())

# main data ----
dat <- read_csv("parameter_compilate.csv") %>% 
    mutate(
        ID = as.factor(ID)
    ) %>% 
  drop_na()

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
    mutate(across(starting_bw:tee_delta, scale)) %>% 
  drop_na()

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
                  data = dat_scaled %>% mutate(null=rnorm(22,0,1)), method="lm",
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

results_option1 <- resamples(list("Null model"=mdl_0_cv,
                                  "Baseline values"=mdl_1_cv,
                                  "ROCV"=mdl_2_cv,
                                  "Full model"=mdl_3_cv,
                                  "Baseline BW"=mdl_4_cv,
                                  "Baseline FI"=mdl_5_cv,
                                  "Baseline AI"=mdl_6_cv,
                                  "Baseline fat"=mdl_7_cv,
                                  "Baseline lean"=mdl_8_cv,
                                  "Δ FI"=mdl_9_cv,
                                  "Δ AI"=mdl_10_cv,
                                  "Δ fat"=mdl_11_cv,
                                  "Δ lean"=mdl_12_cv,
                                  "ΔTEE"=mdl_13_cv
))
dotplot(results_option1)

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


colorblind_m <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#F0E442")

#starting body weight ----

p1 <- dat %>% 
  ggplot(aes(x=starting_bw) )+ 
  geom_density(fill=colorblind_m[1])+
  theme_pubr()
p1

#starting food intake ----
p2<- dat %>% 
  ggplot(aes(x=starting_fi) )+ 
  geom_density(fill=colorblind_m[2])+
  theme_pubr()
p2

#starting adiposity index ----
p3<- dat %>% 
  ggplot(aes(x=starting_ai) )+ 
  geom_density(fill=colorblind_m[3])+
  theme_pubr()
p3

#starting fat mass ----
p4<- dat %>% 
  ggplot(aes(x=starting_fat) )+ 
  geom_density(fill=colorblind_m[4])+
  theme_pubr()
p4

#starting lean mass ----
p5<- dat %>% 
  ggplot(aes(x=starting_lean) )+ 
  geom_density(fill=colorblind_m[5])+
  theme_pubr()
p5

#time FI ----
p6<- dat %>% 
  ggplot(aes(x=time_fi) )+ 
  geom_density(fill=colorblind_m[2])+
  theme_pubr() +
  xlab(latex2exp::TeX(r"($\beta_{Food \ intake}$)"))+
  ylab("Density")
p6 

#time AI ----
p7<- dat %>% 
  ggplot(aes(x=time_ai) )+ 
  geom_density(fill=colorblind_m[3])+
  theme_pubr() +
  xlab(latex2exp::TeX(r"($\beta_{Adiposity \ index}$)"))+
  ylab("Density")
p7 

#time fat ----
p8<- dat %>% 
  ggplot(aes(x=time_fat) )+ 
  geom_density(fill=colorblind_m[4])+
  theme_pubr() +
  xlab(latex2exp::TeX(r"($\beta_{Fat}$)"))+
  ylab("Density")
p8 

#time lean ----
p9<- dat %>% 
  ggplot(aes(x=time_lean) )+ 
  geom_density(fill=colorblind_m[5])+
  theme_pubr() +
  xlab(latex2exp::TeX(r"($\beta_{Lean}$)"))+
  ylab("Density")
p9 

#time EE ----
p10<- dat %>% 
  ggplot(aes(x=tee_delta) )+ 
  geom_density(fill=colorblind_m[7],adjust=5)+
  theme_pubr() +
  xlab(latex2exp::TeX(r"($\beta_{Tee}$)"))+
  ylab("Density")
p10 

#arrange FIGURE 2 Lauren SSIB 2025----


p1 <- p1 + labs(x = "Body weight (g) at baseline", y="Density")
p2 <- p2 + labs(x = "Food intake (kcal) at baseline",y="Density") 
p5 <- p5 + labs(x = "Lean mass (g) at baseline",y="Density") 
p4 <- p4 + labs(x = "Fat mass (g) at baseline",y="Density") 
p3 <- p3 + labs(x = "Adiposity index at baseline",y="Density") 


col_1 <- (p1 + p10 + p2 + p6 +p5 + p9 + p4 + p8 + p3 + p7)+
  plot_layout(ncol = 2, nrow=5) +
  plot_annotation(tag_levels = list(c("A","B","C","D","E","F","G","H","I","J")))
col_1 

