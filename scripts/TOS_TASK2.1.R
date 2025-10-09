# TOS 2025 LL CS TASK 2.1 and 2.2
# AIM task 2.1> Build deviations using modeling (LMER), must fit a linear mixed-effects model (LMM) controlling for sex, strain, and diet, with random intercepts by animal (ID).
# AIM task 2.2> evaluate goodness of fit of each model (generates an AIC ranking for all models using the best per y predicted)

# Libraries####
library(lmerTest)
library(readr)
library(caret)
setwd(this.path::here())


# Y1####
echoMRI_data <- read_csv("../data/echomri.csv")
y1_data <- read_csv("../data/y1_data.csv")
y1_model <- y1_data %>%
    select(ID, STRAIN, SEX, y1) %>%
    left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
    distinct()

# Simple linear models instead of mixed models because I dont have more than one obs per ID
m1_y1 <- lm(y1 ~ 1, data = y1_model) # Apparently this is the best model
m2_y1 <- lm(y1 ~ SEX, data = y1_model)
m3_y1 <- lm(y1 ~ STRAIN, data = y1_model)
m4_y1 <- lm(y1 ~ DIET_FORMULA, data = y1_model)
m5_y1 <- lm(y1 ~ SEX * STRAIN, data = y1_model)
m6_y1 <- lm(y1 ~ SEX * DIET_FORMULA, data = y1_model)
m7_y1 <- lm(y1 ~ STRAIN * DIET_FORMULA, data = y1_model)
m8_y1 <- lm(y1 ~ SEX * STRAIN * DIET_FORMULA, data = y1_model)

# compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
# o term improves the model significantly, the simplest model (Model 1: y1 ~ 1) is the best according to parsimony and significance tests.

# Y2####
y2_data <- read_csv("../data/y2_data.csv")
y2_model <- y2_data %>%
    select(ID, STRAIN, SEX, y2) %>%
    left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
    distinct()

# Simple linear models for y2
m1_y2 <- lm(y2 ~ 1, data = y2_model) # Apparently this is the best model
m2_y2 <- lm(y2 ~ SEX, data = y2_model)
m3_y2 <- lm(y2 ~ STRAIN, data = y2_model)
m4_y2 <- lm(y2 ~ DIET_FORMULA, data = y2_model)
m5_y2 <- lm(y2 ~ SEX * STRAIN, data = y2_model)
m6_y2 <- lm(y2 ~ SEX * DIET_FORMULA, data = y2_model)
m7_y2 <- lm(y2 ~ STRAIN * DIET_FORMULA, data = y2_model)
m8_y2 <- lm(y2 ~ SEX * STRAIN * DIET_FORMULA, data = y2_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
# Does this means adiposity index in our dataset doesn’t differ significantly by SEX, STRAIN, DIET_FORMULA, or their interactions?

# Y3####
y3_data <- read_csv("../data/y3_data.csv")
y3_model <- y3_data %>%
    select(ID, STRAIN, SEX, y3) %>%
    left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
    distinct()

# Simple linear models for y3
m1_y3 <- lm(y3 ~ 1, data = y3_model) # Apparently this is the best model
m2_y3 <- lm(y3 ~ SEX, data = y3_model)
m3_y3 <- lm(y3 ~ STRAIN, data = y3_model)
m4_y3 <- lm(y3 ~ DIET_FORMULA, data = y3_model)
m5_y3 <- lm(y3 ~ SEX * STRAIN, data = y3_model)
m6_y3 <- lm(y3 ~ SEX * DIET_FORMULA, data = y3_model)
m7_y3 <- lm(y3 ~ STRAIN * DIET_FORMULA, data = y3_model)
m8_y3 <- lm(y3 ~ SEX * STRAIN * DIET_FORMULA, data = y3_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# The rate of body weight change during the regain phase does not differ significantly by SEX, STRAIN, DIET_FORMULA, or their interactions in this dataset.

# Y4####
y4_data <- read_csv("../data/y4_data.csv")
y4_model <- y4_data %>%
    select(ID, STRAIN, SEX, y4) %>%
    left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
    distinct()

# Simple linear models for y4
m1_y4 <- lm(y4 ~ 1, data = y4_model)
m2_y4 <- lm(y4 ~ SEX, data = y4_model) # Apparently this is the best model
m3_y4 <- lm(y4 ~ STRAIN, data = y4_model)
m4_y4 <- lm(y4 ~ DIET_FORMULA, data = y4_model)
m5_y4 <- lm(y4 ~ SEX * STRAIN, data = y4_model)
m6_y4 <- lm(y4 ~ SEX * DIET_FORMULA, data = y4_model)
m7_y4 <- lm(y4 ~ STRAIN * DIET_FORMULA, data = y4_model)
m8_y4 <- lm(y4 ~ SEX * STRAIN * DIET_FORMULA, data = y4_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
# SEX effect has p = 0.075, which is not significant at 0.05, but it’s borderline

# aic comparisons
aic_data <- AIC(
    m1_y1, m2_y1, m3_y1, m4_y1, m5_y1, m6_y1, m7_y1, m8_y1,
    m1_y2, m2_y2, m3_y2, m4_y2, m5_y2, m6_y2, m7_y2, m8_y2,
    m1_y3, m2_y3, m3_y3, m4_y3, m5_y3, m6_y3, m7_y3, m8_y3,
    m1_y4, m2_y4, m3_y4, m4_y4, m5_y4, m6_y4, m7_y4, m8_y4
) %>%
    as_tibble(rownames = "model_name") %>%
    mutate(
        parameter = str_extract(string = model_name, pattern = "y[0-9]"),
        model_number = str_extract(string = model_name, pattern = "m[0-9]")
    ) %>%
    group_by(parameter) %>%
    mutate(
        aic_norm = scale(AIC)
    )
aic_data

p1 <- aic_data %>%
    ggplot(aes(
        parameter, model_number,
        fill = aic_norm
    )) +
    geom_tile(
        color = "white",
        lwd = 1.5,
        linetype = 1
    ) +
    geom_text(aes(label = round(AIC, 2)), color = "black", size = 4) +
    coord_fixed() +
    scale_fill_gradient(low = "red", high = "white")
p1


# kfold RMSE-based best model selection -----

## generate data ----
## dependent variables are scales (zscores)
complete_data <- y1_data %>%
    left_join(select(y2_data, ID, y2), by = "ID") %>%
    left_join(select(y3_data, ID, y3), by = "ID") %>%
    left_join(select(y4_data, ID, y4), by = "ID") %>%
    left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
    distinct() %>%
    mutate(
        across(y1:y4, ~ scale(.x))
    )
complete_data

## select best model ----

set.seed(42)
candidate_models <- list(
    "model_1" = cbind(y1, y2, y3, y4) ~ 1,
    "model_2" = cbind(y1, y2, y3, y4) ~ SEX,
    "model_3" = cbind(y1, y2, y3, y4) ~ STRAIN,
    "model_4" = cbind(y1, y2, y3, y4) ~ DIET_FORMULA,
    "model_5" = cbind(y1, y2, y3, y4) ~ SEX * STRAIN,
    "model_6" = cbind(y1, y2, y3, y4) ~ SEX * DIET_FORMULA,
    "model_7" = cbind(y1, y2, y3, y4) ~ STRAIN * DIET_FORMULA,
    "model_8" = cbind(y1, y2, y3, y4) ~ SEX * STRAIN * DIET_FORMULA
)

k_folds <- createFolds(
    complete_data$y1,
    k = 5, list = TRUE, returnTrain = FALSE
)

results_log <- tibble()

for (model_name in names(candidate_models)) {
    formula <- candidate_models[[model_name]]
    fold_scores <- c()

    for (i in 1:5) {
        validation_indices <- k_folds[[i]]
        train_data <- complete_data[-validation_indices, ]
        validation_data <- complete_data[validation_indices, ]

        model <- lm(formula, data = train_data)

        predictions <- predict(model, newdata = validation_data)

        actuals <- validation_data %>%
            select(y1, y2, y3, y4)

        rmses <- map2_dbl(
            actuals,
            as.data.frame(predictions),
            ~ sqrt(mean(mean(.x - .y)^2))
        )

        mean_scaled_rmse <- mean(rmses)

        fold_scores <- c(fold_scores, mean_scaled_rmse)
    }

    results_log <- bind_rows(
        results_log,
        tibble(
            model = model_name,
            mean_rmse = mean(fold_scores),
            std_rmse = sd(fold_scores)
        )
    )
}

final_results <- results_log %>%
    arrange(mean_rmse)
final_results

winner_name <- final_results$model[1]
winning_formula <- candidate_models[[winner_name]]
winning_formula

## get residuals with best model ----

std_residuals_log <- tibble()

target_vars <- c("y1", "y2", "y3", "y4")

for (i in 1:5) {
    validation_indices <- k_folds[[i]]
    train_data <- complete_data[-validation_indices, ]
    validation_data <- complete_data[validation_indices, ]

    model <- lm(winning_formula, data = train_data)

    sd_training_residuals <- apply(residuals(model), 2, sd)

    predictions_validation <- predict(model, newdata = validation_data)
    actuals_validation <- validation_data %>% select(all_of(target_vars))

    raw_residuals_validation <- actuals_validation - predictions_validation

    std_residuals_validation <- sweep(raw_residuals_validation, 2, sd_training_residuals, FUN = "/")

    std_residuals_validation <- std_residuals_validation %>%
        as_tibble() %>%
        rename_with(~ paste0(.x, "_std_resid")) %>%
        mutate(
            ID = validation_data$ID
        )

    std_residuals_log <- bind_rows(std_residuals_log, std_residuals_validation)
}

final_df <- complete_data %>%
    left_join(
        std_residuals_log %>%
            arrange(ID),
        by = "ID"
    )
