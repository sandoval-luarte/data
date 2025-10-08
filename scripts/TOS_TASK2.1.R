#TOS 2025 LL CS TASK 2.1 and 2.2
#AIM task 2.1> Build deviations using modeling (LMER), must fit a linear mixed-effects model (LMM) controlling for sex, strain, and diet, with random intercepts by animal (ID).
#AIM task 2.2> evaluate goodness of fit of each model (generates an AIC ranking for all models using the best per y predicted)

#Libraries####
library(lmerTest)
library(readr)


#Y1####
y1_data <- read_csv("../data/y1_data.csv")
y1_model <- y1_data %>%
  select(ID, STRAIN, SEX, y1) %>%
  left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
  distinct()

# Simple linear models instead of mixed models because I dont have more than one obs per ID
m1 <- lm(y1 ~ 1, data = y1_model)#Apparently this is the best model
m2 <- lm(y1 ~ SEX, data = y1_model)
m3 <- lm(y1 ~ STRAIN, data = y1_model) 
m4 <- lm(y1 ~ DIET_FORMULA, data = y1_model)
m5 <- lm(y1 ~ SEX * STRAIN, data = y1_model)
m6 <- lm(y1 ~ SEX * DIET_FORMULA, data = y1_model)
m7 <- lm(y1 ~ STRAIN * DIET_FORMULA, data = y1_model)
m8 <- lm(y1 ~ SEX * STRAIN * DIET_FORMULA, data = y1_model)

#compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
#o term improves the model significantly, the simplest model (Model 1: y1 ~ 1) is the best according to parsimony and significance tests.

#Y2####
y2_data <- read_csv("../data/y2_data.csv")
y2_model <- y2_data %>%
  select(ID, STRAIN, SEX, y2) %>%
  left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
  distinct()

# Simple linear models for y2
m1 <- lm(y2 ~ 1, data = y2_model) #Apparently this is the best model
m2 <- lm(y2 ~ SEX, data = y2_model)
m3 <- lm(y2 ~ STRAIN, data = y2_model) 
m4 <- lm(y2 ~ DIET_FORMULA, data = y2_model)
m5 <- lm(y2 ~ SEX * STRAIN, data = y2_model)
m6 <- lm(y2 ~ SEX * DIET_FORMULA, data = y2_model)
m7 <- lm(y2 ~ STRAIN * DIET_FORMULA, data = y2_model)
m8 <- lm(y2 ~ SEX * STRAIN * DIET_FORMULA, data = y2_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
#Does this means adiposity index in our dataset doesn’t differ significantly by SEX, STRAIN, DIET_FORMULA, or their interactions?

#Y3####
y3_data <- read_csv("../data/y3_data.csv")
y3_model <- y3_data %>%
  select(ID, STRAIN, SEX, y3) %>%
  left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
  distinct()

# Simple linear models for y3
m1 <- lm(y3 ~ 1, data = y3_model) #Apparently this is the best model
m2 <- lm(y3 ~ SEX, data = y3_model)
m3 <- lm(y3 ~ STRAIN, data = y3_model)
m4 <- lm(y3 ~ DIET_FORMULA, data = y3_model)
m5 <- lm(y3 ~ SEX * STRAIN, data = y3_model)
m6 <- lm(y3 ~ SEX * DIET_FORMULA, data = y3_model)
m7 <- lm(y3 ~ STRAIN * DIET_FORMULA, data = y3_model)
m8 <- lm(y3 ~ SEX * STRAIN * DIET_FORMULA, data = y3_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
#The rate of body weight change during the regain phase does not differ significantly by SEX, STRAIN, DIET_FORMULA, or their interactions in this dataset.

#Y4####
y4_data <- read_csv("../data/y4_data.csv")
y4_model <- y4_data %>%
  select(ID, STRAIN, SEX, y4) %>%
  left_join(select(echoMRI_data, ID, DIET_FORMULA), by = "ID") %>%
  distinct()

# Simple linear models for y4
m1 <- lm(y4 ~ 1, data = y4_model) 
m2 <- lm(y4 ~ SEX, data = y4_model) #Apparently this is the best model
m3 <- lm(y4 ~ STRAIN, data = y4_model)
m4 <- lm(y4 ~ DIET_FORMULA, data = y4_model)
m5 <- lm(y4 ~ SEX * STRAIN, data = y4_model)
m6 <- lm(y4 ~ SEX * DIET_FORMULA, data = y4_model)
m7 <- lm(y4 ~ STRAIN * DIET_FORMULA, data = y4_model)
m8 <- lm(y4 ~ SEX * STRAIN * DIET_FORMULA, data = y4_model)

# Compare models
anova(m1, m2, m3, m4, m5, m6, m7, m8)
#SEX effect has p = 0.075, which is not significant at 0.05, but it’s borderline
