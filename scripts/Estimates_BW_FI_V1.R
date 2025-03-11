
##Readme####
#We aim to obtain the slope of BW and cumulative FI among days. 
#each slope (X= days, Y=Body weight) will represent the rate of weight gain per day for that ID
#each slope (X= days, Y=Cumulative food intake in kcal) will represent the rate of food intake over time (e.g., kcal per day).

##Description####
#The animals used for this analysis correspond to C57BL6J males and females fed with LFD and HFD
#as well as NZO females fed with LFD. 
#The HFD formula is D12451i, research diets and LFD formula is D12450Ki, research diets

##Libraries####
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(writexl)


#Body weight (g) among days for C57 and NZO mice ####

##Data####
bw <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(2,3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study
  group_by(ID) %>% 
  filter(DATE < "2025-02-24") %>% #limited the analysis priod to start the food restriction 
  mutate(rel_days= DATE - min(DATE))
##Graph####
bw %>% 
  ggplot(aes(x = rel_days, y = BW)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  labs(
    x = "Date",
    y = "Body Weight (g)"
  ) +  # White background
  theme_classic()  +
  facet_grid(SEX~STRAIN)+
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) 
##Main model ####
slopes <- bw %>%
  group_by(ID) %>%
  summarize(slope = coef(lm(BW ~ rel_days, data = cur_data()))[2])
# View the results
view(slopes)

write_xlsx(slopes, "slopes.xlsx")  # Save the dataframe as an Excel file

#FI OVER TIME FOR NZO AND C57 ANIMALS####
FI <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(2,3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study
  group_by(ID) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTED ANIMALS FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  drop_na(corrected_intake_kcal) %>% 
  mutate(rel_days= DATE - min(DATE)) %>% 
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal))

##Graph####
FI %>% 
  ggplot(aes(x = rel_days, y = cumsum_kcal)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  labs(
    x = "Days",
    y = "Cumulative Food Intake (kcal)"
  ) +  # White background
  theme_classic()  +
  facet_grid(SEX~STRAIN)
 # geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) 
##Main model ####
slopes_fi <- FI %>%
  group_by(ID) %>%
  summarize(slope = coef(lm(cumsum_kcal ~ rel_days, data = cur_data()))[2])
# View the results
view(slopes_fi)

write_xlsx(slopes_fi, "slopes_fi.xlsx")  # Save the dataframe as an Excel file

## Luis model check BW####
## this uses normalized data (z-scores)
mdl_data <- bw %>% 
    group_by(ID) %>% 
    mutate(
        ID = as.factor(ID),
        time = as.numeric(DATE - min(DATE))
    ) %>% 
    ungroup() %>% 
    mutate(
        time = scale(time) # this is for model convergence
    ) %>% 
    select(
        ID, BW, time, STRAIN, SEX, DIET_CODE
    ) %>% 
    ungroup()


# here I build the statistical model with time as random slope
# and control for relevant variables such as strain, sex, diet 
lmer_bw <- lmerTest::lmer(
    data = mdl_data,
    BW ~ time + STRAIN + SEX + DIET_CODE + (1+time|ID)
)
summary(lmer_bw)
view(coef(lmer_bw))

# here I extract the individual coefficients assigned to time, which is the same as the "slope"
# these are the random slopes
lmer_slopes <- coef(lmer_bw)$ID %>% 
    rownames_to_column("ID") %>% 
    select(ID, time) %>% 
    left_join(mdl_data, ., by = "ID",
              suffix = c("_zscore", "_slope")) %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    slice(1) %>% 
    select(ID, STRAIN, DIET_CODE, SEX, time_slope)
lmer_slopes

write_xlsx(lmer_slopes, "lmer_slopes.xlsx")  # Save the dataframe as an Excel file

# distribution of slopes
# nzo are brutal, look how the even surpass HFD c57 males!
lmer_slopes %>% 
    ggplot(aes(
        interaction(STRAIN, SEX, DIET_CODE), time_slope
    )) +
    geom_boxplot() +
    geom_point()

## Luis model check FI####
## this uses normalized data (z-scores)
mdl_data_fi <- read_csv("../data/FI.csv") %>% 
  group_by(ID) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTED ANIMALS FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  drop_na(corrected_intake_kcal) %>% 
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal)) %>% 
  ungroup() %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  select(
    ID, cumsum_kcal, time, STRAIN, SEX, DIET_CODE
  ) %>% 
  ungroup()

# here I build the statistical model with time as random slope
# and control for relevant variables such as strain, sex, diet 
lmer_fi <- lmerTest::lmer(
  data = mdl_data_fi,
  cumsum_kcal ~ time + STRAIN + SEX + DIET_CODE + (1+time|ID)
)
summary(lmer_fi)
view(coef(lmer_fi))

# here I extract the individual coefficients assigned to time, which is the same as the "slope"
# these are the random slopes
lmer_slopes_fi <- coef(lmer_fi)$ID %>% 
  rownames_to_column("ID") %>% 
  select(ID, time) %>% 
  left_join(mdl_data, ., by = "ID",
            suffix = c("_zscore", "_slope")) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  slice(1) %>% 
  select(ID, STRAIN, DIET_CODE, SEX, time_slope)
lmer_slopes_fi

write_xlsx(lmer_slopes_fi, "lmer_slopes_fi.xlsx")  # Save the dataframe as an Excel file

# distribution of slopes
# nzo are brutal, look how the even surpass HFD c57 males!
lmer_slopes_fi %>% 
  ggplot(aes(
    interaction(STRAIN, SEX, DIET_CODE), time_slope
  )) +
  geom_boxplot() +
  geom_point()





