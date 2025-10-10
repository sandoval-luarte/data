# AIM task 3> take all y₁…y₄ empirical measures (and possibly the standardized residuals) and generate z-scores within each sex×strain×diet group.

# libraries####
library(dplyr)
library(readr)

final_df <- read_csv("../data/final_df_task2.csv")

scaled_empirical <- final_df %>%
  group_by(SEX, STRAIN, DIET_FORMULA) %>%
  mutate(
    y1_z = scale(y1),
    y2_z = scale(y2),
    y3_z = scale(y3),
    y4_z = scale(y4)
  ) %>%
  ungroup()

scaled_empirical <- scaled_empirical %>%
  mutate(across(ends_with("_z"), as.numeric)) #again we have some variables that are matrix


write_csv(scaled_empirical, "../data/scaled_empirical_task3.csv")

