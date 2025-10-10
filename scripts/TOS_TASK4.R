#TASK 4: To create the median split setup based on these z-scores

# libraries ----
library(dplyr)
library(readr)
library(ggplot2)

scaled_empirical <- read_csv("../data/scaled_empirical_task3.csv"). #All variables are on the same scale (mean = 0, SD = 1)

# median split within groups ----
median_split_df <- scaled_empirical %>% #the goal is to turn a continuous variable into a categorical one for classification
  group_by(SEX, STRAIN, DIET_FORMULA) %>% #median is calculated within each combination of SEX, STRAIN, and DIET_FORMULA
  mutate(
    y1_group_median = median(y1_z, na.rm = TRUE),
    y2_group_median = median(y2_z, na.rm = TRUE),
    y3_group_median = median(y3_z, na.rm = TRUE),
    y4_group_median = median(y4_z, na.rm = TRUE),
    y1_split = if_else(y1_z >= y1_group_median, "High", "Low"), #each variable is transformed into a binary categorical variable (y1_split), labeled “High” or “Low”. Do we want that? ask LL
    y2_split = if_else(y2_z >= y2_group_median, "High", "Low"),
    y3_split = if_else(y3_z >= y3_group_median, "High", "Low"),
    y4_split = if_else(y4_z >= y4_group_median, "High", "Low")
  ) %>%
  ungroup()

table(median_split_df$y1_split, median_split_df$SEX)
summary(median_split_df$y1_z)

##plot ----

ggplot(median_split_df, aes(x = y1_z, fill = y1_split)) + #y1_z is an example
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  facet_wrap(~ SEX + STRAIN + DIET_FORMULA) +
  theme_minimal()

write_csv(median_split_df, "../data/median_split_task4.csv")

