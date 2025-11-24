# TASK 4: To create the median split setup based on these z-scores

# libraries ----
library(dplyr)
library(readr)
library(ggplot2)
setwd(this.path::here())

scaled_empirical <- read_csv("../data/scaled_empirical_task3.csv") # All variables are on the same scale (mean = 0, SD = 1)

# median split within groups ----
median_split_df <- scaled_empirical %>% # the goal is to turn a continuous variable into a categorical one for classification
    mutate(
        y1_group_median = median(y1_z, na.rm = TRUE),
        y2_group_median = median(y2_z, na.rm = TRUE),
        y3_group_median = median(y3_z, na.rm = TRUE),
        y4_group_median = median(y4_z, na.rm = TRUE),
        y1_split = if_else(y1_z >= y1_group_median, "high_bw_delta", "low_bw_delta"), # each variable is transformed into a binary categorical variable (y1_split), labeled “High” or “Low”. Do we want that? ask LL
        y2_split = if_else(y2_z >= y2_group_median, "high_adi_delta", "low_adi_delta"),
        y3_split = if_else(y3_z >= y3_group_median, "high_bw_roc", "low_bw_roc"),
        y4_split = if_else(y4_z >= y4_group_median, "high_adi_roc", "low_adi_roc")
    ) %>%
    ungroup() %>%
    mutate(
        median_split_combo = paste(y1_split, y2_split, y3_split, y4_split, sep = "-")
    )

length(unique(median_split_df$median_split_combo))

table(median_split_df$median_split_combo)

table(median_split_df$y1_split)
table(median_split_df$y2_split)
table(median_split_df$y3_split)
table(median_split_df$y4_split)

## plot ----

ggplot(median_split_df, aes(x = y1_z, fill = y1_split)) + # y1_z is an example
    geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
    facet_wrap(~ SEX + STRAIN + DIET_FORMULA) +
    theme_minimal()

write_csv(median_split_df, "../data/median_split_task4.csv")
