# This script aims to explore Megan Sumera Liu Lab
#libraries----
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr) 
library(purrr)
library(Hmisc)
library(lme4)
library(lmerTest)
library(emmeans)
library(pracma)
library(lubridate)
library(stringr)
library(forcats)
library(patchwork)
library(ggpattern)
library(car)
library(broom) 
library(rstatix)

#exploration part1----

data1 <- read_csv("../data/OXYE090MCK.csv") %>% #I believe this data is more fun because contain data from the other files. 
  drop_na()
data1 <- data1 <- data1[-c(1), ]
data1 <- data1 <- data1[-c(2), ]
names(data1) <- as.character(unlist(data1[1, ]))
data1 <- data1 <- data1[-c(1), ] #ok now data1 is clean

#there is many thinks I don't know, for example "Interval". 
# I imagine "Chan" is "channel". Because is always the same number = 0107
# Consequently, I believe this csv file correspond to an unique ID (just 1 animal)
# So the first thing Morgan should help me is to define what the columns names are

#Just for fun I will explore locomotory activity. I imagine correspond to columns 
# "XAMB", "YAMB" and "ZAMB". What is the difference between "XTOT" and "XAMB"?
#I can Imagine the sum of activity in X, Y and Z could be "total locomotion?" let's see

data1_locomotion <- data1 %>% 
select(CHAN, `DATE/TIME`, XTOT, YTOT, ZTOT,XAMB, YAMB)

#After looking this data it seems this experiment started 9:08 am in the morning of 7/6/26 
#and ended 7/10/26 8:26 am so I imagine this data included acclimation period in Oxymax ical

#So long story short this animals stayed 4 complete days in the Oxymax machine
#I noticed , these values represent non cumulative counts. I will focus just in XAMB

#How is the light off/on period in the ical room? 
#let's assume is from 6 am to 6 pm the light on period

data1_locomotion <- data1 %>% 
  select(CHAN, `DATE/TIME`, XAMB)

#We can expect at least the first day this ID will move more in comparison with the following days

data1_locomotion <- data1_locomotion %>%
  mutate(datetime = mdy_hm(`DATE/TIME`), # Convert DATE/TIME
    XAMB = as.numeric(XAMB), # is XAMB numeric? I hope so if not we transformed it
lights = case_when(     # Define lights  
      hour(datetime) < 6 ~ "lights_off",
      hour(datetime) == 6 & minute(datetime) == 0 ~ "lights_off",
      hour(datetime) < 18 ~ "lights_on",
      TRUE ~ "lights_off"
    ),
    experimental_date = case_when(# Define experimental day   
      hour(datetime) < 6 ~ as.Date(datetime) - 1,
      hour(datetime) == 6 & minute(datetime) == 0 ~ as.Date(datetime) - 1,
      TRUE ~ as.Date(datetime)
    )) %>% 
  select(CHAN,XAMB, datetime, lights,experimental_date) %>%
  mutate( # Convert experimental dates to Day 1, Day 2, etc. #day 1 started 6 pm and end 6 am of the following day
    day = as.integer(factor(experimental_date))
  ) %>%
  select(CHAN,XAMB, datetime, lights,day)%>%
  group_by(CHAN, day) %>%   # Calculate cumulative XAMB
  arrange(datetime, .by_group = TRUE) %>%
  mutate(
    cumulative_XAMB = cumsum(XAMB)
  ) %>%
  ungroup()

data1_locomotion <- data1_locomotion %>%
  mutate(
    XAMB = as.numeric(XAMB),
    
    lights = case_when(
      hour(datetime) >= 6 & hour(datetime) < 18 ~ "lights_on",
      TRUE ~ "lights_off"
    ),
    
    experimental_date = as.Date(datetime - hours(6))
  ) %>%
  mutate(
    day = as.integer(factor(experimental_date))
  ) %>%
  arrange(CHAN, datetime) %>%
  
  # cumulative XAMB for the entire experimental day
  group_by(CHAN, day) %>%
  mutate(
    cumulative_XAMB = cumsum(XAMB)
  ) %>%
  ungroup()

data1_locomotion_day_summary <- data1_locomotion %>%
  group_by(CHAN, day) %>%
  summarise(
    cumulative_XAMB = max(cumulative_XAMB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    CHAN == "0107",
    day != 5
  )

#plot 1----

plot1 <- ggplot(
  data1_locomotion_day_summary,
  aes(x = factor(day), y = cumulative_XAMB)
) +
  geom_col() +
  labs(
    x = "Day",
    y = "Cumulative XAMB"
  ) +
  theme_classic()

plot1

## stats ----

# We have only one animal (CHAN = "0107"), so we cannot perform
# a repeated-measures ANOVA or mixed-effects model using animal as
# the repeated/random factor.

# Instead, we can explore whether cumulative XAMB shows a linear
# trend across experimental days for this individual animal.

# Day is treated as a continuous variable here.
# Therefore, the model tests whether cumulative XAMB tends to
# increase or decrease linearly as the number of days in Oxymax increases.

model_trend <- lm(
  cumulative_XAMB ~ day,
  data = data1_locomotion_day_summary
)

summary(model_trend)

#So for this animal, cumulative XAMB decreased by approximately 1,357 counts for each additional day in the Oxymax chamber.
#conclusion: adaptation to the chamber matters

# analysis separated by light----
data1_locomotion_lights <- data1 %>% 
  select(CHAN, `DATE/TIME`, XAMB) %>%
  mutate(
    datetime = mdy_hm(`DATE/TIME`),
    XAMB = as.numeric(XAMB),
    lights = case_when(
      hour(datetime) >= 6 & hour(datetime) < 18 ~ "lights_on",# Define light phase
      TRUE ~ "lights_off"
    ),
 experimental_date = as.Date(datetime - hours(6))  # Experimental day starts at 6 AM
  ) %>% mutate(
    day = as.integer(factor(experimental_date))
  ) %>%
  arrange(CHAN, datetime) %>%
  group_by(CHAN, day, lights) %>%
  mutate(
    cumulative_XAMB = cumsum(XAMB),
    relative_counts = cumulative_XAMB - first(cumulative_XAMB)  # Relativize so every phase starts at 0
  ) %>%
  ungroup()

data1_locomotion_lights_summary <- data1_locomotion_lights %>%
  group_by(CHAN, day, lights) %>%
  summarise(
    total_XAMB = sum(XAMB, na.rm = TRUE),
    .groups = "drop"
  )%>% 
  filter(CHAN == "0107") %>% #what event is?
  filter(day != 5) # Day 5 is incomplete, so we exclude it

#plot 2----
ggplot(
  data1_locomotion_lights_summary,
  aes(
    x = factor(day),
    y = total_XAMB,
    fill = lights
  )
) +
  geom_col(position = "dodge") +
  labs(
    x = "Day",
    y = "Total XAMB counts",
    fill = "Light phase"
  ) +
  theme_classic()

#combined plot----
y_max <- max(
  max(data1_locomotion_day_summary$cumulative_XAMB, na.rm = TRUE),
  max(data1_locomotion_lights_summary$total_XAMB, na.rm = TRUE)
)
plot1 <- ggplot(
  data1_locomotion_day_summary,
  aes(x = factor(day), y = cumulative_XAMB)
) +
  geom_col() +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Day",
    y = "Total XAMB counts",
    title = "Total daily locomotor activity"
  ) +
  theme_classic()
plot2 <- ggplot(
  data1_locomotion_lights_summary,
  aes(
    x = factor(day),
    y = total_XAMB,
    fill = lights
  )
) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Day",
    y = "Total XAMB counts",
    fill = "Light phase",
    title = "Locomotor activity by light phase"
  ) +
  theme_classic()


(plot1 + plot2) +
  plot_annotation(
    tag_levels = "A"
  )
