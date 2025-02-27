#########################
#R script for food restriction in C57BL6J mice cohorts 2
## 2-27-25
#Average food consumption calculations for food restriction

pacman::p_load(
  tidyverse,
  ggplot2,
  ggpubr,
  furrr,
  zoo,
  TTR
)

# helper functions
percent_change <- function(restricted_intake, adlib_intake){
  return(
    ((restricted_intake - adlib_intake) / adlib_intake) * 100
  )
}


#Add Group_assignment_AI.csv to FI.csv so that you have FI and diet group for each mouse 
#Calculate rolling average of daily food intake using all mice (I commented out 
#the step that selected just ad libitum mice)
#FI intake over days####
FI_data_ <- read_csv("../data/FI.csv") %>% 
  filter(COHORT ==2)


df_restrict <-FI_data_ %>%
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, DATE, corrected_intake_gr) %>%
  drop_na(corrected_intake_gr) %>% 
 # filter(ID %in% c(7866,7863, 7872, 7874, 7861, 7865, 7877, 7878)) %>% # ONLY RESTRICTED ANIMALS
  mutate(moving_avg = rollmean(corrected_intake_gr, k=3, fill=NA, align='right')) %>%
  filter(DATE > '2025-01-21') 
  

# confirm the rollmean using spaghetti plots (verify no errors in data entry)
#Need multiple aes because you want multiple lines. 
#Black line shows raw corrected intake values. Red line shows moving average
#facet_wrap gives you multiple graphs (one for each ID rather than showing one compiled graph)
sp1 <- df_restrict %>% 
  mutate(date = lubridate::ymd(DATE)) %>% 
  ggplot(aes(
    date, corrected_intake_gr
  )) +
  geom_line(aes(group=ID)) +
  geom_line(aes(date, moving_avg, group=ID), color="red") +
  facet_wrap(~ID, scale="free_y")
sp1

# take individuaL adlib intake measurements and see their respective percent change
# with the overall 60% restriction
last_intake <- df_restrict %>%
  group_by(ID) %>% 
  mutate(DATE = lubridate::ymd(DATE)) %>% 
  slice_max(order_by = DATE, n = 1) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(
    restricted_intake = mean(moving_avg)*0.6,
    relative_restriction_gr = (restricted_intake - moving_avg),
    relative_restriction_perc = percent_change(restricted_intake, moving_avg)
  )
#Note that this uses the percent_change function that we made at the start of the script
#For percent change function the order of the two parameters (restricted_intake and moving_avg) matters
#because it is plugging the parameters into the places in the function

# data visualization
#graph
last_intake %>% 
  filter(Restriction_group=="Restrict") %>% 
  ggplot(aes(as.factor(ID), relative_restriction_gr)) +
  geom_col() +
  geom_label(aes(label = round(relative_restriction_perc,1))) +
  scale_y_continuous(breaks = seq(-5, 0, 0.5))

#same graph as above, but with labeled bars  
last_intake %>% 
  filter(Restriction_group=="Restrict") %>% 
  ggplot(aes(as.factor(ID), relative_restriction_gr)) +
  geom_col() +
  geom_label(aes(label = paste(round(moving_avg,1), round(restricted_intake,1),
                               round(relative_restriction_perc,1),sep = " -> ")) ) +
  scale_y_continuous(breaks = seq(-5, 0, 0.5))

write_csv(x = last_intake, "/Users/laurenmichels/Documents/GitHub/data/data/Food_restriction.csv")

#Just for practical use to feed mice 
#This section tells me how much to feed per meal for each mouse in the restricted group
Per_meal <-last_intake %>%
  group_by(ID) %>% 
  select(ID,DATE,moving_avg,Restriction_group,restricted_intake,relative_restriction_perc)%>%
  filter(Restriction_group=="Restrict") %>%
  mutate (LFD_per_meal = restricted_intake*0.5)

#_____________________________________________________
#This was used originally when I thought we would be giving all of the mice the same amount of food
#Now that we decided to restrict them on an individual basis it is no longer relevant
#For each mouse, calculate avg daily food consumption based on the most recent 3 measurements
#Calculate the avg of these values across all 23 mice. Calculate 60% of this.
# taking the last rolling mean measurement of each animal 
#Slice_tail with n=1 selects the most recent rolling average value for each mouse

#df_restrict2 <- df_restrict %>% 
#ungroup() %>% 
#group_by(ID) %>% 
#slice_tail(., n=1) %>%
#ungroup() %>%
#summarise(
#last_measurement_mean = mean(moving_avg)
#) %>% 
# mutate(
# restricted_daily_food_gr = last_measurement_mean*0.6,
#restricted_per_meal_gr = restricted_daily_food_gr*0.5
#) %>%
#ungroup()
#df_restrict2

#Running this file will show the spaghetti plot and in the console it will give the mass 
#of food to feed the restricted mice for each meal
  