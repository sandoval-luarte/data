#description####
#This script aims to explore orexin cre basic feeding behavior (LL animals)
#cohort 10
#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)

#bw data import####

weights_copy <- read_csv("weights copy.csv") %>%
  mutate(date_ymd = as.Date(date_ymd, format = "%m/%d/%y")) %>%
  select(date_ymd, id, treatment, weight_gr, pellet_intake, notes) %>%
  drop_na(weight_gr) %>%
  group_by(id) %>%
#  mutate(ref_weight = first(weight_gr[date_ymd == as.Date("2025-02-19")], default = NA),  # Obtener peso de referencia
  #       weight_relative = (weight_gr - min(weight_gr))/min(weight_gr)*100) %>% 
  mutate(comment = if_else(is.na(notes) | notes == "", FALSE, TRUE))

plot <- weights_copy %>% 
  ggplot(aes(x=date_ymd,y=weight_gr, color = as.factor(id), group = id)) +
  geom_point(size = 3) +             # Plot points
  geom_line()   
plot

#FI#####
FI_copy <- read_csv("weights copy.csv") %>%
  mutate(date_ymd = as.Date(date_ymd, format = "%m/%d/%y")) %>%
  select(date_ymd, id, treatment, weight_gr, pellet_intake, notes) %>%
  drop_na(pellet_intake) %>%
  group_by(id) %>%
  #  mutate(ref_weight = first(weight_gr[date_ymd == as.Date("2025-02-19")], default = NA),  # Obtener peso de referencia
  #       weight_relative = (weight_gr - min(weight_gr))/min(weight_gr)*100) %>% 
  mutate(comment = if_else(is.na(notes) | notes == "", FALSE, TRUE))

plot <- FI_copy %>% 
  ggplot(aes(x=date_ymd,y=pellet_intake, color = as.factor(id), group = id)) +
  geom_point(size = 3) +             # Plot points
  geom_line()   
plot
