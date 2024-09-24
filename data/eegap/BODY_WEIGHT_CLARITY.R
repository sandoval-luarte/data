#README####
#This script is to:
#Generate a .csv of Body weight (BW) of male C57BL6J male mice fed with a LFD and HFD 25 weeks
#[MICE INCLUDED n=4 per group (HFD/LFD), TOTAL = 8

#Data description ####

#**MICE_ID_RAR_GROUPED*: Animals initially arrived in groups of 4-5 and shared a common RAR breeding card
#**MICE_ID_RAR_INDIVIDUAL* RAR breeding card number assigned to each mouse 
#**ID_MICE* Short version of RAR breeding card number assigned to each mouse 
#**RACK*: Rack location within 370F room in CCRB
#**LOCATION*: Location in rack
#**COHORT*: Group of animals that are studied or observed together
#**DATE_OF_BIRTH*: Date of birth of each animal
#**SEX*: Sex of included animals (Males:M or Females:F)
#**DIET*: Animals were fed with chow, HFD or LFD
#**AGE_OLD_WEEKS*: Age of the animals (in weeks)
#**WEEK_WITH_DIET*: Weeks that mice have been fed with LFD or HFD
#**FOOD_WEIGHT_START*: Amount of food left in the cage two days before
#**FOOD_WEIGHT_END*: Amount of food measure in the cage
#**DAILY_FOOD_INTAKE*: Amount of food measure in the cage per day in grams
#**DAILY_FOOD_INTAKE_KCAL*: Amount of food measure in the cage per day in kcal
#**BODY_WEIGHT_G*: or notes
#**DATE*: Date of the measurement
#**AIM*: animal destination
#**COMMENTS*: or notes
#**FAT*: Amount of fat of the animal (g)- EchoMRI measurement
#**LEAN*: Amount of lean mass of the animal (g)- EchoMRI measurement
#**FREE_WATER*: Amount of total water of the animal (g)- EchoMRI measurement
#**TOTAL_WATER*: Amount of total water of the animal (g)- EchoMRI measurement
#**NEW_FOOD_END_SESSION*: Amount of total food left to measure the amount of food eaten in the next session
#**RESEARCH_DIET_FORMULA*: Code of the diet formula for animal feeding
#**LOT_FOOD*: Lot of the food
#**COMMENT_TWO*: Notes
#**COMMENT_THREE*: Notes

#Libraries####
rm(list = ls()) # clean R environment
library(lubridate)
library(tidyverse)
library(ggplot2) #to make graphs
library(tidyverse) # to import csv
library(ggpubr)
library(rstatix)
library(lmerTest) #to run lineal models
library(emmeans)
library(dplyr)

#format plot####
format.plot <- theme_pubr() +
  theme(strip.background = element_blank(), 
    #    strip.text = element_blank(),
        plot.margin = unit(rep(0.2,4), "cm"), 
   #legend.position = "none",
        axis.text.x = element_text(family = "Helvetica", size = 14),
        axis.text.y = element_text(family = "Helvetica", size = 14),
        axis.title.y.left =element_text(family = "Helvetica", size = 16),
        axis.title.x.bottom =element_text(family = "Helvetica", size = 16))
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#Data import ####
ds <-read_csv("CLARITY.csv") %>% #ds is the dataframe with all data
  select(ID_MICE,DIET,SEX,LOCATION,DATE_OF_BIRTH,DATE,FOOD_WEIGHT_START_G,FOOD_WEIGHT_END_G,BODY_WEIGHT_G,DATE_EUTHANASIA) %>% 
  filter(DIET!= "CHOW") %>% 
    mutate(
        DATE = lubridate::mdy(DATE)
    ) %>% 
    mutate(
        DATE_EUTHANASIA = lubridate::mdy(DATE_EUTHANASIA)) %>% 
    mutate(
        DATE_OF_BIRTH = lubridate::mdy(DATE_OF_BIRTH)) %>% 
    mutate(KCAL=FOOD_WEIGHT_START_G-FOOD_WEIGHT_END_G) %>% 
    mutate(AGE= DATE_EUTHANASIA - DATE_OF_BIRTH) %>% 
    rename(ID=ID_MICE) %>% 
    rename(BW=BODY_WEIGHT_G) %>% 
    rename(LOC=LOCATION) 
    
#kcal transformation
ds$KCAL <- ifelse(ds$DIET == "HFD",
                              ds$KCAL * 4.73,  # For HFD
                              ds$KCAL * 3.82)  # For LFD

BW <-ds %>% select(ID,BW,DATE) %>% drop_na()
write.csv(BW, "BW.csv", row.names = FALSE)

FI <-ds %>% select(ID,KCAL,DATE)
write.csv(FI, "FI.csv", row.names = FALSE)

META <-ds %>% select(ID,DIET,AGE,SEX,LOC)
write.csv(META, "META.csv", row.names = FALSE)

