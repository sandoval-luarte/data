#aim####
#This script aim to explore fisetin data sent by Morgan ward to CS on march 21-2025
#libraries####
library(readxl)
library(readr)
library(ggplot2)
library(dplyr)
library(emmeans)
library(lme4)
library(lmerTest) 

#Data import
rm(list = ls())
fisetin <- read_excel("Documents/GitHub/data/scripts/Fisetin_project/Fisetin qPCR data All.xlsx")
View(fisetin)
#Exploration
sapply(fisetin, class) #which class is each column
colnames(fisetin) #what we have in the kitchen
#Here I believe Age, deltaCT and Fold change should be class numeric
colnames(fisetin) <- make.names(colnames(fisetin))
fisetin <- fisetin %>% 
  rename(
    ID = "Mouse.ID",
    delta_ct = "X..Ct",
    fold_change = "Fold.Change",
    SEX = "Sex"
  ) %>% 
  mutate(delta_ct= as.numeric(delta_ct), fold_change=as.numeric(fold_change)) %>% 
  filter(Cohort==5) #Just to check animals that I know
#Left joint to cohort 1 to add complete the lack of info in SEX column
COHORT_1 <- read_csv("Documents/GitHub/data/data/COHORT_1.csv") %>% 
  distinct(ID, .keep_all = TRUE) %>%  # Keep only the first row for each unique ID
  select(ID, SEX) #Extract only sex to do a left joint with fisetin data
merged_df <- left_join(fisetin, COHORT_1, by = "ID", suffix = c("", "_unique")) # Merge fisetin dataframe with the COHORT_1 based on 'ID'
merged_df <- merged_df %>% # Update the 'SEX' column in fisetin with the unique 'SEX' values of COHORT_1
  mutate(SEX = ifelse(is.na(SEX), SEX_unique, SEX))
# Remove the redundant 'SEX_unique' column
merged_df <- merged_df %>%
  filter(!(is.na(delta_ct) & is.na(fold_change))) %>% 
  select(-SEX_unique)

#exploratory plot####
merged_df %>%
  group_by(SEX, Age, Diet, Treatment, Gene) %>%
  ggplot(aes(x = Gene, y = fold_change)) +
  geom_point() +
  facet_wrap(~Treatment)
# 3-way ANOVA####
# Filter for the specific gene
gene_data <- merged_df %>%
  filter(Gene == "p21")  # gene of interest
# Fit the ANOVA model
anova_model <- aov(fold_change ~ Diet * SEX * Treatment, data = merged_df)
# View the ANOVA table
summary(anova_model)

#cohort 1 EchoMRIdata ####
echomri <- read_csv("Documents/GitHub/data/data/echomri.csv") %>% 
  filter(COHORT ==1)
#exploratory plot for Adiposity index ####
echomri %>%
  group_by(SEX, DIET_FORMULA, ID) %>%  # Ensure ID is included for grouping
  ggplot(aes(x = as.factor(n_measurement), y = adiposity_index, group = as.factor(ID), color = as.factor(ID))) + 
  geom_point() + 
  geom_line() +  # Add lines to connect points by ID
  facet_wrap(SEX ~ DIET_FORMULA)
#Here these mice should have 3 measurement days and not two as show the graph
#We must retrieve the echoMRI data from 4/1/24 
#CS will check on 03/24/25 if these data is within the backup disk

#LMER
model <- lmer(adiposity_index ~ DIET_FORMULA * SEX + (1 | ID), data = echomri)
summary(model)
emmeans(model, pairwise ~ DIET_FORMULA * SEX)
#CS also will go deep in the raw data in the morning 
#CS will do a diagram of the experiment (what were done)