
#Libraries####
library(dplyr) #to open a RDS and use pipe
library(tidyr) #to use cumsum
library(ggplot2)
library(readr)
library(lmerTest)
library(emmeans)

#Change in bw after ~11 weeks with LFD in NZO mice ####

#bw data import
mdl_data_bw <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study 
  #filter(DIET_FORMULA == "D12450Ki") %>% #We exclude chow data I dont know why still appear the three 3 dates that correspond to chow
  filter(DATE <="2025-02-24") %>% #we want data prior to diet restriction
  group_by(ID) %>% 
  mutate(BW_rel=((BW - min(BW))/min(BW)*100)) %>% 
  #mutate(DIET_STATUS = if_else(ID %in% c(3708,3714,3720,3721,3710,3722,3723,3724,3725,3727,3728,3729), "restricted", "ad_lib")) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  select(
    ID, BW, time, DIET_FORMULA, BW_rel, DATE
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols="BW") %>% 
  select(
    ID, time, name, value
  ) 

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_bw <- lmerTest::lmer(
  data = mdl_data_bw,
  value ~ time + (1+time|ID)
)
summary(lmer_bw)
view(coef(lmer_bw))
coef_bw <- coef(lmer_bw)$ID %>% 
rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_bw= `(Intercept)`)

write_csv(coef_bw, "coef_bw.csv")  # Save the dataframe as an csv file



#plot 1 bw over days ####
mdl_data_bw %>% 
  ggplot(aes(x = time, y = BW_rel)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  labs(
    x = "Day",
    y = "BW (% of change)"
  ) +  # White background
  theme_classic()  
#does this means that the x axis means the different days in which the IDs started with the LFD

#plot 2 bw over days but with echoMRI data ####

mdl_data_bw  %>% 
  ggplot(aes(x = DATE, y = BW)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  geom_vline(xintercept = as.numeric(echomri_dates), 
             linetype = "dashed", 
            color = "maroon") + #added echoMRI days
  labs(
    x = "Day",
    y = "BW (g)"
  ) +  # White background
  theme_classic()  

#FI for NZO####

mdl_data_fi <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study 
 # filter(DIET == "LFD") %>%  #We exclude chow data I dont know why still appear the three 3 dates that correspond to chow
  group_by(ID) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTION PERIOD FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  drop_na(corrected_intake_kcal) %>% 
  mutate(cumsum_kcal= cumsum(corrected_intake_kcal)) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  select(
    ID, cumsum_kcal, time, DIET
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols="cumsum_kcal") %>% 
  select(
    ID, time, name, value
  ) 

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_fi <- lmerTest::lmer(
  data = mdl_data_fi,
  value ~ time + (1+time|ID)
)
summary(lmer_fi)
view(coef(lmer_fi))
coef_fi <- coef(lmer_fi)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_fi= `(Intercept)`)

write_csv(coef_fi, "coef_fi.csv")  # Save the dataframe as an csv file

# echoMRI####

#Data import of echoMRI
echomri_data <- read_csv("../data/echomri.csv") %>% 
  filter(COHORT %in% c(3, 4, 5)) %>%
  filter(ID != 3715) %>% #Died during the study
  select(Date,ID, adiposity_index, Fat, Lean) %>% 
  mutate(Date = as.Date(Date, format="%m/%d/%Y")) %>% 
  rename(DATE = Date) %>% 
  filter(DATE <="2025-02-24") %>% #we want data prior to diet restriction
  group_by(ID) %>% 
  mutate(
    ID = as.factor(ID),
    time = as.numeric(DATE - min(DATE))
  ) %>% 
  mutate(
    time = scale(time) # this is for model convergence
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols=c ("adiposity_index", "Fat", "Lean")) %>% 
  select(
    ID, time, name, value
  ) 
# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_adi <- lmerTest::lmer(
  data = echomri_data %>% filter(name=="adiposity_index"),
  value ~ time + (1+time|ID)
)
summary(lmer_adi)
view(coef(lmer_adi))
coef_adi<- coef(lmer_adi)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_ai= `(Intercept)`)

write_csv(coef_adi, "coef_adi.csv")  # Save the dataframe as an csv file

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_fat <- lmerTest::lmer(
  data = echomri_data %>% filter(name=="Fat"),
  value ~ time + (1+time|ID)
)
summary(lmer_fat)
view(coef(lmer_fat))
coef_fat<- coef(lmer_fat)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_fat= `(Intercept)`)

write_csv(coef_fat, "coef_fat.csv")  # Save the dataframe as an csv file

# here I build the statistical model with time as random slope
# and control for time and ID (not all animals start the same day with LFD)
lmer_lean <- lmerTest::lmer(
  data = echomri_data %>% filter(name=="Lean"),
  value ~ time + (1+time|ID)
)
summary(lmer_lean)
view(coef(lmer_lean))
coef_lean<- coef(lmer_lean)$ID %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(ID=rowname,starting_lean= `(Intercept)`)

write_csv(coef_lean, "coef_lean.csv")  # Save the dataframe as an csv file

