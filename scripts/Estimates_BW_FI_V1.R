library(dplyr)
#BW OVER TIME FOR NZO AND C57 ANIMALS####
#Our aim is obtaining the slope of the BW of NZO and C57 over the weeks
##Data####
bw <- read_csv("../data/BW.csv") %>% 
  filter(COHORT %in% c(2,3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study
  group_by(ID) %>% 
  filter(DATE < "2025-02-24") %>% #ELIMINATE RESTRICTED ANIMALS FROM THE ANALYSIS
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
  summarize(slope = coef(lm(BW ~ as.numeric(DATE), data = cur_data()))[2])
# View the results
view(slopes)
##Alternative model####
lm(BW ~ DATE + (1|ID), data = bw) #IS THIS OK (asking LL)?
#FI OVER TIME FOR NZO AND C57 ANIMALS####
FI <- read_csv("../data/FI.csv") %>% 
  filter(COHORT %in% c(2,3, 4, 5)) %>% 
  filter(ID != 3715) %>%  # Exclude animals that died during study
  group_by(ID) %>% 
  filter(DATE < "2025-02-24") %>%  #ELIMINATE RESTRICTED ANIMALS FROM THE ANALYSIS
  filter(corrected_intake_kcal < 50) %>% #Eliminate weird data that is probably typing mistake
  drop_na(corrected_intake_kcal) %>% 
  mutate(rel_days= DATE - min(DATE))

##Graph####
FI %>% 
  ggplot(aes(x = rel_days, y = corrected_intake_kcal)) + 
  geom_smooth()+
  geom_point(alpha=0.2) +
  geom_line(aes(group = ID))+
  labs(
    x = "Date",
    y = "FI (kcal)"
  ) +  # White background
  theme_classic()  +
  facet_grid(SEX~STRAIN)+
  geom_text(aes(label = ID), hjust = 0.5, vjust = -0.5, size = 3) 
##Main model ####
slopes <- FI %>%
  group_by(ID) %>%
  summarize(slope = coef(lm(corrected_intake_kcal ~ as.numeric(DATE), data = cur_data()))[2])
# View the results
view(slopes)
##Alternative model####
lm(FI ~ DATE + (1|ID), data = FI) #IS THIS OK?
