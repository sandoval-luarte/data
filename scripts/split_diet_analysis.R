
echomri_data_2 <- echomri_data %>% 
  filter(COHORT ==7) %>% 
  select(ID, adiposity_index,SEX) %>% 
  mutate(DIET = if_else(ID %in% c(2322, 2321,2319,2330,2328,2320,2316,2313,2312,2317,2314,2307), "HFD", "LFD")) %>% 
  ggplot(aes(DIET, adiposity_index), color = SEX) +
  geom_point() +
  geom_line() +
  facet_wrap(~SEX)

#check if the adiposity index between HFD and LFD are different using t-test
t_test_result <- lm(adiposity_index ~ DIET, data = echomri_data_2)
# View the result
print(t_test_result)
# Summary
summary(t_test_result)
# read sable data ----