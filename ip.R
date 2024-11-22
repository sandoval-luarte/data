pacman::p_load(
    tidyverse,
    ggplot2
)

raw_data <- readxl::read_xlsx("~/tmp/SourceData_Fig2.xlsx")


raw_data %>% 
    ggplot(aes(x = factor(State, level = c("Wheel", "HPFWheel")), MetersMaze)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_line(aes(group = Mouse), alpha = 0.5) + 
    facet_grid(~factor(Inj, levels=c("Veh", "Alm"))) +
    scale_y_continuous(limits = c(0, 50), expand = c(0,0))

raw_data %>% 
    ggplot(aes(Inj, MetersMaze)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    facet_wrap(~State)

mdl1 <- lmerTest::lmer(
    data = raw_data,
    MetersMaze ~ Inj + State + Sex + (1|Mouse)
)
summary(mdl1)

mdl2 <- lmerTest::lmer(
    data = raw_data %>% filter(State == "HPFWheel"),
    MetersMaze ~ Inj + Sex + (1|Mouse)
)
summary(mdl2)