library(dplyr)
library(readr)
library(ggplot2)

# --- Import CSV files ---
pedmeters_data <- read_csv("../data/sable_pedmeters_data.csv")
total_act_data <- read_csv("../data/sable_locomotion_data.csv")

# --- Left join by ID and SABLE, then calculate SPA ---
spa_data <- total_act_data %>%
  left_join(
    pedmeters_data %>% select(ID, SABLE, total_pedmeters),
    by = c("ID", "SABLE")
  ) %>%
  mutate(SPA = total_act - total_pedmeters)

# --- Inspect result ---
head(spa_data)


# --- Collapse into PlotGroup by STRAIN ---
spa_plotdata <- spa_data %>%
  group_by(STRAIN) %>%  # ensure grouping by strain
  mutate(
    PlotGroup = case_when(
      SABLE %in% c("baseline", "peak obesity") ~ "all",
      SABLE %in% c("BW loss", "BW maintenance") ~ GROUP,
      SABLE == "BW regain" ~ paste(GROUP, DRUG, sep = "_")
    ),
    # reorder SABLE as factor
    SABLE = factor(SABLE, levels = c("baseline", "peak obesity", 
                                     "BW loss", "BW maintenance", "BW regain"))
  ) %>%
  ungroup()  # remove grouping for plotting



# --- Custom color palette ---
custom_colors <- c(
  "all" = "gray70",
  "ad lib" = "#E67E22",
  "restricted" = "#3498DB",
  "ad lib_vehicle" = "#E67E22",
  "ad lib_RTIOXA_47" = "#F39C12",
  "restricted_vehicle" = "#3498DB",
  "restricted_RTIOXA_47" = "#5DADE2"
)

# --- Plot SPA ---
spa_plot <- spa_plotdata %>%
  ggplot(aes(x = SABLE, y = SPA, fill = PlotGroup)) +
  
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7, alpha = 0.7) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.3) +
  
  geom_point(aes(color = PlotGroup),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  theme_minimal() +
  labs(y = "Spontaneous Physical Activity (m/day)", fill = "Group", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~STRAIN)

# --- Display plot ---
spa_plot

