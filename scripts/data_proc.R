pacman::p_load(
  tidyverse
)

fname <- rstudioapi::selectDirectory()

# bodyweight and food intake ----

filter_files <- tibble(
  filepath = list.files(fname, full.names = TRUE)
) %>% 
  filter(
    grepl("COHORT_[0-9]+.csv", filepath)
  )
filter_files

open_files <- filter_files %>% 
  mutate(r = row_number()) %>% 
  group_by(r) %>% 
  group_split() %>% 
  map_dfr(
    ., function(X){
      read_csv(X$filepath) %>% 
        select(ID, FOOD_WEIGHT_START_G, FOOD_WEIGHT_END_G, DATE, DIET, BODY_WEIGHT_G, DIET_FORMULA) %>% 
        mutate(
          INTAKE_GR = (FOOD_WEIGHT_START_G - FOOD_WEIGHT_END_G),
          DATE = lubridate::mdy(DATE)
        ) %>% 
        select(ID, INTAKE_GR, DATE, BODY_WEIGHT_G, DIET_FORMULA) %>% 
        rename(
          BW = BODY_WEIGHT_G
        )
    }
  )

# load food description

food_desc <- read_csv(paste(fname, "/food_description.csv", sep = ""))

# load metadata

metadata <- read_csv(paste(fname, "/META.csv", sep = "")) %>% 
    select(ID, SEX, COHORT, STRAIN, AIM, DIET_CODE)

# output food-intake file

FI <- open_files %>% 
  select(ID, DIET_FORMULA, INTAKE_GR, DATE) %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(
    delta_measurement = DATE - lag(DATE)
  ) %>% 
  drop_na(delta_measurement) %>% 
  mutate(
    corrected_intake_gr = INTAKE_GR / as.numeric(delta_measurement)
  ) %>% 
    left_join(., food_desc, by = "DIET_FORMULA") %>% 
    mutate(
        corrected_intake_kcal = corrected_intake_gr * KCAL_G
    ) %>% 
    left_join(., metadata, by = "ID")

# output bodyweight file
BW <- open_files %>% 
  group_by(ID) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  select(ID, BW, DATE) %>% 
  drop_na(BW) %>% 
    left_join(., metadata, by = "ID")

write_csv(x = FI, "../data/FI.csv")
write_csv(x = BW, "../data/BW.csv")

# echo MRI ----

fname2 <- rstudioapi::selectDirectory()

filter_files2 <- tibble(
    filepath = list.files(fname2, full.names = TRUE)
) %>% 
    filter(
        grepl("*.xlsx", filepath)
    )
filter_files2

open_files2 <- filter_files2 %>% 
    mutate(r = row_number()) %>% 
    group_by(r) %>% 
    group_split() %>% 
    map(., function(X){
        readxl::read_xlsx(X$filepath) %>% 
            select(Label, Fat, Lean, Weight, TimeDateDura) %>% 
            rename(ID = Label) %>% 
            separate_wider_delim(TimeDateDura, delim = ";", names = c("Date", "A", "B")) %>% 
            select(-A, -B) %>% 
            separate_wider_delim(Date, delim = " ", names = c("hms", "month", "day", "year")) %>% 
            mutate(day = gsub(",", "", day),
                   Date = paste(year, month, day, sep = "-"),
                   Date = lubridate::ymd(Date),
                   ID =  as.numeric(ID)) %>% 
            select(-hms, -month, -day, -year)
    }) %>% 
    bind_rows() %>% 
    left_join(., metadata, by = "ID")

open_files2

# compare adiposity index = fat / lean ----

echomri_data <- open_files2 %>% 
    mutate(adiposity_index = Fat / Lean) %>% 
    group_by(ID) %>% 
    mutate(
        n_measurement = as.numeric(as.factor(Date))
    )
echomri_data

write_csv(x = echomri_data, "../data/echomri.csv")