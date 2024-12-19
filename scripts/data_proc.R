pacman::p_load(
  tidyverse,
  googledrive
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

# read sable data ----

sable_csv_files <- list.files(
    path = "~/Downloads/sable/sable_csv",
    full.names = TRUE)
sable_csv_files

metadata <- read_csv("../data/META.csv")

x <- read_csv(sable_csv_files[[1]])

# read





# we also need the META file
META <- read_csv("data/META.csv") %>% 
    pivot_longer(cols = starts_with("SABLE_DAY_"),
                 names_to = "sable_idx",
                 values_to = "metadata_code")
# then we create a column with the sable code m/d/y-cage
sable_data_tmp <- sable_data %>% 
    mutate(
        date = lubridate::as_date(lubridate::mdy_hms(DateTime))
    ) %>% 
    group_by(date) %>% 
    group_split()
library(furrr)
plan(multisession, workers = 8)
match_id <- sable_data_tmp %>% 
    future_map(., function(X){
        file_name <- paste0("data/sable/",X$date[1], ".rds")
        d <- X %>%
            pivot_longer(cols = matches("*_[0-9]+")) %>% 
            mutate(
                cage_number = str_extract(str_extract(name, "_[0-9]+"), "[0-9]+"),
                metadata_code = paste(format(date, "%m/%d/%Y"),
                                      cage_number, sep = "-")
            )
        tmp <- left_join(d, META, by = "metadata_code")
        saveRDS(tmp, file = file_name, compress = TRUE)
        d <- NULL
        tmp <- NULL
    })
    
