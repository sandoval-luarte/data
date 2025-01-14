# libs ----

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo
)

fname <- rstudioapi::selectDirectory() #route folder, the repository itself

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
    corrected_intake_gr = INTAKE_GR / as.numeric(delta_measurement) #DAILY INTAKE
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

## first we need to download the files from google drive

## step 1 auth yourself with your UMN acc

drive_auth()

## step 2 search for the correct folder, in this case "csv_files"
## this folder contains all the raw sable data in csv files

folder_id <- as_id(drive_find("csv_files"))
folder_contents <- drive_ls(folder_id)

## step 3 download the folder contents, remember to set the path to soruce file location

## select the place where you want to download the csv files
csv_dir <- rstudioapi::selectDirectory()

## this loops the folder contents and download the csv files one by one
folder_contents$id %>% 
  imap(., possibly(function(X, idx){
    path <- paste(csv_dir, "/", folder_contents$name[idx], sep = "")
    drive_download(X, path = path, overwrite = FALSE)
  }))

## get all the csv file names
sable_csv_files <- list.files(
    path = csv_dir,
    full.names = TRUE,
    pattern = "*.csv")
sable_csv_files

## read the metadata file and set it for merge with sables csv files
# set directory to source file location
metadata <- read_csv("../data/META.csv") %>% 
    pivot_longer(cols = starts_with("SABLE_DAY_"),
                 names_to = "sable_idx",
                 values_to = "metadata_code")

## merge metadata with sable data in an hourly fashion
detect_bouts <- function(X){
    X %>% 
    mutate(
      rolling_sd_raw = zoo::rollapply(value, 10, fill = 0, FUN = sd, align = "right"),
      rolling_sd = as.numeric((zoo::rollapply(value, 10, fill = 0, FUN = sd, align = "right")) > 0.02),
      bouts = as.numeric(data.table::rleid(rolling_sd))
      ) %>% 
  group_by(bouts) %>% 
  mutate(bout_duration = n())
}

compute_t_tests <- function(X){
    bout_length <- length(unique(X$bouts))
    if(bout_length > 1){
    X %>% 
    filter(rolling_sd == 0) %>% 
    group_by(bouts) %>% 
    group_split() %>% {
        a <- head(., -1)
        b <- tail(., -1)
        map2_dfr(a, b, possibly(function(X, Y){
            bouts <- Y$bouts[1]
            test_t <- t.test(X$value, Y$value, alternative = c("two.sided"))
            t_test_pval <- test_t$p.value
            diff <- mean(X$value) - mean(Y$value)
            return(
                tibble(
                    bouts = bouts,
                    t_test_pval = t_test_pval,
                    diff = diff
                )
            )
        }))
    }}
    else{
        return(
            tibble(
                bouts = NA_real_,
                t_test_pval = NA_real_,
                diff = NA_real_
            ))
    }
}


corrected_intake <- function(bouts, t_tests){
  left_join(bouts, t_tests, by = "bouts") %>%
  ungroup() %>% 
  group_by(bouts) %>%
  mutate(
    value = if_else(rolling_sd == 0, mean(value), NA_real_), # para el resumen en 1h, considerar solo los momentos donde sd = 0, tomar promedio
    corrected_diff = if_else(diff < 0 | t_test_pval > 0.05, 0, diff)
  )
}


plan(multisession, workers = availableCores())
sable_hr_data <- sable_csv_files %>% 
    future_map_dfr(
        ., function(X){
            proc_data <- X %>% 
                read_csv(., show_col_types = FALSE) %>% 
                select(c(DateTime, matches(c("kcal_hr_",
                                           "VO2_",
                                           "VCO2_",
                                           "RQ_",
                                           "AllMeters_",
                                           "PedMeters_",
                                           "PedSpeed_",
                                           "FoodA_",
                                           "Water_",
                                           "BodyMass")))) %>% 
                pivot_longer(
                    cols = matches(c("kcal_hr_",
                                     "VO2_",
                                     "VCO2_",
                                     "RQ_",
                                     "AllMeters_",
                                     "PedMeters_",
                                     "PedSpeed_",
                                     "AllMeters_",
                                     "FoodA_",
                                     "Water_",
                                     "BodyMass")),
                    names_to = "parameter",
                    values_to = "value"
                ) %>% 
                mutate(
                    DateTime = lubridate::mdy_hms(DateTime),
                    cage_number = str_extract(str_extract(parameter, "_[0-9]+"), "[0-9]+"),
                    date = lubridate::as_date(DateTime),
                    hr = lubridate::hour(DateTime),
                    metadata_code = paste(gsub('(?<=\\/)0|^0', '', format(date, "%m/%d/%Y"), perl=TRUE),
                                          cage_number, sep = "-")
                )
            food_intake <- proc_data %>% 
                filter(grepl("FoodA_", parameter)) %>% 
                group_by(cage_number, date, parameter) %>% 
                group_split() %>% 
                map_dfr(., function(X){
                    return(corrected_intake(detect_bouts(X), compute_t_tests(detect_bouts(X))))
                })
            non_food_intake <- proc_data %>% 
                filter(!grepl("FoodA_", parameter))
            complete_data <- bind_rows(food_intake, non_food_intake) %>% 
                left_join(., metadata, by = "metadata_code") %>% # aqui poner el codigo para el food y water intake, tomar bout length (mean), # of bouts
                group_by(date, hr, ID, cage_number, DIET, DIET_CODE, KCAL_PER_GR,
                         SEX, COHORT, STRAIN, AIM, sable_idx, metadata_code,
                         parameter) %>% 
                group_split() %>% 
                map_dfr(., function(X){
                  patterns <- c("AllMeters_", "PedMeters_")
                  regex_pattern <- paste(patterns, collapse = "|")
                  if (any(grepl(regex_pattern, X$parameter[1]))){
                    out <- X %>% 
                      group_by(date, hr, ID, cage_number, DIET, DIET_CODE, KCAL_PER_GR,
                               SEX, COHORT, STRAIN, AIM, sable_idx, metadata_code,
                               parameter) %>% 
                      summarise(
                        value = max(value), .groups = "drop_last"
                      )}
                    else{
                      out <- X %>% 
                        group_by(date, hr, ID, cage_number, DIET, DIET_CODE, KCAL_PER_GR,
                                 SEX, COHORT, STRAIN, AIM, sable_idx, metadata_code,
                                 parameter) %>% 
                        summarise(
                          value = mean(value, na.rm = TRUE), .groups = "drop_last"
                        )
                    }
                    return(out)
                })
        }, .progress = TRUE
    )
saveRDS(sable_hr_data, file = "../data/sable/sable_hr_data.rds", compress = TRUE)
sable_hr_data <- readRDS("../data/sable/sable_hr_data.rds")
# before/after injection ----

## helper functions ----

# this function works for cumulative data such as allmeters
create_corrected_updata <- function(X){
    corrected_values <- X %>% 
        ungroup() %>% 
        mutate(
            lag_series = replace_na(value - dplyr::lag(value, n=1), 0),
            flag_event = cumsum(replace_na(if_else(lag_series < 0, 1, 0), 0)),
            corrected_value = value + cumsum(if_else(lag_series < 0, abs(lag_series), 0))
            )
    return(corrected_values)
}
# this functions works for data that always goes down, such as food intake
create_corrected_downdata <- function(X){
    corrected_values <- X %>% 
        ungroup() %>% 
        mutate(
            lag_series = replace_na(value - dplyr::lag(value, n=1), 0),
            flag_event = cumsum(replace_na(if_else(lag_series > 0, 1, 0), 0)),
            corrected_value = value - cumsum(if_else(lag_series > 0, abs(lag_series), 0))
            )
    return(corrected_values)
}
# selects the data before and after an injection time
time_window <- function(injection_time, window_size_hr, data){
    create_datetime <- lubridate::ymd_hms(injection_time)
    time_before_injection <- create_datetime - lubridate::hours(window_size_hr)
    time_after_injection <- create_datetime + lubridate::hours(window_size_hr)
    data_with_dttm <- data %>% 
        mutate(datetime = lubridate::ymd_hms(paste(date, " ", hr, ":00:00", sep = "")),
               event_flag = case_when(
                   datetime >= time_before_injection & datetime <=create_datetime ~ "before",
                   datetime >= create_datetime & datetime <= time_after_injection ~ "after",
                   TRUE ~ NA
               )) %>% 
        drop_na(event_flag)
    return(data_with_dttm)
}

## before/after analysis ----

# get injections metadata
injection_time <- read_csv("../data/META_INJECTIONS.csv") %>% 
    # set injection time into dttm
    mutate(
        INJECTION_TIME = lubridate::mdy_hm(INJECTION_TIME)
    ) %>% 
    group_by(row_number()) %>% 
    group_split()

# get food intake, body weight and locomotion
# sable hr data is the hourly data needed for this analysis
before_after_data <- sable_hr_data %>% 
    filter(grepl("FoodA_*|BodyMass_*|AllMeters_*", parameter)) %>% 
    group_by(str_remove(parameter, "_[0-9]+")) %>% 
    group_split()
before_after_data

# this grid maps each dataset (food, bw and meter) to each one of the injections
data_injection_grid <- expand_grid(
    before_after_idx = 1:length(before_after_data),
    injection_time_idx = 1:length(injection_time)
)

# feed the grid into the map
before_after_analysis <- data_injection_grid %>% 
    group_by(row_number()) %>% 
    group_split() %>% 
    map_dfr(., function(X){
        # filters the data by the id in the injection metadata
        filtered_data <- before_after_data[[X$before_after_idx]] %>% 
            filter(ID == injection_time[[X$injection_time_idx]]$ID[1])
        # if there's no match just return the filtered data
        if (nrow(filtered_data) * ncol(filtered_data) == 0){
            return(filtered_data)
        }
        else{
            # find out which kind of data this is
            data_type <- str_remove(before_after_data[[X$before_after_idx]]$parameter[1], "_[0-9]+")
            # do different correction according to data type
            corrected_data <- switch(data_type,
                FoodA = filtered_data %>% create_corrected_downdata(),
                AllMeters = filtered_data %>% create_corrected_updata(),
                BodyMass = filtered_data %>% mutate(corrected_value = value),
                print("ERROR")
            )
            # select here the number of hours for the time window
            time_window_data <- time_window(
                injection_time[[X$injection_time_idx]]$INJECTION_TIME[1],
                12, # change this
                corrected_data
            ) %>% mutate(drug = injection_time[[X$injection_time_idx]]$DRUG[1])
            # return the corrected data values
            return(time_window_data)
        }
    })
saveRDS(before_after_analysis, file = "../data/sable/before_after_analysis.rds", compress = TRUE)
before_after_analysis <- readRDS("../data/sable/before_after_analysis.rds")

before_after_analysis %>% 
  filter(ID == 1006) %>% 
  mutate(parameter = str_remove(parameter, "_[0-9]+")) %>% 
  ggplot(aes(datetime, corrected_value, color = event_flag)) +
  geom_point() +
  geom_line() +
  facet_wrap(parameter~ID*drug, scales = "free")

before_after_analysis %>% 
  group_by(parameter, event_flag, drug, ID) %>%
  summarise(
    delta = abs(max(corrected_value)-min(corrected_value))
  ) %>% 
  filter(grepl("FoodA_", parameter)) %>% 
  ggplot(aes(interaction(drug, event_flag), delta, color = as.factor(ID))) +
  geom_point() +
  geom_line(aes(group = ID)) +
  facet_wrap(~drug, scale = "free_x")
