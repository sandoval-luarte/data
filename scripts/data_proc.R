# libs ----

pacman::p_load(
  tidyverse,
  googledrive,
  furrr,
  zoo
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
metadata <- read_csv("../data/META.csv") %>% 
    pivot_longer(cols = starts_with("SABLE_DAY_"),
                 names_to = "sable_idx",
                 values_to = "metadata_code")

## esto es solo para trabajar en la funcion para detectar bouts
x <- read_csv(sable_csv_files[[1]])


## merge metadata with sable data in an hourly fashion

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
                    metadata_code = paste(format(date, "%m/%d/%Y"),
                                          cage_number, sep = "-")
                ) %>% 
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
                          value = mean(value), .groups = "drop_last"
                        )
                    }
                    return(out)
                })
        }, .progress = TRUE
    )

saveRDS(sable_hr_data, file = "../data/sable/sable_hr_data.rds", compress = TRUE)

## esto es solo para probar en un solo csv
test <- x %>%
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
        metadata_code = paste(format(date, "%m/%d/%Y"),
                              cage_number, sep = "-")
    ) %>% 
    filter(grepl("FoodA_", parameter)) %>% 
    group_by(cage_number, date, parameter) %>% 
    group_split()

## bt es 1 csv al cual le hice el procesamiento de el rolling sd y la deteccion de bouts
bt <- test[[1]] %>% 
    mutate(
      rolling_sd_raw = zoo::rollapply(value, 10, fill = 0, FUN = sd, align = "right"),
      rolling_sd = as.numeric((zoo::rollapply(value, 10, fill = 0, FUN = sd, align = "right")) > 0.02),
      bouts = as.numeric(data.table::rleid(rolling_sd))
      ) %>% 
  group_by(bouts) %>% 
  mutate(bout_duration = n())

## esto es para computar las diferencias entre los periodos quietos
## saco el t test y todas las manos
bt_filter <- bt %>% 
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
    }


## aca devuelvo el analisis entre periodos quietos a la data original
## mediante un left join
## aqui aplico las correcciones para las ingestas
output <- bt %>% 
  left_join(., bt_filter, by = "bouts") %>%
  ungroup() %>% 
  group_by(bouts) %>% 
  mutate(
    corrected_value = if_else(rolling_sd == 0, mean(value), value), # para el resumen en 1h, considerar solo los momentos donde sd = 0, tomar promedio
    corrected_diff = if_else(diff < 0 | t_test_pval > 0.05, 0, diff)
  )


output %>% 
  filter(rolling_sd == 0) %>% 
  ggplot(aes(DateTime, corrected_value)) +
  geom_line() +
  geom_smooth(span = 0.8)

