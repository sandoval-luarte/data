#The aim of this script is to make quality control of the raw data sable
#This means that the number of exp files within the folder called backup_allsable should be the same as sable_raw
#The folder backup_allsable contains all historical data in .exp .exd and .exp.bak whereas sable_raw only has .exp data 


# Define the paths to your specific folders####
folder1 <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/backup_allsable"  # Folder 1 (backup_allsable)
folder2 <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/sable_raw"  # Folder 2 (sable_raw)

# Get the list of files in both folders####
files_folder1 <- list.files(folder1)
files_folder2 <- list.files(folder2)

# Find missing files in folder1 (files that are in folder2 but not in folder1)
missing_in_folder1 <- setdiff(files_folder2, files_folder1)

# Find missing files in folder2 (files that are in folder1 but not in folder2)
missing_in_folder2 <- setdiff(files_folder1, files_folder2)

# Print the results####
cat("Files missing in folder1 (backup_allsable):\n")
print(missing_in_folder1)

cat("\nFiles missing in folder2 (sable_raw):\n")
print(missing_in_folder2)

#The second step is to know if we have the same amount of files in sable_raw and csv files####

# Define the paths to both folders
sable_raw_path <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/Sable_raw"
csv_files_path <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/Sable_raw/csv_files"

# Get the number of files in each folder
num_sable_raw <- length(list.files(sable_raw_path, pattern = "\\.exp$", full.names = TRUE))
num_csv_files <- length(list.files(csv_files_path, pattern = "\\.csv$", full.names = TRUE))

# Print the results
cat("Number of .exp files in Sable_raw:", num_sable_raw, "\n")
cat("Number of .csv files in csv_files:", num_csv_files, "\n")

# Check if they have the same amount
if (num_sable_raw == num_csv_files) {
  cat("✅ Both folders have the same number of files.\n")
} else {
  cat("⚠️ The number of files is different!\n")
}
# The third step is to identify which files are not in csv_files folder####
# Define folder paths
sable_raw_path <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/Sable_raw"
csv_files_path <- "/Users/carosandovalcab/Library/CloudStorage/GoogleDrive-clsandov@umn.edu/My Drive/PostDoc_Data/Sable_raw/csv_files"

# Get list of .exp files (remove extension)
exp_files <- list.files(sable_raw_path, pattern = "\\.exp$", full.names = FALSE)
exp_files_base <- sub("\\.exp$", "", exp_files)  # Remove .exp extension

# Get list of .csv files (remove suffix and extension)
csv_files <- list.files(csv_files_path, pattern = "-extd_m\\.csv$", full.names = FALSE)
csv_files_base <- sub("-extd_m\\.csv$", "", csv_files)  # Remove -extd_m.csv suffix

# Find missing files in csv_files
missing_in_csv <- setdiff(exp_files_base, csv_files_base)

# Print missing files
cat("Files missing in csv_files:\n")
print(missing_in_csv)

