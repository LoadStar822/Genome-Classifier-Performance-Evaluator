install.packages("plyr")
library(plyr)

# Get the list of all RDS files in the current directory
files <- list.files(pattern = "\\.rds$")

# Loop over the files
for(file in files) {
  # Define the name of the CSV file
  csv_file <- sub("\\.rds$", ".csv", file)

  # Check if the CSV file already exists
  if(!file.exists(csv_file)) {
    # Read the RDS file
    data <- readRDS(file)

    # Convert the data to a data frame and write to CSV
    df_list <- lapply(data, function(x) data.frame(t(unlist(x))))
    df <- plyr::rbind.fill(df_list)
    write.csv(df, file = csv_file, row.names = FALSE)

    print(paste("Created", csv_file))
  } else {
    print(paste("Skipped", csv_file, "because it already exists"))
  }
}