
# Load necessary libraries
library(tidyverse)
library(here)

# Function to load CSV, skip metadata rows, remove first three data rows, and fill down blank cells
fill_down_csv <- function(input_file, output_file, columns_to_fill) {
  # Read the CSV file, skipping the first few metadata rows (adjust `skip` as needed)
  df <- read_csv(input_file, skip = 3, show_col_types = FALSE)  # Adjust skip value if needed
  
  # Remove the first three rows of the actual data if necessary
  df <- df[-c(1:3), ]
  
  # Check and print column names for verification
  cat("Column names in the data:\n")
  print(names(df))
  
  # Adjust column names in `columns_to_fill` if necessary based on the loaded data
  # columns_to_fill <- c("...")  # Uncomment and adjust if actual column names differ
  
  # Apply fill-down to specified columns
  df <- df %>%
    fill(all_of(columns_to_fill), .direction = "down")
  
  # Save the modified data frame to a new CSV file
  write_csv(df, output_file)
  
  cat("\nFilled data saved to:", output_file, "\n")
}

# Define the input and output file paths
input_file <- here("data/master_ioc_list_v11.2.csv")    # Replace with your input CSV path
output_file <- here("data/master_ioc_filled_v11.2.csv")  # Replace with your output CSV path

# Specify the columns to fill down - adjust these based on the actual column names in your data
columns_to_fill <- c("Order", "Family (Scientific)", "Family (English)", "Genus", "Species (Scientific)")

# Run the fill-down function
fill_down_csv(input_file, output_file, columns_to_fill)
view(df)
