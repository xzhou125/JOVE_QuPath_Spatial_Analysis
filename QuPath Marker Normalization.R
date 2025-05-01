# Script for processing proteomics data (from QuPath)

# Load libraries
library(tidyverse)

# Set working directory and define list of PhenoCycler panel lineage markers (i.e., markers for clustering)
setwd("...")
data_path = getwd()
marker_list <- c("SOX10", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD11C", "F4/80", "CD68", "CD86", "CD163", "CD206", "NK1.1", "CD31")

# Define function to import csv files
ReadCSV <- function(file) {
  data = read_csv(file, show_col_types = FALSE)
  return(data)
}

# Load in csv file(s)
raw_data = ReadCSV(paste0(data_path, "/measurements.csv"))

# Filter cells according to nuclear area
raw_data <- raw_data[-which(raw_data$`Nucleus: Area µm^2` < 8),]
raw_data <- raw_data[-which(raw_data$`Nucleus: Area µm^2` > 135),]

# Define the function for min-max normalization with clipping
normalize_with_clipping <- function(df, columns) {
  # Loop through the specified columns
  for (col in columns) {
    # Calculate the 99.7th percentile of the column
    percentile_997 <- quantile(df[[col]], 0.997, na.rm = TRUE)
    
    # Clip the values at the 99.7th percentile
    df[[col]] <- pmin(df[[col]], percentile_997)
    
    # Perform min-max normalization
    min_val <- min(df[[col]], na.rm = TRUE)
    max_val <- max(df[[col]], na.rm = TRUE)
    
    # Normalize the column (scaled to [0, 1])
    df[[col]] <- (df[[col]] - min_val) / (max_val - min_val)
  }
  
  return(df)
}

# Assign markers to normalize and run noramlization function
columns_to_normalize <- marker_list
df_normalized <- normalize_with_clipping(raw_data, columns_to_normalize)

# Add column for numeric id assignment
raw_data$obj_id <- seq.int(nrow(raw_data))

# Export min-max normalized data
write.csv(raw_data, file = paste0(getwd(), "/MinMax Normalized.csv"))
