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
raw_data = ReadCSV(paste0(data_path, "/1.5_Expansion_Cytoplasm.csv"))

# Filter cells according to nuclear area
raw_data <- raw_data[-which(raw_data$`Nucleus: Area µm^2` < 8),]
raw_data <- raw_data[-which(raw_data$`Nucleus: Area µm^2` > 135),]

# Normalize marker MFIs from 0-1, clip maximum MFI at 99.7th percentile
SOX10_quar = as.numeric(quantile(raw_data$SOX10, c(.997)))
SOX10_min = as.numeric(min(raw_data$SOX10))
for(row in 1:nrow(raw_data)){
  SOX10 <- raw_data$SOX10[row]
  if(SOX10 > SOX10_quar) {
    raw_data$SOX10[row] <- as.numeric(1.0)
  } else if(SOX10 < SOX10_quar) {
    raw_data$SOX10[row] <- as.numeric((raw_data$SOX10[row] - SOX10_min) / (SOX10_quar - SOX10_min))
  }
}

CD45_quar = as.numeric(quantile(raw_data$CD45, c(.997)))
CD45_min = as.numeric(min(raw_data$CD45))
for(row in 1:nrow(raw_data)){
  CD45 <- raw_data$CD45[row]
  if(CD45 > CD45_quar) {
    raw_data$CD45[row] <- as.numeric(1.0)
  } else if(CD45 < CD45_quar) {
    raw_data$CD45[row] <- as.numeric((raw_data$CD45[row] - CD45_min) / (CD45_quar - CD45_min))
  }
}

CD3_quar = as.numeric(quantile(raw_data$CD3, c(.997)))
CD3_min = as.numeric(min(raw_data$CD3))
for(row in 1:nrow(raw_data)){
  CD3 <- raw_data$CD3[row]
  if(CD3 > CD3_quar) {
    raw_data$CD3[row] <- as.numeric(1.0)
  } else if(CD3 < CD3_quar) {
    raw_data$CD3[row] <- as.numeric((raw_data$CD3[row] - CD3_min) / (CD3_quar - CD3_min))
  }
}

CD4_quar = as.numeric(quantile(raw_data$CD4, c(.997)))
CD4_min = as.numeric(min(raw_data$CD4))
for(row in 1:nrow(raw_data)){
  CD4 <- raw_data$CD4[row]
  if(CD4 > CD4_quar) {
    raw_data$CD4[row] <- as.numeric(1.0)
  } else if(CD4 < CD4_quar) {
    raw_data$CD4[row] <- as.numeric((raw_data$CD4[row] - CD4_min) / (CD4_quar - CD4_min))
  }
}

CD8_quar = as.numeric(quantile(raw_data$CD8, c(.997)))
CD8_min = as.numeric(min(raw_data$CD8))
for(row in 1:nrow(raw_data)){
  CD8 <- raw_data$CD8[row]
  if(CD8 > CD8_quar) {
    raw_data$CD8[row] <- as.numeric(1.0)
  } else if(CD8 < CD8_quar) {
    raw_data$CD8[row] <- as.numeric((raw_data$CD8[row] - CD8_min) / (CD8_quar - CD8_min))
  }
}

FOXP3_quar = as.numeric(quantile(raw_data$FOXP3, c(.997)))
FOXP3_min = as.numeric(min(raw_data$FOXP3))
for(row in 1:nrow(raw_data)){
  FOXP3 <- raw_data$FOXP3[row]
  if(FOXP3 > FOXP3_quar) {
    raw_data$FOXP3[row] <- as.numeric(1.0)
  } else if(FOXP3 < FOXP3_quar) {
    raw_data$FOXP3[row] <- as.numeric((raw_data$FOXP3[row] - FOXP3_min) / (FOXP3_quar - FOXP3_min))
  }
}

CD20_quar = as.numeric(quantile(raw_data$CD20, c(.997)))
CD20_min = as.numeric(min(raw_data$CD20))
for(row in 1:nrow(raw_data)){
  CD20 <- raw_data$CD20[row]
  if(CD20 > CD20_quar) {
    raw_data$CD20[row] <- as.numeric(1.0)
  } else if(CD20 < CD20_quar) {
    raw_data$CD20[row] <- as.numeric((raw_data$CD20[row] - CD20_min) / (CD20_quar - CD20_min))
  }
}

CD11C_quar = as.numeric(quantile(raw_data$CD11C, c(.997)))
CD11C_min = as.numeric(min(raw_data$CD11C))
for(row in 1:nrow(raw_data)){
  CD11C <- raw_data$CD11C[row]
  if(CD11C > CD11C_quar) {
    raw_data$CD11C[row] <- as.numeric(1.0)
  } else if(CD11C < CD11C_quar) {
    raw_data$CD11C[row] <- as.numeric((raw_data$CD11C[row] - CD11C_min) / (CD11C_quar - CD11C_min))
  }
}

F4_80_quar = as.numeric(quantile(raw_data$`F4/80`, c(.997)))
F4_80_min = as.numeric(min(raw_data$`F4/80`))
for(row in 1:nrow(raw_data)){
  F4_80 <- raw_data$`F4/80`[row]
  if(F4_80 > F4_80_quar) {
    raw_data$`F4/80`[row] <- as.numeric(1.0)
  } else if(F4_80 < F4_80_quar) {
    raw_data$`F4/80`[row] <- as.numeric((raw_data$`F4/80`[row] - F4_80_min) / (F4_80_quar - F4_80_min))
  }
}

CD68_quar = as.numeric(quantile(raw_data$CD68, c(.997)))
CD68_min = as.numeric(min(raw_data$CD68))
for(row in 1:nrow(raw_data)){
  CD68 <- raw_data$CD68[row]
  if(CD68 > CD68_quar) {
    raw_data$CD68[row] <- as.numeric(1.0)
  } else if(CD68 < CD68_quar) {
    raw_data$CD68[row] <- as.numeric((raw_data$CD68[row] - CD68_min) / (CD68_quar - CD68_min))
  }
}

CD86_quar = as.numeric(quantile(raw_data$CD86, c(.997)))
CD86_min = as.numeric(min(raw_data$CD86))
for(row in 1:nrow(raw_data)){
  CD86 <- raw_data$CD86[row]
  if(CD86 > CD86_quar) {
    raw_data$CD86[row] <- as.numeric(1.0)
  } else if(CD86 < CD86_quar) {
    raw_data$CD86[row] <- as.numeric((raw_data$CD86[row] - CD86_min) / (CD86_quar - CD86_min))
  }
}

CD163_quar = as.numeric(quantile(raw_data$CD163, c(.997)))
CD163_min = as.numeric(min(raw_data$CD163))
for(row in 1:nrow(raw_data)){
  CD163 <- raw_data$CD163[row]
  if(CD163 > CD163_quar) {
    raw_data$CD163[row] <- as.numeric(1.0)
  } else if(CD163 < CD163_quar) {
    raw_data$CD163[row] <- as.numeric((raw_data$CD163[row] - CD163_min) / (CD163_quar - CD163_min))
  }
}

CD206_quar = as.numeric(quantile(raw_data$CD206, c(.997)))
CD206_min = as.numeric(min(raw_data$CD206))
for(row in 1:nrow(raw_data)){
  CD206 <- raw_data$CD206[row]
  if(CD206 > CD206_quar) {
    raw_data$CD206[row] <- as.numeric(1.0)
  } else if(CD206 < CD206_quar) {
    raw_data$CD206[row] <- as.numeric((raw_data$CD206[row] - CD206_min) / (CD206_quar - CD206_min))
  }
}

NK1.1_quar = as.numeric(quantile(raw_data$NK1.1, c(.997)))
NK1.1_min = as.numeric(min(raw_data$NK1.1))
for(row in 1:nrow(raw_data)){
  NK1.1 <- raw_data$NK1.1[row]
  if(NK1.1 > NK1.1_quar) {
    raw_data$NK1.1[row] <- as.numeric(1.0)
  } else if(NK1.1 < NK1.1_quar) {
    raw_data$NK1.1[row] <- as.numeric((raw_data$NK1.1[row] - NK1.1_min) / (NK1.1_quar - NK1.1_min))
  }
}

CD31_quar = as.numeric(quantile(raw_data$CD31, c(.997)))
CD31_min = as.numeric(min(raw_data$CD31))
for(row in 1:nrow(raw_data)){
  CD31 <- raw_data$CD31[row]
  if(CD31 > CD31_quar) {
    raw_data$CD31[row] <- as.numeric(1.0)
  } else if(CD31 < CD31_quar) {
    raw_data$CD31[row] <- as.numeric((raw_data$CD31[row] - CD31_min) / (CD31_quar - CD31_min))
  }
}

# Add column for numeric id assignment
raw_data$obj_id <- seq.int(nrow(raw_data))

# Export min-max normalized data
write.csv(raw_data, file = paste0(getwd(), "/MinMax Normalized.csv"))
