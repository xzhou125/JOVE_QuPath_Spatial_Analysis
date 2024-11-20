# Script for phenotyping proteomics data using FlowSOM

# Load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(janitor)
library(flowCore)
library(FlowSOM)
library(pheatmap)

# Set working directory and define list of PhenoCycler panel lineage markers (i.e., markers for clustering)
setwd("...")
data_path = getwd()
marker_list <- c("SOX10", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD20", "CD11C", "F4/80", "CD68", "CD86", "CD163", "CD206", "NK1.1", "CD31")

# Function to import csv files
ReadCSV <- function(file) {
  data = read_csv(file, show_col_types = FALSE)
  return(data)
}

# Load in csv file(s)
data = ReadCSV(paste0(data_path, "/MinMax Normalized.csv"))

# Define columns and metadata for FlowSOM input
columns <- marker_list
metadata <- list(
  "id" = data$obj_id,
  "x_c" = data$`Centroid X µm`,
  "y_c" = data$`Centroid Y µm`
)

# Define FlowSOM Function
run_and_plot_flowSOM <- function(data, columns, markers, meta_list, maxM, xdim, ydim, comp = FALSE, trans = FALSE, scl = FALSE) {
  # Create flowFrame object
  ff <- flowFrame(as.matrix(data %>% select(all_of(columns))))
  
  # Add metadata to the FlowFrame
  for (key in names(meta_list)) {
    metadata_vector <- as.character(meta_list[[key]])
    keyword(ff)[[key]] <- metadata_vector
  }
  
  # Perform FlowSOM analysis
  fSOM <- FlowSOM(
    ff,
    colsToUse = columns,
    compensate = comp,
    transform = trans,
    scale = scl,
    maxMeta = maxM,
    xdim = xdim,
    ydim = ydim
  )
  
  # Plot minimal spanning tree
  plot <- PlotStars(
    fSOM,
    markers = markers,
    backgroundValues = fSOM$metaclustering,
    equalNodeSize = TRUE
  )
  
  return(list("fSOM" = fSOM, "plot" = plot))
}

# Run FlowSOM function
fSOM <- run_and_plot_flowSOM(
  data = data,
  columns = columns,
  markers = columns,
  meta_list = metadata,
  maxM = 50,
  xdim = 10,
  ydim = 10
)

# Generate FlowSOM clusters heatmap AND EXPORT RESULTS
### (IMPORTANT - Rerunning FlowSOM will generate different clustering)
node_means <- fSOM$fSOM$map$codes
rownames(node_means) <- as.character(1:nrow(node_means))
pheatmap(t(node_means))

# Export summary of FlowSOM results
FlowSOMmary(
  fsom = fSOM$fSOM,
  plotFile = paste0(data_path, "/FlowSOMmary.pdf")
)

# Add FlowSOM cluster and metacluster classifications to original csv
fsom_nodes = fSOM$fSOM$map$mapping[,1]
fsom_metaclusters = GetMetaclusters(fsom = fSOM$fSOM)
data_output <- data
data_output$cluster <- fsom_nodes
data_output$metaclusters <- fsom_metaclusters

# Export csv file with cluster classifications
### (IMPORTANT - Rerunning FlowSOM will generate different clustering)
filename <- paste0(data_path, "/With FlowSOM Clusters.csv")
write.csv(data_output, filename, row.names = FALSE)

# To preserve FlowSOM cluster assignments for subsequent phenotyping runs, load csv file containing cluster information
# data_cluster = ReadCSV(paste0(data_path, "/With FlowSOM Clusters.csv"))

# Assign a Phenotype column to the dataframe
data_cluster = data_output
data_cluster$Phenotype = "Remove"

# Assign clusters to phenotypes (consult heatmap generated from FlowSOM)
for(row in 1:nrow(data_cluster)){
  cluster <- data_cluster$cluster[row]
  if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Tumor Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "CD4 T Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "CD8 T Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Tregs"
  } else if(cluster == 3 | cluster == 22) {
    data_cluster$Phenotype[row] = "Other T Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "M1 Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "M2 Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Other Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Dendritic Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "NK Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Endothelial Cells"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Other Immune"
  } else if(cluster == 0 | cluster == 0) {
    data_cluster$Phenotype[row] = "Remove"
  }
}

# Export phenotype results in csv
filename <- paste0(data_path, "/FlowSOM Phenotyping.csv")
write.csv(data_cluster, filename, row.names = FALSE)
