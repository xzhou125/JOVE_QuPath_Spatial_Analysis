# Script for phenotyping proteomics data using Seurat

# Install dplyr, tidyverse, and ggplot2 packages
# Install Seurat (specifically v.4.4.0) package

# Load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Seurat)

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
raw_data = ReadCSV(paste0(data_path, "/MinMax Normalized.csv"))

# Format data for Seurat entry
CustomReadAkoya <- function(mtx){
  cell_intensity_df <- data.frame(matrix(ncol = length(marker_list), nrow = length(mtx$obj_id)))
  colnames(cell_intensity_df)<-marker_list
  rownames(cell_intensity_df)<-mtx$obj_id
  
  for(marker in marker_list) {
    message(marker)
    cell_intensity_df[marker]<-mtx[[paste0(marker)]]
  }
  
  expression_matrix <- t(cell_intensity_df)
  expression_matrix <- data.frame(expression_matrix)
  colnames(expression_matrix)<-mtx$obj_id
  
  expression_matrix <- as.sparse(x = expression_matrix)
  
  centroids <- data.frame(
    x = mtx["Centroid X µm"],
    y = mtx["Centroid Y µm"],
    cell = as.character(x = mtx["obj_id"]),
    stringAsFactors = FALSE
  )
  
  meta <- centroids
  rownames(meta) <- mtx$obj_id
  
  list(matrix = expression_matrix, centroids = centroids, metadata = meta)
}

data <- CustomReadAkoya(mtx = raw_data)

# Define Seurat object
seurat_obj <- CreateSeuratObject(
  counts = data$matrix,
  assay = "Akoya",
  meta.data = data$metadata,
  project = '...',
  min.cells = 3,
  min.features = 1,
)

save_name = paste0(getwd(), "/Cytoplasm_Seurat.Robj")
save(seurat_obj, file = save_name)

# Begin clustering process and UMAP/dot plot generation
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), vars.to.regress = "nFeature_Akoya")
features = marker_list

# Define npcs as the greatest value in which the results are still valid (i.e., converge)
seurat_obj <- RunPCA(object = seurat_obj, npcs = 13, features = features, verbose = FALSE)
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:13)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.95) # Lower resolution gives fewer clusters
seurat_obj <- RunUMAP(seurat_obj, dims = 1:13, n.neighbors = 50L, min.dist = 0.3)

# Plot clustering results
Idents(seurat_obj) = "seurat_clusters"
DimPlot(seurat_obj, reduction = "umap", label = TRUE, raster = FALSE)
DotPlot(NormalizeData(seurat_obj), features = features) + RotatedAxis()
averages <- AverageExpression(NormalizeData(seurat_obj), features = features, return.seurat = TRUE)
DoHeatmap(averages, features = features) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

# Create new dataframe for phenotype assignments
phenotypes = data.frame(matrix(NA, nrow = nrow(raw_data), ncol = 2))
df_headings <- c("seurat_cluster", "phenotype")
colnames(phenotypes) <- df_headings
phenotypes$seurat_cluster <- seurat_obj$seurat_clusters

# Phenotype seurat clusters using the heatmap of marker expression profiles
for(row in 1:nrow(phenotypes)){
  cluster <- phenotypes$seurat_cluster[row]
  if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Tumor Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "CD4 T Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "CD8 T Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Tregs"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Other T Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "M1 Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "M2 Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Other Macrophages"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Dendritic Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "NK Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Endothelial Cells"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Other Immune"
  } else if(cluster == 0 | cluster == 0) {
    phenotypes$phenotype[row] = "Remove"
  }
}

# Add cluster and phenotype information to original data table
raw_data$cluster <- phenotypes$seurat_cluster
raw_data$phenotype <- phenotypes$phenotype

# Export new data table as csv
filename <- paste0(data_path, "/Seurat Phenotyping.csv")
write.csv(raw_data, filename, row.names = FALSE)