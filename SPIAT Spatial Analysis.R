# Script for implementation of SPIAT tools in spatial analysis pipeline

# Load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(SPIAT)

# Set working directory
setwd("...")
data_path = getwd()

# Load centroid, intensity, and phenotype data
raw_data = read.csv(paste0(data_path, "/FlowSOM Phenotyping.csv"))

# Define dataframe for intensity matrix
marker_list <- c("SOX10", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD11C", "F4/80", "CD68", "CD86", "CD163", "CD206", "NK1.1", "CD31")
intensity <- data.frame(matrix(NA, nrow = nrow(raw_data), ncol = length(marker_list)))

colnames(intensity) <- marker_list
intensity$SOX10 = raw_data$SOX10
intensity$CD45 = raw_data$CD45
intensity$CD3 = raw_data$CD3
intensity$CD4 = raw_data$CD4
intensity$CD8 = raw_data$CD8
intensity$FOXP3 = raw_data$FOXP3
intensity$CD20 = raw_data$CD20
intensity$CD11C = raw_data$CD11C
intensity$`F4/80` = raw_data$F4.80
intensity$CD68 = raw_data$CD68
intensity$CD86 = raw_data$CD86
intensity$CD163 = raw_data$CD163
intensity$CD206 = raw_data$CD206
intensity$NK1.1 = raw_data$NK1.1
intensity$CD31 = raw_data$CD31

# Define intensity matrix (transpose), xy coordinates, and phenotypes
intensity_matrix <- data.frame(t(intensity))
coord_x <- c(raw_data$Centroid.X.µm)
coord_y <- c(raw_data$Centroid.Y.µm)
phenotypes <- c(raw_data$Phenotype)

############### SPIAT BEGINS ###############

spiat_object <- format_image_to_spe(format = "general", 
                                    intensity_matrix = intensity_matrix,
                                    phenotypes = phenotypes, 
                                    coord_x = coord_x, coord_y = coord_y)

# Optional: Plot heatmap for different markers (note y-axis is inverted, so image is flipped)
plot_marker_level_heatmap(spiat_object, num_splits = 100, "SOX10")

# Define cell types
formatted_image <- define_celltypes(
  spiat_object, 
  categories = c("Tumor Cells", "CD4 T Cells", "CD8 T Cells", "Tregs", "M1 Macrophages", "M2 Macrophages", "Other Macrophages", "Dendritic Cells", "NK Cells", "Endothelial Cells", "Other Immune"),
  category_colname = "Phenotype", 
  names = c("Tumor Cells", "CD4 T Cells", "CD8 T Cells", "Tregs", "M1 Macrophages", "M2 Macrophages", "Other Macrophages", "Dendritic Cells", "NK Cells", "Endothelial Cells", "Other Immune"),
  new_colname = "Cell.Type")

# Average minimum cell distance calculation
min_dist <- calculate_minimum_distances_between_celltypes(
  spe_object = formatted_image, 
  cell_types_of_interest = c("Tumor Cells", "CD4 T Cells", "CD8 T Cells", "Tregs", "M1 Macrophages", "M2 Macrophages", "Other Macrophages", "Dendritic Cells", "NK Cells", "Endothelial Cells", "Other Immune"),
  feature_colname = "Cell.Type")

# Export average minimum distance data
min_summary_dist <- calculate_summary_distances_between_celltypes(min_dist)
write.csv(min_summary_dist, "Minimum Distances.csv")

# Calculate percentage breakdowns of cell types around one phenotype; radius can be adjusted
average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells", 
                                          target_celltype = "Tumor Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells", 
                                          target_celltype = "CD4 T Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells", 
                                          target_celltype = "CD8 T Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "Tregs", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells", 
                                          target_celltype = "M1 Macrophages", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "M2 Macrophages", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "Other Macrophages", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "Dendritic Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "NK Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "Other Immune Cells", 
                                          radius = 100, feature_colname="Cell.Type")

average_percentage_of_cells_within_radius(spe_object = formatted_image, 
                                          reference_celltype = "CD8 T Cells",
                                          target_celltype = "Endothelial Cells", 
                                          radius = 100, feature_colname="Cell.Type")

# ... Repeat for all relevant cell types

# Calculate normalized mixing scores (NMS); radius can be adjusted
mixing_score_summary(spe_object = formatted_image, reference_celltype = "Tumor Cells", 
                     target_celltype = "CD4 T Cells", radius=100, feature_colname ="Cell.Type")

mixing_score_summary(spe_object = formatted_image, reference_celltype = "Tumor Cells", 
                     target_celltype = "CD8 T Cells", radius=100, feature_colname ="Cell.Type")

mixing_score_summary(spe_object = formatted_image, reference_celltype = "CD8 T Cells", 
                     target_celltype = "M1 Macrophages", radius=100, feature_colname ="Cell.Type")

mixing_score_summary(spe_object = formatted_image, reference_celltype = "CD8 T Cells", 
                     target_celltype = "M2 Macrophages", radius=100, feature_colname ="Cell.Type")

mixing_score_summary(spe_object = formatted_image, reference_celltype = "CD8 T Cells", 
                     target_celltype = "Other Macrophages", radius=100, feature_colname ="Cell.Type")

# ... Repeat for all relevant cell types