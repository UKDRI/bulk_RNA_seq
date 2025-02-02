#!/usr/bin/env Rscript

library(ggplot2)
library(pheatmap)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
samplesheet_file <- args[2]
output_pca <- args[3]
output_corr <- args[4]

# Ensure output directories exist
dir.create(dirname(output_pca), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_corr), recursive = TRUE, showWarnings = FALSE)

# Load expression data
expression_data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Load samplesheet and ensure it contains 'sample_id' and 'group' columns
samplesheet <- read.csv(samplesheet_file)
if (!all(c("sample_id", "group") %in% colnames(samplesheet))) {
    stop("Samplesheet must contain 'sample_id' and 'group' columns.")
}

# Remove rows with all zeros
expression_data <- expression_data[rowSums(expression_data != 0) > 0, ]

# Transpose data for PCA
expression_data_t <- t(expression_data)

# Perform PCA
pca <- prcomp(expression_data_t, scale. = TRUE)

# Create PCA data frame
pca_df <- as.data.frame(pca$x)
pca_df$sample_id <- rownames(pca_df)

# Merge PCA data with group information from the samplesheet
pca_df <- merge(pca_df, samplesheet, by = "sample_id", all.x = TRUE)

# PCA Plot using the 'group' column from the samplesheet
# Generate a dynamic color palette for groups
group_colors <- setNames(
  RColorBrewer::brewer.pal(n = length(unique(pca_df$group)), name = "Set2"),
  unique(pca_df$group)
)

# PCA Plot using the dynamically assigned colors
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample_id)) +
    geom_point(size = 5, alpha = 0.8) +  # Solid, filled points
    geom_text(
        aes(label = sample_id),
        hjust = 1.1, vjust = -0.5, size = 4, 
        color = "black", fontface = "bold"
    ) +  # Improved label placement
    theme_minimal(base_size = 16) +
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),  # Larger axis tick labels
        axis.title = element_text(size = 16),  # Larger axis titles
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels
        axis.text.y = element_text(size = 14),  # Larger y-axis labels
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.margin = margin(10, 10, 10, 10)  # Add margin space
    ) +
    labs(
        title = "Principal Component Analysis (PCA)",
        x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 2), "% variance)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 2), "% variance)")
    ) +
    scale_color_manual(values = group_colors)

# Save PCA plot
ggsave(output_pca, plot = pca_plot, width = 12, height = 8)  # Adjust plot dimensions for better readability



# Correlation Matrix
cor_matrix <- cor(expression_data_t)

# Professional correlation matrix heatmap
pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 10,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Sample Correlation Matrix",
  fontsize = 12,
  filename = output_corr
)
