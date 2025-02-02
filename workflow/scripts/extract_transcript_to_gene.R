#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: extract_transcript_to_gene.R <input_gtf> <output_csv>")
}

input_gtf <- args[1]
output_csv <- args[2]

# Install and load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(rtracklayer)
library(dplyr)

# Load GTF file using rtracklayer::import
gtf <- rtracklayer::import(input_gtf, format = "gtf")

# Check if 'transcript_id' and 'gene_id' exist in metadata columns
if (!"transcript_id" %in% colnames(mcols(gtf)) || !"gene_id" %in% colnames(mcols(gtf))) {
  stop("Error: The required columns 'transcript_id' or 'gene_id' are not present in the GTF file.")
}

# Extract transcript and gene information
transcript_to_gene <- as.data.frame(mcols(gtf)) %>%
  select(transcript_id, gene_id) %>%
  distinct()

# Remove any rows with NA values
transcript_to_gene <- na.omit(transcript_to_gene)

# Save the transcript-to-gene mapping to CSV
write.csv(transcript_to_gene, output_csv, row.names = FALSE)
cat(paste("Transcript-to-gene mapping saved to", output_csv, "\n"))
