# Install and load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(rtracklayer)
library(dplyr)

# Load GTF file using rtracklayer::import
gtf <- rtracklayer::import("data/genome/mouse/mouse_annotation.gtf", format = "gtf")

# Check if 'transcript_id' and 'gene_id' exist in metadata columns
if (!"transcript_id" %in% colnames(mcols(gtf)) || !"gene_id" %in% colnames(mcols(gtf))) {
  stop("The required columns 'transcript_id' or 'gene_id' are not present in the GTF file.")
}

# Extract transcript and gene information from the metadata columns
transcript_to_gene <- as.data.frame(mcols(gtf)) %>%
  dplyr::select(transcript_id, gene_id) %>%
  distinct()

# Remove any rows with NA values
transcript_to_gene <- na.omit(transcript_to_gene)

# Save the transcript-to-gene mapping to CSV
write.csv(transcript_to_gene, "../test-dataset/data/Quant/transcript_to_gene_mapping.csv", row.names = FALSE)
cat("Transcript-to-gene mapping saved to ../test-dataset/data/Quant/transcript_to_gene_mapping.csv\n")
