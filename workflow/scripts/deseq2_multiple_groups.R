#!/usr/bin/env Rscript

is_snakemake <- FALSE
# Check if we are running the script in Snakemake
if (is.object(snakemake)) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(file = log, type = "output")
    sink(file = log, type = "message")
    is_snakemake <- TRUE
}

# Load required libraries
required_packages <- c("DESeq2", "data.table", "optparse", "ggplot2", "dplyr", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "ggrepel")

install_libraries <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (pkg %in% c("DESeq2", "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db")) {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager", repos = "https://cloud.r-project.org")
                }
                BiocManager::install(pkg)
            } else {
                install.packages(pkg, repos = "https://cloud.r-project.org")
            }
        }
    }
}

# Parse command-line options
parse_args <- function() {
    option_list <- list(
        make_option(c("--expression"), type = "character", help = "Path to expression file"),
        make_option(c("--metadata"), type = "character", help = "Path to metadata file"),
        make_option(c("--comparisons"), type = "character", help = "Path to group comparisons CSV"),
        make_option(c("--canonicals"), type = "character", help = "Path to list of canonical transcripts IDs"),
        make_option(c("--species"), type = "character", help = "Species to get the gene names from, human or mouse"),
        make_option(c("--output_dir"), type = "character", help = "Directory to save results")
    )
    opt <- parse_args(OptionParser(option_list = option_list))

    # Check input files
    if (!file.exists(opt$expression)) stop("Error: Expression file not found.")
    if (!file.exists(opt$metadata)) stop("Error: Metadata file not found.")
    if (!file.exists(opt$comparisons)) stop("Error: Comparisons file not found.")
    if (!file.exists(opt$canonicals)) stop("Error: Canonical transcripts file not found.")
    if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
    if (!(opt$species == "human" | opt$species == "mouse")) stop("Only human or mouse are allowed for --species")

    return(opt)
}

deseq2_multiple_groups <- function(expression, metadata, comparisons, canonicals, species, output_dir) {
    # Load input data
    counts <- fread(expression, data.table = FALSE)
    metadata <- fread(metadata, data.table = FALSE)
    comparisons <- fread(comparisons, data.table = FALSE)
    canonical_transcript_ids <- fread(canonicals, data.table = FALSE, col.names=c("gene_id", "transcript_id"), quote="")

    # Ensure consistent samples
    colnames(counts)[1] <- "transcript_id"
    sample_ids <- intersect(colnames(counts)[-1], metadata$sample_id)
    if (length(sample_ids) < 2) stop("Error: Not enough matching samples between counts and metadata.")
    cat("Matching sample IDs:", paste(sample_ids, collapse = ", "), "\n")

    # Removing the version of the transcript if present
    counts$unversioned_id <- sapply(strsplit(counts$transcript_id,"\\."), `[`, 1)
    # Keeping only the canonical transcripts
    counts <- counts[counts$unversioned_id %in% canonical_transcript_ids$transcript_id,]
    mapping_transcript_gene <- setNames(canonical_transcript_ids$gene_id, canonical_transcript_ids$transcript_id)
    counts$gene_id <- mapping_transcript_gene[counts$unversioned_id]

    # Map Ensembl IDs to Gene Names
    if (species == "human") {
        gene_names <- AnnotationDbi::mapIds(
            org.Hs.eg.db,
            keys = counts$gene_id,
            column = "SYMBOL",
            keytype = "ENSEMBL",
            multiVals = "first"
        )
    } else {
        gene_names <- AnnotationDbi::mapIds(
            org.Mm.eg.db,
            keys = counts$gene_id,
            column = "SYMBOL",
            keytype = "ENSEMBL",
            multiVals = "first"
        )
    }

    counts$gene_name <- ifelse(is.na(gene_names), "Unknown", gene_names)

    # Subset counts and metadata
    counts <- counts[, c("gene_id", "gene_name", sample_ids), drop = FALSE]
    rownames(counts) <- counts$gene_id
    counts <- counts[, -c(1, 2)]  # Remove gene_id and gene_name columns for DESeq2 input
    metadata <- metadata[metadata$sample_id %in% sample_ids, ]
    rownames(metadata) <- metadata$sample_id

    # Ensure integer counts
    if (!all(apply(counts, 2, is.numeric))) stop("Error: Counts contain non-numeric values.")
    if (!all(counts == floor(counts))) {
        warning("Non-integer counts detected. Rounding to nearest integers.")
        counts <- round(counts)
    }

    # Perform DESeq2 analysis for each comparison
    for (i in 1:nrow(comparisons)) {
        treatment <- comparisons$treatment[i]
        control <- comparisons$control[i]
        comparison_name <- paste(treatment, "vs", control, sep = "_")

        # Filter metadata for the current comparison
        metadata_subset <- metadata[metadata$group %in% c(treatment, control), ]
        counts_subset <- counts[, rownames(metadata_subset), drop = FALSE]

        # Convert treatment column to factor
        metadata_subset$group <- factor(metadata_subset$group, levels = c(control, treatment))

        # DESeq2 analysis
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_subset), colData = metadata_subset, design = ~ group)
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("group", treatment, control))
        res <- as.data.frame(res)
        res$gene_id <- rownames(res)
        res$gene_name <- gene_names[res$gene_id]

        # Filter results: Remove NA or infinite values and order by padj
        res <- res %>%
            filter(!is.na(padj) & !is.infinite(log2FoldChange)) %>%
            mutate(
                regulation_status = case_when(
                    log2FoldChange > 0.6 & padj < 0.05 ~ "Upregulated",
                    log2FoldChange < -0.6 & padj < 0.05 ~ "Downregulated",
                    TRUE ~ "Not Significant"
                )
            ) %>%
            arrange(padj)  # Order by padj in ascending order

        # Save ordered results
        output_file <- file.path(output_dir, paste0(comparison_name, "_deg_results.csv"))
        write.csv(res, output_file, row.names = FALSE)

        # Volcano plot
    volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = regulation_status)) +
        geom_point(alpha = 0.7, size = 2) +  # Slightly larger points and better transparency
        geom_text_repel(
            data = res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1),  # Label only highly significant genes
            aes(label = gene_name),
            size = 3.5,
            max.overlaps = 10,  # Avoid too many overlapping labels
            box.padding = 0.3,
            point.padding = 0.3,
            segment.color = "grey50"
        ) +
        scale_color_manual(
            values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5", "Not Significant" = "grey70")
        ) +  # Use more vibrant and distinguishable colors
        theme_classic(base_size = 14) +  # Use a cleaner, larger base font size
        theme(
            axis.title = element_text(face = "bold"),  # Bold axis titles
            legend.position = "top",  # Move legend to the top for better space utilization
            legend.title = element_blank()  # Remove legend title for simplicity
        ) +
        labs(
            title = paste("Volcano Plot:", treatment, "vs", control),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"
        ) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +  # Dashed line for p=0.05
        geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", alpha = 0.7)  # Dotted lines for fold-change threshold

    # Save the plot
    ggsave(file.path(output_dir, paste0(comparison_name, "_volcano_plot.png")), plot = volcano, width = 8, height = 6)

    }
}

if (!is_snakemake) {
    install_libraries(required_packages)
}

for (pkg in required_packages) {
    suppressMessages(library(pkg, character.only = TRUE))
}

if (is_snakemake) {
    expression <- snakemake@input[["expression"]]
    canonicals <- snakemake@input[["canonicals"]]
    metadata <- snakemake@input[["metadata"]]
    comparisons <- snakemake@input[["comparisons_file"]]
    species <- snakemake@params[["species"]]
    output_dir <- snakemake@params[["out_dir"]]
} else {
    opt <- parse_args()
    expression <- opt$expression
    canonicals <- opt$canonicals
    metadata <- opt$metadata
    comparisons_file <- opt$comparisons
    species <- opt$species
    output_dir <- opt$output_dir
}

dir.create(c(output_dir), recursive = TRUE)
deseq2_multiple_groups(expression, metadata, comparisons, canonicals, species, output_dir)
cat("DESeq2 analysis completed successfully for all comparisons.\n")
