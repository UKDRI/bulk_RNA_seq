#!/usr/bin/env Rscript

is_snakemake <- FALSE
# Check if we are running the script in Snakemake
if (is.object(snakemake)) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(file = log, type = "output")
    sink(file = log, type = "message")
    is_snakemake <- TRUE
}

# Define necessary packages
required_packages <- c("clusterProfiler", "gprofiler2", "org.Hs.eg.db", "org.Mm.eg.db", "dplyr", "readr", "igraph", "optparse")

install_libraries <- function(packages) {
    # Load or install necessary packages
    options(repos = c(CRAN = "https://cran.r-project.org"))

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }

    # Install any missing packages
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            BiocManager::install(pkg, ask = FALSE)
        }
    }
}

# Load libraries after installation
for (pkg in required_packages) {
    suppressMessages(library(pkg, character.only = TRUE))
}

parse_args <- function() {
    # Parse arguments
    option_list <- list(
        make_option(c("-i"), type = "character", help = "Input file"),
        make_option(c("-o_bp"), type = "character", help = "Path to metadata file"),
        make_option(c("-o_mf"), type = "character", help = "Path to group comparisons CSV"),
        make_option(c("-o_cc"), type = "character", help = "Path to list of canonical transcripts IDs"),
        make_option(c("-org"), type = "character", help = "Species to get the gene names from, human or mouse"),
        make_option(c("-p"), type = "float", help = "p-value for cut-offs")
    )
    opt <- parse_args(OptionParser(option_list = option_list))
    if (!file.exists(opt$i)) {
        stop("The input file does not exist")
    }

    return(opt)
}

go_enrichment_analysis <- function(input_file, output_bp, output_cc, output_mf, organism, pvalue_cutoff) {
    # Read DESeq2 results; if missing or empty, create an empty data frame
    de_genes <- read_csv(input_file)
    if (nrow(de_genes) == 0) {
        warning("Input file is empty: ", input_file, ". Proceeding with empty DE results.")
        de_genes <- data.frame(gene_id = character(), log2FoldChange = numeric(), padj = numeric())
    } else {
        required_columns <- c("gene_id", "log2FoldChange", "padj")
        missing_columns <- setdiff(required_columns, colnames(de_genes))
        if (length(missing_columns) > 0) {
            stop("The input file is missing the following required columns: ",
                paste(missing_columns, collapse = ", "))
        }
    }

    # Filter upregulated and downregulated genes based on log2 fold change and p-value
    upregulated_genes <- de_genes %>%
        filter(log2FoldChange > 0 & padj < 0.05) %>%
        pull(gene_id)

    downregulated_genes <- de_genes %>%
        filter(log2FoldChange < 0 & padj < 0.05) %>%
        pull(gene_id)

    # Function to convert ENSG IDs to gene symbols using bitr
    convert_ensg_to_symbols <- function(genes, organism) {
        # If no genes are provided, return an empty data frame with the expected columns.
        if (length(genes) == 0) return(data.frame(ENSEMBL = character(), SYMBOL = character()))

        tryCatch({
            symbol_df <- bitr(genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = organism)
            if (nrow(symbol_df) == 0) stop("No gene symbols could be mapped.")
            return(symbol_df)
        }, error = function(e) {
            warning("Error converting ENSG IDs to gene symbols using bitr: ", e$message)
            return(data.frame(ENSEMBL = genes, SYMBOL = NA))
        })
    }

    # Convert gene IDs for both upregulated and downregulated genes
    upregulated_symbols   <- convert_ensg_to_symbols(upregulated_genes, organism)
    downregulated_symbols <- convert_ensg_to_symbols(downregulated_genes, organism)

    # Function to perform GO enrichment analysis
    perform_enrichment <- function(symbols_df, organism, ont, pvalue_cutoff) {
        valid_genes <- na.omit(symbols_df$SYMBOL)
        if (length(valid_genes) > 0) {
            enrich_result <- tryCatch({
                enrichGO(
                    gene = valid_genes,
                    OrgDb = get(organism),
                    keyType = "SYMBOL",
                    ont = ont,
                    pvalueCutoff = pvalue_cutoff
                )
            }, error = function(e) {
                warning("Error during GO enrichment for ontology ", ont, ": ", e$message)
                NULL
            })

            if (!is.null(enrich_result) && nrow(as.data.frame(enrich_result)) > 0) {
                enrich_result_df <- as.data.frame(enrich_result)
                # Restore gene symbols in the geneID column
                enrich_result_df$geneID <- sapply(enrich_result_df$geneID, function(x) {
                    ids <- unlist(strsplit(x, "/"))
                    symbols <- symbols_df$SYMBOL[symbols_df$SYMBOL %in% ids]
                    paste(symbols, collapse = "/")
                })
                return(enrich_result_df)
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
    }

    # Define the GO ontologies and corresponding output file paths
    ontologies <- c("BP", "MF", "CC")
    outputs <- list(
        BP = output_bp,
        MF = output_mf,
        CC = output_cc
    )

    # Define an empty data frame template to write out if no enrichment is detected
    empty_df <- data.frame(
        Term = character(),
        GeneRatio = character(),
        BgRatio = character(),
        pvalue = numeric(),
        padj = numeric(),
        Status = character(),
        geneID = character(),
        stringsAsFactors = FALSE
    )

    # Loop over each ontology, perform enrichment, and write the output CSV file
    for (ont in ontologies) {
        up_enrich   <- perform_enrichment(upregulated_symbols, organism, ont, pvalue_cutoff)
        down_enrich <- perform_enrichment(downregulated_symbols, organism, ont, pvalue_cutoff)

        if (!is.null(up_enrich)) {
            up_enrich$Status <- "Upregulated"
        }
        if (!is.null(down_enrich)) {
            down_enrich$Status <- "Downregulated"
        }

        # Combine up and down results; if both are NULL, use the empty data frame
        if (!is.null(up_enrich) && !is.null(down_enrich)) {
            combined <- bind_rows(up_enrich, down_enrich)
        } else if (!is.null(up_enrich)) {
            combined <- up_enrich
        } else if (!is.null(down_enrich)) {
            combined <- down_enrich
        } else {
            combined <- empty_df
        }

        # Write the combined results to the specified output CSV file
        write.csv(combined, file = outputs[[ont]], row.names = FALSE)
    }
}

if (is_snakemake) {
    input_file <- snakemake@input[["de_genes"]]
    output_bp <- snakemake@output[["bp"]]
    output_cc <- snakemake@output[["cc"]]
    output_mf <- snakemake@output[["mf"]]
    organism <- snakemake@params[["organism"]]
    pvalue_cutoff <- snakemake@params[["pvalue_cutoff"]]
} else {
    opt <- parse_args()
    input_file <- opt$i
    output_bp <- opt$o_bp
    output_cc <- opt$o_cc
    output_mf <- opt$o_mf
    organism <- opt$organism
    pvalue_cutoff <- opt$p
}

# We have three file to write and they can be in different folders
dir.create(c(dirname(output_bp)), recursive = TRUE)
dir.create(c(dirname(output_cc)), recursive = TRUE)
dir.create(c(dirname(output_mf)), recursive = TRUE)
go_enrichment_analysis(input_file, output_bp, output_cc, output_mf, organism, pvalue_cutoff)
cat("GO enrichment analysis for all ontologies completed successfully.\n")
