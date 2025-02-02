# Save this as `workflow/scripts/convert_transcript_to_gene.py`
import pandas as pd
import sys

# Read input arguments
expression_file = sys.argv[1]
mapping_file = sys.argv[2]
output_file = sys.argv[3]

# Load data
expression_df = pd.read_csv(expression_file, index_col=0)  # Load expression data with transcript IDs as index
mapping_df = pd.read_csv(mapping_file)  # Load transcript-to-gene mapping file

# Check if expected columns are present in the mapping file
if 'transcript_id' not in mapping_df.columns or 'gene_id' not in mapping_df.columns:
    raise ValueError("Mapping file must contain 'transcript_id' and 'gene_id' columns")

# Ensure 'transcript_id' is a column in expression_df for merging
expression_df = expression_df.rename_axis("transcript_id").reset_index()

# Merge expression data with transcript-to-gene mapping
merged_df = expression_df.merge(mapping_df, on="transcript_id", how="inner")

# Group by 'gene_id' and aggregate (sum) expression values for each gene
gene_expression_df = merged_df.groupby("gene_id").sum()

# Save output
gene_expression_df.to_csv(output_file)
print(f"Gene expression matrix saved to {output_file}")
