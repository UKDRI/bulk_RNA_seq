"""Script to merge samtools stat files"""

import logging
import sys
from typing import List

import pandas as pd

if snakemake is not None:
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(msg)s",
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(msg)s",
    )

logger = logging.getLogger()


def combine_expression_matrix(quant_files: List[str], output_file: str):
    """Function to combine the expression matrix

    Args:
        quant_files: List of path to the quantification files
        output_file: Path to the outpur file in CSV format
    """
    combined_data = []

    # Iterate through each quant file and extract data
    logger.info(quant_files)
    for quant_file in quant_files:
        sample_name = quant_file.split("/")[
            -2
        ]  # Assuming the folder name is the sample name
        try:
            # Read the first file to infer column names
            df = pd.read_csv(quant_file, sep="\t")
            # Identify the identifier and TPM columns
            id_col = "Name" if "Name" in df.columns else df.columns[0]
            tpm_col = "TPM" if "TPM" in df.columns else "TPM"
            df = df[[id_col, tpm_col]].rename(
                columns={id_col: "gene_id", tpm_col: sample_name}
            )
        except Exception as error:
            logger.warning("Error reading %s: %s", quant_file, error)
            continue

        # Store data
        combined_data.append(df.set_index("gene_id"))

    # Check if combined data is populated
    if not combined_data:
        logger.error("Error: No valid quantification files were loaded.")
        sys.exit(1)

    # Concatenate all sample dataframes into a single dataframe along columns
    combined_df = pd.concat(combined_data, axis=1)
    combined_df.reset_index(inplace=True)  # Flatten index, with id as a column

    # Ensure combined_df has data
    if combined_df.empty:
        logger.error("Error: Combined expression data is empty after concatenation.")
        sys.exit(1)

    # Save combined data
    combined_df.to_csv(output_file, index=False)
    logger.info("Expression matrix saved to %s", output_file)


# Input arguments
if snakemake is not None:
    combine_expression_matrix(snakemake.input, snakemake.output[0])
elif __name__ == "__main__":
    quantification_files = sys.argv[1:-1]  # All quantification files
    output_csv = sys.argv[-1]
    combine_expression_matrix(quantification_files, output_csv)
