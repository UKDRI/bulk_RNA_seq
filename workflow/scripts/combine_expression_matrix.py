import pandas as pd
import sys
import re

# Input arguments
quant_files = sys.argv[1:-1]  # All quantification files
output_file = sys.argv[-1]

# Initialize a dictionary to store each sample's data
combined_data = {}

# Iterate through each quant file and extract data
for quant_file in quant_files:
    sample_name = quant_file.split("/")[-2]  # Assuming the folder name is the sample name
    try:
        # Read the first file to infer column names
        df = pd.read_csv(quant_file, sep="\t")
        # Identify the identifier and TPM columns
        id_col = "Name" if "Name" in df.columns else df.columns[0]
        tpm_col = "TPM" if "TPM" in df.columns else "TPM"
        df = df[[id_col, tpm_col]].rename(columns={id_col: "gene_id", tpm_col: sample_name})
    except Exception as e:
        print(f"Error reading {quant_file}: {e}")
        continue

    # Store data
    combined_data[sample_name] = df.set_index("gene_id")

# Check if combined data is populated
if not combined_data:
    print("Error: No valid quantification files were loaded.")
    sys.exit(1)

# Concatenate all sample dataframes into a single dataframe along columns
combined_df = pd.concat(combined_data.values(), axis=1)
combined_df.reset_index(inplace=True)  # Flatten index, with id as a column

# Ensure combined_df has data
if combined_df.empty:
    print("Error: Combined expression data is empty after concatenation.")
    sys.exit(1)

# Save combined data
combined_df.to_csv(output_file, index=False)
print(f"Expression matrix saved to {output_file}")
