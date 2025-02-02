import os

# Snakemake input and output variables
r1 = snakemake.input.r1
r2 = snakemake.input.r2
adapter_summary = snakemake.input.adapter_summary
r1_trimmed = snakemake.output.r1_trimmed
r2_trimmed = snakemake.output.r2_trimmed
log_file = snakemake.output.log_file

# Parse the detected adapters from the summary file
adapter_r1 = None
adapter_r2 = None
with open(adapter_summary, 'r') as summary_file:
    headers = summary_file.readline()  # Read the headers
    for line in summary_file:
        parts = line.split()
        if len(parts) > 5:  # Ensure we have enough parts in the line
            # Extract adapter for R1 (from 'best_adapter1' column) and R2 (from 'best_adapter2' column)
            adapter_r1 = parts[1] if parts[1] != '0.0' else None  # Check for adapter in best_adapter1
            adapter_r2 = parts[5] if parts[5] != '0.0' else None  # Check for adapter in best_adapter2

# Construct the Atria command with adapters if detected
atria_cmd = f"tools/Atria/atria-4.0.3/bin/atria -r {r1} -R {r2} -o data/trimmed --threads 4"

# Add adapters only if they are detected
if adapter_r1:
    print(f"Detected adapter for R1: {adapter_r1}")
    atria_cmd += f" -a {adapter_r1}"
else:
    print("No adapter detected for R1")

if adapter_r2:
    print(f"Detected adapter for R2: {adapter_r2}")
    atria_cmd += f" -A {adapter_r2}"
else:
    print("No adapter detected for R2")

# Run Atria to trim the reads using the detected adapters (if found)
os.system(atria_cmd)

# Rename the output logs
log_path = "data/trimmed/SRR7657872_1.atria.log"
if os.path.exists(log_path):
    os.rename(log_path, log_file)
else:
    raise FileNotFoundError(f"Atria log file not found at {log_path}")
