import pandas as pd
import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm  # For progress bar
import logging


# Set up logging
logging.basicConfig(
    filename="generate_align_pct.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)


def run_command(cmd):
    """
    Run a shell command and return its output.
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"Error executing command: {cmd}\n{result.stderr}")
    return result.stdout.strip()


def extract_metrics(bam_file):
    """
    Extract required metrics from a BAM file using samtools.
    """
    sample = os.path.basename(bam_file).replace(".bam", "")
    logging.info(f"Processing sample: {sample}...")

    try:
        # Combine commands to minimize repeated calls
        stats = run_command(
            f"samtools view {bam_file} | "
            f"awk 'BEGIN {{ OFS=\"\\t\" }} {{ "
            f"if (and($2, 64)) read1_map++; "
            f"if (and($2, 128)) read2_map++; "
            f"if (and($2, 16)) negative_map++; "
            f"else positive_map++; "
            f"if ($6 ~ /N/) splice_map++; "
            f"if (and($2, 2)) proper_map++; "
            f"}} END {{ print read1_map, read2_map, negative_map, positive_map, splice_map, proper_map }}'"
        ).split()

        total_reads = int(run_command(f"samtools view -c {bam_file}"))
        total_map = int(run_command(f"samtools view -c -F 4 {bam_file}"))
        unique_map = int(run_command(f"samtools view -c -q 1 {bam_file}"))
        multi_map = total_map - unique_map

        read1_map, read2_map, negative_map, positive_map, splice_map, proper_map = map(int, stats)
        unsplice_map = total_map - splice_map

        # Return metrics as a dictionary
        return {
            "sample": sample,
            "total_reads": total_reads,
            "total_map": total_map,
            "unique_map": unique_map,
            "multi_map": multi_map,
            "read1_map": read1_map,
            "read2_map": read2_map,
            "positive_map": positive_map,
            "negative_map": negative_map,
            "splice_map": splice_map,
            "unsplice_map": unsplice_map,
            "proper_map": proper_map,
        }
    except Exception as e:
        logging.error(f"Error processing sample {sample}: {e}")
        return None


def generate_align_pct(bam_files, output_file, threads):
    """
    Generate align_pct.xlsx from BAM files.
    """
    data = []

    print("Starting alignment metrics extraction...")
    logging.info("Starting alignment metrics extraction...")

    # Use multiprocessing for parallel file processing with a progress bar
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for metrics in tqdm(executor.map(extract_metrics, bam_files), total=len(bam_files)):
            if metrics:
                data.append(metrics)

    print("Processing complete. Saving results to Excel...")
    logging.info("Processing complete. Saving results to Excel...")

    # Convert to DataFrame
    df = pd.DataFrame(data)

    # Save to Excel
    df.to_excel(output_file, index=False, engine="openpyxl")
    print(f"align_pct saved to {output_file}")
    logging.info(f"align_pct saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate align_pct.xlsx from BAM files")
    parser.add_argument("--bam_files", nargs="+", required=True, help="List of BAM files")
    parser.add_argument("--output", required=True, help="Path to save align_pct.xlsx")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    args = parser.parse_args()

    generate_align_pct(args.bam_files, args.output, args.threads)
