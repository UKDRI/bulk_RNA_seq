import pandas as pd
import os
import argparse


def parse_stats(file_path):
    """
    Parse a samtools stats file and extract relevant metrics.
    """
    data = {}
    sample_name = os.path.basename(file_path).replace('_samtools_stats.txt', '')
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('SN'):
                parts = line.split('\t')
                if len(parts) < 3:
                    continue  # Skip malformed lines
                key = parts[1].strip(':')
                value = parts[2].strip()
                data[key] = value
    data['Sample'] = sample_name
    return data


def merge_stats(input_files, output_file):
    """
    Merge multiple samtools stats files into a single Excel file.
    """
    all_data = []
    for file_path in input_files:
        print(f"Parsing file: {file_path}")
        all_data.append(parse_stats(file_path))

    # Convert to DataFrame and save to Excel
    df = pd.DataFrame(all_data)
    df.to_excel(output_file, index=False)
    print(f"Merged stats saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge samtools stats files into an Excel file.")
    parser.add_argument("--input", nargs='+', required=True, help="List of samtools stats files to merge.")
    parser.add_argument("--output", required=True, help="Output Excel file path.")
    args = parser.parse_args()

    # Merge stats and save
    merge_stats(args.input, args.output)
