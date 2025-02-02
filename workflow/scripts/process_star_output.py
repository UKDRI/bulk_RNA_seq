import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_log(log_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()
    # Extract relevant metrics (customize based on STAR log format)
    stats = {"Unique Reads": 0, "Multi-mapping Reads": 0}  # Example placeholders
    for line in lines:
        if "UNIQUE READS" in line:
            stats["Unique Reads"] = int(line.split()[-1])
        elif "MULTI-MAPPING READS" in line:
            stats["Multi-mapping Reads"] = int(line.split()[-1])
    return stats

def generate_pie_chart(data, output_file):
    labels = data.keys()
    sizes = data.values()
    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title("Read Classification")
    plt.savefig(output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--mapstat", required=True)
    parser.add_argument("--read_class", required=True)
    parser.add_argument("--pie_reads", required=True)
    args = parser.parse_args()

    stats = parse_log(args.log)
    pd.DataFrame([stats]).to_csv(args.mapstat, sep="\t", index=False)
    pd.DataFrame(stats, index=[0]).to_csv(args.read_class, sep="\t", index=False)
    generate_pie_chart(stats, args.pie_reads)

if __name__ == "__main__":
    main()
