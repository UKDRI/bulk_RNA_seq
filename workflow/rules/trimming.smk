import csv
import os

def get_adapters(sample):
    """
    Function to read adapters from a file.
    If the file is missing or improperly formatted, raises an error during execution.
    """
    adapter_file = f"../test-dataset/data/adapters/{sample}_adapters.txt"
    try:
        with open(adapter_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)
            if len(rows) < 1:
                raise ValueError(f"Adapter file {adapter_file} is improperly formatted. "
                                 "Expected at least 1 data row.")
            # Extract adapters from columns 'best_adapter1' and 'best_adapter2'
            adapter1 = rows[0].get('best_adapter1', '').strip()
            adapter2 = rows[0].get('best_adapter2', '').strip()
            if not adapter1 or not adapter2:
                raise ValueError(f"Missing adapter sequences in {adapter_file}.")
        return adapter1, adapter2
    except FileNotFoundError:
        raise FileNotFoundError(f"Adapter file {adapter_file} not found for sample {sample}. "
                                "Ensure detect_adapters rule has been executed.")
    except KeyError as e:
        raise ValueError(f"Expected column {e} missing in {adapter_file}. Ensure the file format is correct.")
    except Exception as e:
        raise RuntimeError(f"Error reading adapter file {adapter_file}: {e}")


rule trimming:
    input:
        fastq1="../test-dataset/data/raw2/{sample}_1.fq.gz",
        fastq2="../test-dataset/data/raw2/{sample}_2.fq.gz",
        adapters="../test-dataset/data/adapters/{sample}_adapters.txt"
    output:
        trimmed1="../test-dataset/data/trimmed/{sample}_1.atria.fq.gz",
        trimmed2="../test-dataset/data/trimmed/{sample}_2.atria.fq.gz"
    threads: 8
    conda:
        "../envs/trim-env.yaml"
    params:
        adapters=lambda wildcards: get_adapters(wildcards.sample)
    shell:
        """
        set -euo pipefail

        adapter1="{params.adapters[0]}"
        adapter2="{params.adapters[1]}"

        echo "Using adapter1: $adapter1"
        echo "Using adapter2: $adapter2"

        # Ensure the output directory exists
        mkdir -p $(dirname {output.trimmed1})

        # Run Atria
        tools/Atria/atria-4.1.1/bin/atria -r {input.fastq1} -R {input.fastq2} \
            -a "$adapter1" -A "$adapter2" \
            -o $(dirname {output.trimmed1}) --threads {threads} \
            --pcr-dedup --pcr-dedup-count --polyG --polyT --polyA --polyC \
            --poly-length 10 --quality-kmer 5 --length-range 30:500 \
            --max-n 5 --tail-length 12 --check-identifier

        echo "Trimming completed successfully for sample {wildcards.sample}."
        """
