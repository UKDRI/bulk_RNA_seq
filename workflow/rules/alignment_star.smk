# Global variable definitions
import os
import time

TMP_BASE_DIR = "../test-dataset/data/aligned/tmp/star_jobs"
THREADS = 40
STAR_RAM = 10000000000  # Limit for STAR's RAM usage in bytes

# Ensure TMP_BASE_DIR exists
os.makedirs(TMP_BASE_DIR, exist_ok=True)

# Function to generate a unique temporary directory
def get_tmp_dir(sample):
    return os.path.join(TMP_BASE_DIR, f"STARtmp_{sample}_{int(time.time() * 1000)}")

# Function to check and adjust file descriptor limits
def adjust_ulimit():
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    new_limit = min(hard, 65536)  # Set to 65536 or the maximum allowed
    if soft < new_limit:
        resource.setrlimit(resource.RLIMIT_NOFILE, (new_limit, hard))

# Adjust ulimit
adjust_ulimit()

rule alignment_star:
    """
    Align paired-end reads with STAR for a sample that includes {sample} and {genome}.
    The genome index is computed from wildcards: data/genome/{genome}/star_index_{genome}.
    """
    input:
        # Trimmed FASTQs
        r1 = "../test-dataset/data/trimmed/{sample}_1.atria.fq.gz",
        r2 = "../test-dataset/data/trimmed/{sample}_2.atria.fq.gz",
        # Use a lambda to derive the star index path from the genome wildcard
        star_index=f"data/genome/{selected_genome}/star_index_{selected_genome}"
    output:
        bam = "../test-dataset/data/aligned/{sample}.bam",
        log = "../test-dataset/data/aligned/{sample}_Log.out",
        sj  = "../test-dataset/data/aligned/{sample}_SJ.out.tab"
    params:
        tmp_dir = lambda wildcards: get_tmp_dir(wildcards.sample),
        star_ram = STAR_RAM
    conda:
        "../envs/align.yaml"
    threads: THREADS
    shell:
        r"""
        set -euo pipefail

        OUTPUT_DIR=$(dirname "{output.bam}")
        TMP_DIR="{params.tmp_dir}"

        mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

        # Clean up temporary directory on exit
        trap 'rm -rf "$TMP_DIR"' EXIT

        echo "Running STAR alignment for sample: {wildcards.sample}"

        STAR \
          --genomeDir "{input.star_index}" \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --runThreadN {threads} \
          --outFileNamePrefix "$TMP_DIR/" \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM {params.star_ram} \
          --outTmpDir "$TMP_DIR/_STARtmp"

        # Check that STAR outputs exist and are non-empty
        if [[ ! -s "$TMP_DIR/Aligned.sortedByCoord.out.bam" || \
              ! -s "$TMP_DIR/Log.final.out" || \
              ! -s "$TMP_DIR/SJ.out.tab" ]]; then
          echo "Error: STAR alignment outputs are missing or empty." >&2
          exit 1
        fi

        mv "$TMP_DIR/Aligned.sortedByCoord.out.bam" "{output.bam}"
        mv "$TMP_DIR/Log.final.out" "{output.log}"
        mv "$TMP_DIR/SJ.out.tab" "{output.sj}"

        echo "STAR alignment completed for sample: {wildcards.sample}"
        """
