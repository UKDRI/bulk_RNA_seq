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
    input:
        r1="../test-dataset/data/trimmed/{sample}_1.atria.fq.gz",
        r2="../test-dataset/data/trimmed/{sample}_2.atria.fq.gz",
        genome_index="data/genome/mouse/star_index_genome"
    output:
        bam="../test-dataset/data/aligned/{sample}.bam",
        log="../test-dataset/data/aligned/{sample}_Log.out",
        sj="../test-dataset/data/aligned/{sample}_SJ.out.tab"
    params:
        tmp_dir=lambda wildcards: get_tmp_dir(wildcards.sample)
    conda:
        "../envs/align.yaml"
    threads: THREADS
    shell:
        """
        set -euo pipefail

        # Define directories
        output_dir=$(dirname "{output.bam}")
        tmp_dir="{params.tmp_dir}"

        # Ensure output and temp directories exist
        mkdir -p "$output_dir" "$tmp_dir"

        # Check and clean up if temp directory already exists
        if [[ -d "$tmp_dir" ]]; then
            echo "Temporary directory already exists. Cleaning up: $tmp_dir"
            rm -rf "$tmp_dir"
        fi

        # Create fresh temp directory
        mkdir -p "$tmp_dir"

        # Clean up temp directory on exit
        trap "rm -rf $tmp_dir" EXIT

        # Run STAR alignment
        echo "Running STAR alignment for sample: {wildcards.sample}"
        STAR --genomeDir "{input.genome_index}" \
             --readFilesIn "{input.r1}" "{input.r2}" \
             --readFilesCommand zcat \
             --runThreadN {threads} \
             --outFileNamePrefix "$tmp_dir/" \
             --outSAMtype BAM SortedByCoordinate \
             --limitBAMsortRAM {STAR_RAM} \
             --outTmpDir "$tmp_dir/_STARtmp"

        # Check if STAR outputs are generated
        if [[ ! -f "$tmp_dir/Aligned.sortedByCoord.out.bam" || ! -f "$tmp_dir/Log.final.out" || ! -f "$tmp_dir/SJ.out.tab" ]]; then
            echo "Error: STAR alignment outputs are missing." >&2
            exit 1
        fi

        # Move STAR outputs
        mv "$tmp_dir/Aligned.sortedByCoord.out.bam" "{output.bam}"
        mv "$tmp_dir/Log.final.out" "{output.log}"
        mv "$tmp_dir/SJ.out.tab" "{output.sj}"

        echo "STAR alignment completed for sample: {wildcards.sample}"
        """

