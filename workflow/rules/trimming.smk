rule trimming:
    input:
        fastq1="../rawdata/{sample}_1.fq.gz",
        fastq2="../rawdata/{sample}_2.fq.gz",
        adapters="../results/adapters/{sample}_adapters.txt",
        atria="../results/tools/Atria/atria-4.1.1/bin/atria"
    output:
        trimmed1="../results/trimmed/{sample}_1.atria.fq.gz",
        trimmed2="../results/trimmed/{sample}_2.atria.fq.gz"
    threads: 8
    conda:
        "../envs/trim-env.yaml"
    shell:
        """
        set -euo pipefail

        # Read adapter sequences at runtime from the adapter file.
        # Adjust the extraction command based on your file format.
        adapter1=$(awk -F'\t' 'NR==2 {{print $1}}' {input.adapters})
        adapter2=$(awk -F'\t' 'NR==2 {{print $2}}' {input.adapters})
        
        echo "Using adapter1: $adapter1"
        echo "Using adapter2: $adapter2"

        # Ensure the output directory exists.
        mkdir -p $(dirname {output.trimmed1})

        # Run Atria.
        {input.atria} -r {input.fastq1} -R {input.fastq2} \
            -a "$adapter1" -A "$adapter2" \
            -o $(dirname {output.trimmed1}) --threads {threads} \
            --pcr-dedup --pcr-dedup-count --polyG --polyT --polyA --polyC \
            --poly-length 10 --quality-kmer 5 --length-range 30:500 \
            --max-n 5 --tail-length 12 --check-identifier

        echo "Trimming completed successfully for sample {wildcards.sample}."
        """
