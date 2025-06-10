rule trimming:
    input:
        fastq1=lookup(query="sample_id == '{sample}'", cols="fastq_1", within=samplesheet_df),
        fastq2=lookup(query="sample_id == '{sample}'", cols="fastq_2", within=samplesheet_df),
        adapters="results/adapters/{sample}_adapters.txt",
    output:
        trimmed1="results/trimmed/{sample}_1.atria.fastq.gz",
        trimmed2="results/trimmed/{sample}_2.atria.fastq.gz"
    log:
        "logs/trimming_{sample}.log"
    container:
        config["atria_image"]
    threads: 8
    shell:
        """
        # Read adapter sequences at runtime from the adapter file.
        # Adjust the extraction command based on your file format.
        adapter1=$(awk -F'\\t' 'NR==2 {{print $2}}' {input.adapters})
        adapter2=$(awk -F'\\t' 'NR==2 {{print $7}}' {input.adapters})
        
        echo "Using adapter1: $adapter1" > {log}
        echo "Using adapter2: $adapter2" >> {log}

        # Ensure the output directory exists.
        mkdir -p $(dirname {output.trimmed1})

        # Run Atria.
        atria -r {input.fastq1} -R {input.fastq2} \
            -a "$adapter1" -A "$adapter2" \
            -o $(dirname {output.trimmed1}) --threads {threads} \
            --pcr-dedup --pcr-dedup-count --polyG --polyT --polyA --polyC \
            --poly-length 10 --quality-kmer 5 --length-range 30:500 \
            --max-n 5 --tail-length 12 --check-identifier --force \
            >> {log} 2>&1
        """
