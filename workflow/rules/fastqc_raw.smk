rule fastqc_raw:
    input:
        r1=lookup(query="sample_id == '{sample}'", cols="fastq_1", within=samplesheet_df),
        r2=lookup(query="sample_id == '{sample}'", cols="fastq_2", within=samplesheet_df)
    output:
        r1_fastqc_html="results/fastqc/raw/{sample}_1_fastqc.html",
        r1_fastqc_zip="results/fastqc/raw/{sample}_1_fastqc.zip",
        r2_fastqc_html="results/fastqc/raw/{sample}_2_fastqc.html",
        r2_fastqc_zip="results/fastqc/raw/{sample}_2_fastqc.zip"
    log:
        "logs/fastqc_raw_{sample}.log"
    conda:
        "../envs/qc.yaml"
    threads: 2
    resources:
        mem_gb = 1
        runtime_m = 10
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.r1_fastqc_html})

        # Run FastQC for R1 and R2
        fastqc -t {threads} \
               -o $(dirname {output.r1_fastqc_html}) \
               {input.r1} {input.r2} > {log} 2>&1
        """
