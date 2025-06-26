rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_1.atria.fastq.gz",
        r2="results/trimmed/{sample}_2.atria.fastq.gz"
    output:
        r1_fastqc_html="results/fastqc/trimmed/{sample}_1.atria_fastqc.html",
        r1_fastqc_zip="results/fastqc/trimmed/{sample}_1.atria_fastqc.zip",
        r2_fastqc_html="results/fastqc/trimmed/{sample}_2.atria_fastqc.html",
        r2_fastqc_zip="results/fastqc/trimmed/{sample}_2.atria_fastqc.zip"
    log:
        "logs/fastqc_trimmed_{sample}.log"
    conda:
        "../envs/qc.yaml"
    threads: 2
    resources:
        mem_gb = 1,
        runtime_m = 10
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.r1_fastqc_html})

        # Run FastQC for both paired-end reads
        fastqc -t {threads} \
               -o $(dirname {output.r1_fastqc_html}) \
               {input.r1} {input.r2} > {log} 2>&1
        """
