rule fastqc_trimmed:
    input:
        r1="../test-dataset/data/trimmed/{sample}_1.atria.fq.gz",
        r2="../test-dataset/data/trimmed/{sample}_2.atria.fq.gz"
    output:
        r1_fastqc_html="../test-dataset/data/fastqc/trimmed/{sample}_1.atria_fastqc.html",
        r1_fastqc_zip="../test-dataset/data/fastqc/trimmed/{sample}_1.atria_fastqc.zip",
        r2_fastqc_html="../test-dataset/data/fastqc/trimmed/{sample}_2.atria_fastqc.html",
        r2_fastqc_zip="../test-dataset/data/fastqc/trimmed/{sample}_2.atria_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    threads: 20
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.r1_fastqc_html})

        # Run FastQC for both paired-end reads
        fastqc -t {threads} \
               -o $(dirname {output.r1_fastqc_html}) \
               {input.r1} {input.r2}
        """
