rule fastqc_raw:
    input:
        r1="../rawdata/{sample}_1.fq.gz",
        r2="../rawdata/{sample}_2.fq.gz"
    output:
        r1_fastqc_html="../results/fastqc/raw/{sample}_1_fastqc.html",
        r1_fastqc_zip="../results/fastqc/raw/{sample}_1_fastqc.zip",
        r2_fastqc_html="../results/fastqc/raw/{sample}_2_fastqc.html",
        r2_fastqc_zip="../results/fastqc/raw/{sample}_2_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    threads: 2
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.r1_fastqc_html})

        # Run FastQC for R1 and R2
        fastqc -t {threads} \
               -o $(dirname {output.r1_fastqc_html}) \
               {input.r1} {input.r2}
        """
