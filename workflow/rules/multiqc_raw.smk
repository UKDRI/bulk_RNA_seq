rule multiqc_raw:
    input:
        fq1=expand("results/fastqc/raw/{sample}_1_fastqc.zip", sample=samples),
        fq2=expand("results/fastqc/raw/{sample}_2_fastqc.zip", sample=samples)
    output:
        html="results/multiqc/raw/multiqc_report.html"
    log:
        "logs/multiqc_raw.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.html})
        multiqc results/fastqc/raw/ -o $(dirname {output.html}) > {log} 2>&1
        """
