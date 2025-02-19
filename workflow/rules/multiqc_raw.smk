rule multiqc_raw:
    input:
        fq1=expand("../results/fastqc/raw/{sample}_1_fastqc.zip", sample=samples),
        fq2=expand("../results/fastqc/raw/{sample}_2_fastqc.zip", sample=samples)
    output:
        html="../results/multiqc/raw/multiqc_report.html"
    conda:
        "../envs/multiqc_new.yaml"
    shell:
        """
        multiqc ../results/fastqc/raw/ -o ../results/multiqc/raw/
        """
