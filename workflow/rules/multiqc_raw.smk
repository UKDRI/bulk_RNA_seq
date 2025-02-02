rule multiqc_raw:
    input:
        fq1=expand("../test-dataset/data/fastqc/raw/{sample}_1_fastqc.zip", sample=samples),
        fq2=expand("../test-dataset/data/fastqc/raw/{sample}_2_fastqc.zip", sample=samples)
    output:
        html="../test-dataset/data/multiqc/raw/multiqc_report.html"
    conda:
        "../envs/multiqc_new.yaml"
    shell:
        """
        multiqc ../test-dataset/data/fastqc/raw/ -o ../test-dataset/data/multiqc/raw/
        """
