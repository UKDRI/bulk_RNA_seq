rule multiqc_trimmed:
    input:
        fastqc_r1=expand("results/fastqc/trimmed/{sample}_1.atria_fastqc.zip", sample=samples),
        fastqc_r2=expand("results/fastqc/trimmed/{sample}_2.atria_fastqc.zip", sample=samples)
    output:
        "results/multiqc/trimmed/multiqc_report.html"
    log:
        "logs/multiqc_trimmed.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        # Run MultiQC if all files are present
        multiqc results/fastqc/trimmed/ -o $(dirname {output}) > {log} 2>&1
        """