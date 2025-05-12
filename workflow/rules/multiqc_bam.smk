rule run_multiqc:
    priority: 30
    input:
        # RSeQC outputs
        infer_experiment=expand("../results/qc/rseqc/{sample}_infer_experiment.txt", sample=samples),
        #gene_body_coverage=expand(".../results/qc/rseqc/{sample}.geneBodyCoverage.txt", sample=samples),
        read_distribution=expand("../results/qc/rseqc/{sample}_read_distribution.txt", sample=samples),
        junction_saturation=expand("../results/qc/rseqc/{sample}.junctionSaturation_plot.pdf", sample=samples),
        dup_rate=expand("../results/qc/rseqc/{sample}.DupRate_plot.pdf", sample=samples),
        saturation_pdf=expand("../results/qc/rseqc/{sample}.saturation.pdf", sample=samples),
        bam_stat=expand("../results/qc/rseqc/{sample}_bam_stat.txt", sample=samples),
        # Samtools outputs
        samtools_stats=expand("../results/qc/samtools/{sample}_samtools_stats.txt", sample=samples),
        samtools_flagstat=expand("../results/qc/samtools/{sample}_flagstat.txt", sample=samples)
    output:
        html="../results/multiqc/bam/multiqc_report.html"
    conda:
        "../envs/multiqc_new.yaml"
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.html})
        
        # Run MultiQC and store the report in the output directory
        multiqc -o $(dirname {output.html}) ../results/qc/
        """
