rule run_multiqc:
    input:
        # RSeQC outputs
        infer_experiment=expand("../test-dataset/data/qc/rseqc/{sample}_infer_experiment.txt", sample=samples),
        #gene_body_coverage=expand("../test-dataset/data/qc/rseqc/{sample}.geneBodyCoverage.txt", sample=samples),
        read_distribution=expand("../test-dataset/data/qc/rseqc/{sample}_read_distribution.txt", sample=samples),
        junction_saturation=expand("../test-dataset/data/qc/rseqc/{sample}.junctionSaturation_plot.pdf", sample=samples),
        dup_rate=expand("../test-dataset/data/qc/rseqc/{sample}.DupRate_plot.pdf", sample=samples),
        saturation_pdf=expand("../test-dataset/data/qc/rseqc/{sample}.saturation.pdf", sample=samples),
        bam_stat=expand("../test-dataset/data/qc/rseqc/{sample}_bam_stat.txt", sample=samples),
        # Samtools outputs
        samtools_stats=expand("../test-dataset/data/qc/samtools/{sample}_samtools_stats.txt", sample=samples),
        samtools_flagstat=expand("../test-dataset/data/qc/samtools/{sample}_flagstat.txt", sample=samples)
    output:
        html="../test-dataset/data/multiqc/bam/multiqc_report.html"
    conda:
        "../envs/multiqc_new.yaml"
    shell:
        """
        # Ensure the output directory exists
        mkdir -p $(dirname {output.html})
        
        # Run MultiQC and store the report in the output directory
        multiqc -o $(dirname {output.html}) ../test-dataset/data/qc/
        """
