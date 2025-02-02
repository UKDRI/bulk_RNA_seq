rule generate_pca_correlation_plots:
    input:
        combined_expression="../test-dataset/data/Quant/Count/combined_expression.csv",
        samplesheet=config["sample_sheet"]
    output:
        pca="../test-dataset/data/plots/pca_plot.png",
        correlation="../test-dataset/data/plots/correlation_matrix.png"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        Rscript workflow/scripts/generate_pca_correlation.R \
        "{input.combined_expression}" "{input.samplesheet}" \
        "{output.pca}" "{output.correlation}"
        """


rule annotate_genomic_regions:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        gtf="data/genome/mouse/mouse_annotation.gtf"
    output:
        annotated="../test-dataset/data/stats/{sample}_genomic_counts.txt"
    log:
        "../test-dataset/data/stats/{sample}_featureCounts.log"
    conda:
        "../envs/plots.yaml"  # Use a dedicated environment for featureCounts
    shell:
        """
        featureCounts -a "{input.gtf}" -o "{output.annotated}" "{input.bam}" \
        2> "{log}"
        """

rule process_star_output:
    input:
        log="../test-dataset/data/aligned/{sample}_Log.out",
        genomic_counts="../test-dataset/data/stats/{sample}_genomic_counts.txt"
    output:
        mapstat="../test-dataset/data/stats/{sample}_mapstat.tsv",
        read_class="../test-dataset/data/stats/{sample}_read_class.tsv",
        pie_reads="../test-dataset/data/plots/{sample}_read_classification.png"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        python workflow/scripts/process_star_output.py \
        --sample "{wildcards.sample}" \
        --log "{input.log}" \
        --genomic_counts "{input.genomic_counts}" \
        --mapstat "{output.mapstat}" \
        --read_class "{output.read_class}" \
        --pie_reads "{output.pie_reads}"
        """
