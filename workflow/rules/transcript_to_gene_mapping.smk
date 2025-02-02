rule transcript_to_gene_mapping:
    input:
        gtf="data/genome/mouse/mouse_annotation.gtf"
    output:
        map_csv="../test-dataset/data/Quant/transcript_to_gene_mapping.csv"
    conda:
        "../envs/transcript_to_gene_env.yaml"
    shell:
        """
        Rscript workflow/scripts/extract_transcript_to_gene.R
        """
