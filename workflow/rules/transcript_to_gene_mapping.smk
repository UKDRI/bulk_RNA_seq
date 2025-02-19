rule transcript_to_gene_mapping:
    input:
        gtf=lambda wildcards: f"../results/genome/{selected_genome}/{selected_genome}_annotation.gtf"
    output:
        map_csv="../results/Quant/transcript_to_gene_mapping.csv"
    conda:
        "../envs/transcript_to_gene_env.yaml"
    shell:
        """
        set -e  # Exit on error

        # Run Rscript to extract transcript-to-gene mapping
        Rscript workflow/scripts/extract_transcript_to_gene.R "{input.gtf}" "{output.map_csv}"

        # Verify output file exists and is not empty
        if [ ! -s {output.map_csv} ]; then
            echo "Error: extract_transcript_to_gene.R did not generate expected output." >&2
            exit 1
        fi
        """
