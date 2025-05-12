rule create_rseqc_annotation:
    input:
        gtf=f"../results/genome/{selected_genome}/{selected_genome}_annotation.gtf",
    output:
        bed=f"../results/genome/{selected_genome}/{selected_genome}.bed"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on error

        # Ensure input GTF file exists
        if [ ! -f {input.gtf} ]; then
            echo "Error: GTF file {input.gtf} not found!" >&2
            exit 1
        fi

        # Convert GTF to genePred format
        gtfToGenePred {input.gtf} temp_{selected_genome}.genePred

        # Convert genePred to BED12 format
        genePredToBed temp_{selected_genome}.genePred {output.bed}

        # Clean up intermediate files
        rm -f temp_{selected_genome}.genePred
        """
