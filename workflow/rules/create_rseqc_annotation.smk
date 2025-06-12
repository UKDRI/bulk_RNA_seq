rule create_rseqc_annotation:
    input:
        gtf=expand(config["genes_gtf"], selected_genome=config["selected_genome"])
    output:
        bed="resources/{selected_genome}.bed"
    log:
        "logs/create_rseqc_annotation_{selected_genome}.log"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Convert GTF to genePred format
        gtfToGenePred {input.gtf} temp_{wildcards.selected_genome}.genePred 2> {log}

        # Convert genePred to BED12 format
        genePredToBed temp_{wildcards.selected_genome}.genePred {output.bed} 2>> {log}

        # Clean up intermediate files
        rm -f temp_{wildcards.selected_genome}.genePred
        """
