rule convert_gtf_to_bed:
    input:
        gtf="data/genome/{selected_genome}/{selected_genome}_annotation.gtf"
    output:
        bed="data/genome/{selected_genome}/{selected_genome}.bed"
    conda:
        "../envs/gtf2bed_env.yaml"
    shell:
        """
        gtf2bed < {input.gtf} > {output.bed}
        """
