rule convert_gtf_to_bed:
    input:
        gtf="../results/genome/{selected_genome}/{selected_genome}_annotation.gtf"
    output:
        bed_ops="../results/genome/{selected_genome}/{selected_genome}.bed"
    conda:
        "../envs/gtf2bed_env.yaml"
    shell:
        """
        gtf2bed < {input.gtf} > {output.bed_ops}
        """
