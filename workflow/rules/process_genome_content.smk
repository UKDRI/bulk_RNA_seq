rule process_genome_content:
    input:
        bam="../results/aligned/{sample}.bam",
        bai="../results/aligned/{sample}.bam.bai"
    output:
        genome_content="../results/qc/genome_content/{sample}_genome_content.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.genome_content})

        # Extract genome content metrics using samtools idxstats
        samtools idxstats {input.bam} > {output.genome_content}
        """
