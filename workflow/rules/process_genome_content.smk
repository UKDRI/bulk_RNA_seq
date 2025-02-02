rule process_genome_content:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        bai="../test-dataset/data/aligned_CM/{sample}.bam.bai"
    output:
        genome_content="../test-dataset/data/qc_CM/genome_content/{sample}_genome_content.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.genome_content})

        # Extract genome content metrics using samtools idxstats
        samtools idxstats {input.bam} > {output.genome_content}
        """
