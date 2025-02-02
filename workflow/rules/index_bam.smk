rule index_bam:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam"
    output:
        index="../test-dataset/data/aligned/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        """
