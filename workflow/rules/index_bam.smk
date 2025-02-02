rule index_bam:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam"
    output:
        index="../test-dataset/data/aligned_CM/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        """
