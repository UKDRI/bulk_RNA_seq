rule index_bam:
    input:
        bam="../results/aligned/{sample}.bam"
    output:
        index="../results/aligned/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        """
