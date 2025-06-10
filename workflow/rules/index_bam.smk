rule index_bam:
    input:
        bam="results/aligned/{sample}.bam"
    output:
        index="results/aligned/{sample}.bam.bai"
    log:
        "logs/index_bam_{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam} 2> {log}
        """
