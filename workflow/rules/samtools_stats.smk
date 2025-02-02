rule samtools_stats:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam"
    output:
        stats="../test-dataset/data/qc/samtools/{sample}_samtools_stats.txt",
        flagstat="../test-dataset/data/qc/samtools/{sample}_flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.stats})
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

