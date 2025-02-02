rule samtools_stats:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam"
    output:
        stats="../test-dataset/data/qc_CM/samtools/{sample}_samtools_stats.txt",
        flagstat="../test-dataset/data/qc_CM/samtools/{sample}_flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.stats})
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

