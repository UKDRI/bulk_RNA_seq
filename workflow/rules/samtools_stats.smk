rule samtools_stats:
    input:
        bam="../results/aligned/{sample}.bam"
    output:
        stats="../results/qc/samtools/{sample}_samtools_stats.txt",
        flagstat="../results/qc/samtools/{sample}_flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.stats})
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

