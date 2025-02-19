rule extract_bam_stats:
    input:
        bam="../results/aligned/{sample}.bam"
    output:
        stats="../results/qc/samtools/{sample}_bam_stats.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats}
        """
