rule extract_bam_stats:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam"
    output:
        stats="../test-dataset/data/qc_CM/samtools/{sample}_bam_stats.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats}
        """
