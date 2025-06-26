rule extract_bam_stats:
    input:
        bai="results/aligned/{sample}.bam.bai",
        bam="results/aligned/{sample}.bam"
    output:
        stats="results/qc/samtools/{sample}_bam_stats.txt"
    log:
        "logs/extract_bam_stats_{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_gb = 1,
        runtime_m = 10
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """
