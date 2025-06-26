rule samtools_stats:
    input:
        bai="results/aligned/{sample}.bam.bai",
        bam="results/aligned/{sample}.bam"
    output:
        stats="results/qc/samtools/{sample}_samtools_stats.txt",
        flagstat="results/qc/samtools/{sample}_flagstat.txt"
    log:
        "logs/samtools_stats_{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb = 500
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output.stats})
        samtools stats {input.bam} > {output.stats} 2> {log}
        samtools flagstat {input.bam} > {output.flagstat} 2>> {log}
        """

