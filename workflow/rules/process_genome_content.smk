rule process_genome_content:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai"
    output:
        genome_content="results/qc/genome_content/{sample}_genome_content.txt"
    log:
        "logs/process_genome_content_{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_gb = 1,
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output.genome_content})

        # Extract genome content metrics using samtools idxstats
        samtools idxstats {input.bam} > {output.genome_content} 2> {log}
        """
