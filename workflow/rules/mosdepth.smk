rule mosdepth:
    input:
        bam="results/aligned/{sample}.bam",
        index="results/aligned/{sample}.bam.bai"
    output:
        summary="results/qc/mosdepth/{sample}.mosdepth.summary.txt"
    log:
        "logs/mosdepth_{sample}.log"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        # Create the output directory
        mkdir -p $(dirname {output.summary})

        # Run mosdepth with specified parameters
        mosdepth --threads 10 --fast-mode $(dirname {output.summary})/{wildcards.sample} {input.bam} > {log}
        """
