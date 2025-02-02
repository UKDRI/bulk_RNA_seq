rule mosdepth:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        index="../test-dataset/data/aligned/{sample}.bam.bai"
    output:
        summary="../test-dataset/data/qc/mosdepth/{sample}.mosdepth.summary.txt"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        # Create the output directory
        mkdir -p $(dirname {output.summary})

        # Run mosdepth with specified parameters
        mosdepth --threads 10 --fast-mode $(dirname {output.summary})/{wildcards.sample} {input.bam}
        """
