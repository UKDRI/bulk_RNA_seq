# Rule for running RSeQC bam statistics
rule bam_stat:
    input:
        bai="results/aligned/{sample}.bam.bai",
        bam="results/aligned/{sample}.bam"
    output:
        txt="results/qc/rseqc/{sample}_bam_stat.txt"
    log:
        "logs/bam_stat_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.txt})

        # Run bam_stat.py to generate BAM file statistics
        bam_stat.py -i {input.bam} > {output.txt} 2> {log}
        """
