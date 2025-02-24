# Rule for running RSeQC bam statistics
rule bam_stat:
    input:
        bam="../results/aligned/{sample}.bam"
    output:
        txt="../results/qc/rseqc/{sample}_bam_stat.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.txt})

        # Run bam_stat.py to generate BAM file statistics
        bam_stat.py -i {input.bam} > {output.txt}
        """
