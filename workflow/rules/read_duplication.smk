rule read_duplication:
    input:
        bam="../results/aligned/{sample}.bam"
    output:
        pdf="../results/qc/rseqc/{sample}.DupRate_plot.pdf",
        r="../results/qc/rseqc/{sample}.DupRate_plot.r",
        pos_xls="../results/qc/rseqc/{sample}.pos.DupRate.xls",
        seq_xls="../results/qc/rseqc/{sample}.seq.DupRate.xls"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run read_duplication.py
        read_duplication.py -i {input.bam} -o $(dirname {output.pdf})/{wildcards.sample}
        """
