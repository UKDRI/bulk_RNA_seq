rule read_duplication:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam"
    output:
        pdf="../test-dataset/data/qc/rseqc/{sample}.DupRate_plot.pdf",
        r="../test-dataset/data/qc/rseqc/{sample}.DupRate_plot.r",
        pos_xls="../test-dataset/data/qc/rseqc/{sample}.pos.DupRate.xls",
        seq_xls="../test-dataset/data/qc/rseqc/{sample}.seq.DupRate.xls"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run read_duplication.py
        read_duplication.py -i {input.bam} -o $(dirname {output.pdf})/{wildcards.sample}
        """
