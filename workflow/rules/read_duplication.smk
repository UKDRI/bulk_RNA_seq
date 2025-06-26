rule read_duplication:
    input:
        bai="results/aligned/{sample}.bam.bai",
        bam="results/aligned/{sample}.bam"
    output:
        pdf="results/qc/rseqc/{sample}.DupRate_plot.pdf",
        r="results/qc/rseqc/{sample}.DupRate_plot.r",
        pos_xls="results/qc/rseqc/{sample}.pos.DupRate.xls",
        seq_xls="results/qc/rseqc/{sample}.seq.DupRate.xls"
    log:
        "logs/read_duplication_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 15,
        runtime_m = 30
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run read_duplication.py
        read_duplication.py -i {input.bam} -o $(dirname {output.pdf})/{wildcards.sample} > {log} 2>&1
        """
