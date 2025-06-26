rule junction_saturation:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai",
        annotation=expand("resources/{selected_genome}.bed", selected_genome=config["selected_genome"])
    output:
        pdf=ensure("results/qc/rseqc/{sample}.junctionSaturation_plot.pdf", non_empty=True),
        r=ensure("results/qc/rseqc/{sample}.junctionSaturation_plot.r", non_empty=True)
    log:
        "logs/junction_saturation_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 2,
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output.pdf})

        # Run junction_saturation.py
        junction_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample} > {log} 2>&1
        """
