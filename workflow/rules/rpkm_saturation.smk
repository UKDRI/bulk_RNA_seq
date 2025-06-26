rule rpkm_saturation:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai",
        annotation=expand("resources/{selected_genome}.bed", selected_genome=config["selected_genome"])
    output:
        eRPKM=ensure("results/qc/rseqc/{sample}.eRPKM.xls", non_empty=True),
        rawCount=ensure("results/qc/rseqc/{sample}.rawCount.xls", non_empty=True),
        saturation_pdf=ensure("results/qc/rseqc/{sample}.saturation.pdf", non_empty=True),
        saturation_r=ensure("results/qc/rseqc/{sample}.saturation.r", non_empty=True)
    log:
        "logs/rpkm_saturation_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 34
        runtime_h = 2
    shell:
        """
        # Ensure output directory exists
        mkdir -p $(dirname {output.saturation_pdf})

        # Run RPKM_saturation.py
        RPKM_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.saturation_pdf})/{wildcards.sample} > {log} 2>&1
        """
