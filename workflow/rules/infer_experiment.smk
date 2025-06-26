rule infer_strand:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai",
        annotation=expand("resources/{selected_genome}.bed", selected_genome=config["selected_genome"])
    output:
        txt=ensure("results/qc/rseqc/{sample}_infer_experiment.txt", non_empty=True)
    log:
        "logs/infer_strand_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 1,
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output.txt})

        # Run strand inference
        infer_experiment.py -i {input.bam} -r {input.annotation} > {output.txt} 2> {log}
        """
