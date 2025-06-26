rule read_distribution:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai",
        annotation=expand("resources/{selected_genome}.bed", selected_genome=config["selected_genome"])
    output:
        txt=ensure("results/qc/rseqc/{sample}_read_distribution.txt", non_empty=True)
    log:
        "logs/read_distribution_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 2
        runtime_m = 30
    shell:
        """
        # Ensure output directory exists
        mkdir -p $(dirname {output.txt})

        # Run read_distribution.py
        read_distribution.py -i {input.bam} -r {input.annotation} > {output.txt} 2> {log}
        """
