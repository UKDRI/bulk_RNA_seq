rule gene_body_coverage:
    input:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai",
        annotation=expand("resources/{selected_genome}.bed", selected_genome=config["selected_genome"])
    output:
        pdf="results/qc/rseqc/{sample}.geneBodyCoverage.curves.pdf",
        txt="results/qc/rseqc/{sample}.geneBodyCoverage.txt"
    log:
        "logs/gene_body_coverage_{sample}.log"
    conda:
        "../envs/rseqc_env.yaml"
    threads: 1
    resources:
        mem_gb = 2,
        runtime_h = 1
    shell:
        """
        mkdir -p $(dirname {output.pdf})

        # Run geneBody_coverage.py once and generate both PDF and TXT outputs
        geneBody_coverage.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample} > {log} 2>&1
        """
