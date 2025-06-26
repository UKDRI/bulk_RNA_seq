rule salmon_quant_reads:
    input:
        r1="results/trimmed/{sample}_1.atria.fastq.gz",
        r2="results/trimmed/{sample}_2.atria.fastq.gz",
        index=expand(config["salmon_index"], selected_genome=config["selected_genome"])
    output:
        quant="results/Quant/Count/quant/{sample}/quant.sf"
    log:
        "logs/salmon_quant_reads_{sample}.log"
    conda:
        "../envs/quant.yaml"
    params:
        libtype=config["salmon_libtype"],
        extra=config["salmon_extra_args"]
    threads: 50
    resources:
        mem_gb = 1,
        runtime_h = 1
    shell:
        """
        mkdir -p $(dirname {output.quant})
        salmon quant -i {input.index} \
                     -l {params.libtype} \
                     -1 {input.r1} \
                     -2 {input.r2} \
                     -p {threads} \
                     -o $(dirname {output.quant}) \
                     {params.extra} > {log} 2>&1
        """
