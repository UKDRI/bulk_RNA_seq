rule salmon_quant_reads:
    input:
        r1="../test-dataset/data/trimmed/{sample}_1.fq.gz",
        r2="../test-dataset/data/trimmed/{sample}_2.fq.gz",
        index=rules.build_salmon_index.output.index_dir  # Ensure dependency on the built index
    output:
        quant="../test-dataset/data/Quant/Count/quant/{sample}/quant.sf"
    log:
        "../test-dataset/data/Quant/Count/quant/{sample}/{sample}.log"
    params:
        libtype="A",
        extra=""
    threads: 50
    conda:
        "../envs/quant.yaml"
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
