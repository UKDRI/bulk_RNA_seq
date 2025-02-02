rule salmon_quant_reads:
    input:
        r1="../test-dataset/data/trimmed/{sample}_1.atria.fq.gz",  # Trimmed first paired-end read
        r2="../test-dataset/data/trimmed/{sample}_2.atria.fq.gz",  # Trimmed second paired-end read
        index="data/genome/mouse/salmon/salmon_index_mouse"  # Salmon index
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
        # Ensure output directory exists
        mkdir -p $(dirname {output.quant})

        # Run Salmon quantification
        salmon quant -i {input.index} \
                     -l {params.libtype} \
                     -1 {input.r1} \
                     -2 {input.r2} \
                     -p {threads} \
                     -o $(dirname {output.quant}) \
                     {params.extra} > {log} 2>&1
        """
