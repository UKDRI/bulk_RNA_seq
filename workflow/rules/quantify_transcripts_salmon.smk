rule salmon_quant_reads:
    input:
        r1="../results/trimmed/{sample}_1.atria.fq.gz",
        r2="../results/trimmed/{sample}_2.atria.fq.gz",
        index=f"../results/genome/{selected_genome}/salmon/salmon_index_{selected_genome}"
    output:
        quant="../results/Quant/Count/quant/{sample}/quant.sf"
    log:
        "../results/Quant/Count/quant/{sample}/{sample}.log"
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
