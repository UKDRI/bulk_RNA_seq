rule combine_expression_matrix:
    priority: 60
    input:
        quant_files=expand("../results/Quant/Count/quant/{sample}/quant.sf", sample=samples),
    output:
        gene_expression_file="../results/Quant/Count/combined_expression.csv"
    conda:
        "../envs/pandas.yaml"
    priority: 60
    script:
        "../scripts/combine_expression_matrix.py"
