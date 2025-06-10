rule combine_expression_matrix:
    input:
        quant_files=expand("results/Quant/Count/quant/{sample}/quant.sf", sample=samples),
    output:
        gene_expression_file="results/Quant/Count/combined_expression.csv"
    log:
        "logs/combine_expression_matrix.log"
    conda:
        "../envs/pandas.yaml"
    priority: 60
    script:
        "../scripts/combine_expression_matrix.py"
