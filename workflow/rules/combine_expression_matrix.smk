rule combine_expression_matrix:
    input:
        quant_files=expand("../test-dataset/data/Quant/Count/quant/{sample}/quant.sf", sample=samples),
    output:
        gene_expression_file="../test-dataset/data/Quant/Count/combined_expression.csv"
    conda:
        "../envs/quant.yaml"
    shell:
        """
        python workflow/scripts/combine_expression_matrix.py {input.quant_files} {output.gene_expression_file}
        """
