rule deseq2_analysis_multiple:
    input:
        expression="../test-dataset/data/Quant/Count/combined_expression.csv",
        metadata="test-data/samplesheet_full.csv",
        comparisons="test-data/comparison.csv"
    output:
        deg_results=expand("../test-dataset/data/Differential/deglist/{comparison}_deg_results.csv", comparison=comparisons)
    conda:
        "../envs/deseq2_multiple_comparisons.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p ../test-dataset/data/Differential/deglist/
        Rscript workflow/scripts/deseq2_multiple_groups.R \
            --expression {input.expression} \
            --metadata {input.metadata} \
            --comparisons {input.comparisons} \
            --output_dir ../test-dataset/data/Differential/deglist/
        """
