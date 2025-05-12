rule deseq2_analysis_multiple:
    input:
        expression="../results/Quant/Count/combined_expression.csv",
        metadata="data/samplesheet_full.csv",
        canonicals="../results/genome/hg38/hg38_canonical_transcripts.txt",
        comparisons="data/comparison.csv"
    output:
        deg_results=expand("../results/Differential/deglist/{comparison}_deg_results.csv", comparison=comparisons)
    conda:
        "../envs/deseq2_multiple_comparisons.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p ../results/data/Differential/deglist/
        Rscript workflow/scripts/deseq2_multiple_groups.R \
            --expression {input.expression} \
            --metadata {input.metadata} \
            --comparisons {input.comparisons} \
            --canonicals {input.canonicals} \
            --output_dir ../results/Differential/deglist/
        """
