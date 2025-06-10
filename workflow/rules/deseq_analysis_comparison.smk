rule deseq2_analysis_multiple:
    priority: 50
    input:
        expression="../results/Quant/Count/combined_expression.csv",
        metadata="data/samplesheet_full.csv",
        canonicals="../results/genome/hg38/hg38_canonical_transcripts.txt",
        comparisons="data/comparison.csv"
    output:
        deg_results=expand("../results/Differential/deglist/{comparison}_deg_results.csv", comparison=comparisons)
    conda:
        "../envs/deseq2_multiple_comparisons.yaml"
    params:
        # Choose organism database based on selected genome.
        species = "human" if config["selected_genome"] == "hg38" else "mouse",
        out_dir = "results/Differential/deglist"
    priority: 50
    script:
        "../scripts/deseq2_multiple_groups.R"
