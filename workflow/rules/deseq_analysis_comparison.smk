rule deseq2_analysis_multiple:
    input:
        expression="results/Quant/Count/combined_expression.csv",
        metadata=config["sample_sheet"],
        canonicals=expand("resources/{selected_genome}_canonical_transcripts.txt", selected_genome=config["selected_genome"]),
        comparisons_file=config["comparison_sheet"]
    output:
        deg_results=expand("results/Differential/deglist/{comparison}_deg_results.csv", comparison=comparisons)
    log:
        "logs/deseq2_analysis_multiple.log"
    conda:
        "../envs/deseq2_multiple_comparisons.yaml"
    params:
        # Choose organism database based on selected genome.
        species = "human" if config["selected_genome"] == "hg38" else "mouse",
        out_dir = "results/Differential/deglist"
    priority: 50
    script:
        "../scripts/deseq2_multiple_groups.R"
