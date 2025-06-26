rule go_enrichment:
    input:
        # DESeq2 results file exists, but might have few or no significant genes.
        de_genes = "results/Differential/deglist/{comparison}_deg_results.csv"
    output:
        bp = "results/enrichment/GO/Biological_Process/{comparison}_biological_process_enrichment.csv",
        mf = "results/enrichment/GO/Molecular_Function/{comparison}_molecular_function_enrichment.csv",
        cc = "results/enrichment/GO/Cellular_Component/{comparison}_cellular_component_enrichment.csv"
    log:
        "logs/go_enrichment_{comparison}.log"
    conda:
        "../envs/go_enrichment.yaml"
    params:
        # Choose organism database based on selected genome.
        organism = "org.Hs.eg.db" if config["selected_genome"] == "GRCh38" else "org.Mm.eg.db",
        pvalue_cutoff = config["go_enrichment_pvalue"]
    script:
        "../scripts/go_enrichment_analysis.R"
