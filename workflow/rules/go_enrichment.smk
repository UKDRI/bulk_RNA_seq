rule go_enrichment:
    input:
        # DESeq2 results file exists, but might have few or no significant genes.
        de_genes = "../results/Differential/deglist/{comparison}_deg_results.csv"
    output:
        bp = "../results/enrichment/GO/Biological_Process/{comparison}_biological_process_enrichment.csv",
        mf = "../results/enrichment/GO/Molecular_Function/{comparison}_molecular_function_enrichment.csv",
        cc = "../results/enrichment/GO/Cellular_Component/{comparison}_cellular_component_enrichment.csv"
    params:
        # Choose organism database based on selected genome.
        organism = "org.Hs.eg.db" if selected_genome == "hg38" else "org.Mm.eg.db",
        pvalue_cutoff = 0.05
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        set -e
        mkdir -p $(dirname {output.bp})
        Rscript --vanilla workflow/scripts/go_enrichment_analysis.R \
            -i {input.de_genes} \
            -o_bp {output.bp} \
            -o_mf {output.mf} \
            -o_cc {output.cc} \
            -org {params.organism} \
            -p {params.pvalue_cutoff}
        """
