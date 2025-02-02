rule go_enrichment:
    input:
        de_genes="../test-dataset/data/Differential/deglist/{comparison}_deg_results.csv"
    params:
        organism="org.Mm.eg.db",
        pvalue_cutoff=0.05
    output:
        bp="../test-dataset/data/enrichment/GO/Biological_Process/{comparison}_biological_process_enrichment.csv",
        mf="../test-dataset/data/enrichment/GO/Molecular_Function/{comparison}_molecular_function_enrichment.csv",
        cc="../test-dataset/data/enrichment/GO/Cellular_Component/{comparison}_cellular_component_enrichment.csv"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p $(dirname "{output.bp}")
        Rscript --vanilla workflow/scripts/go_enrichment_analysis.R \
            -i "{input.de_genes}" \
            -o_bp "{output.bp}" \
            -o_mf "{output.mf}" \
            -o_cc "{output.cc}" \
            -org "{params.organism}" \
            -p "{params.pvalue_cutoff}"
        """
