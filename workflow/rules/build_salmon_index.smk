rule build_salmon_index:
    input:
        transcriptome=f"data/genome/{selected_genome}/{selected_genome}_transcriptome.fa"
    output:
        index_dir=directory(f"data/genome/{selected_genome}/salmon/salmon_index_{selected_genome}")
    log:
        f"data/genome/{selected_genome}/salmon/salmon_index_{selected_genome}/genome_index.log"
    conda:
        "../envs/salmon.yaml"
    threads: 20
    shell:
        """
        set -e  # Exit on first error

        # Ensure index directory exists
        mkdir -p {output.index_dir}

        # Run Salmon indexing
        salmon index -t {input.transcriptome} -i {output.index_dir} -p {threads} > {log} 2>&1

        # Check if indexing was successful
        if [ $? -ne 0 ]; then
            echo "Salmon indexing failed. Please check the log for details." >> {log}
            exit 1
        fi
        """
