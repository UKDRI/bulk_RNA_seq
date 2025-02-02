rule build_salmon_index_human:
    input:
        sequences="data/genome/human/human_transcriptome.fa"
    output:
        index_dir=directory("data/genome/human/salmon/salmon_index_human")
    log:
        "data/genome/human/salmon/salmon_index_human/genome_index.log"
    conda:
        "../envs/salmon.yaml"
    threads: 10  # Reduced number of threads for stability
    shell:
        """
        # Create the output directory
        mkdir -p {output.index_dir}

        # Run Salmon indexing, redirecting stdout and stderr to the log file
        salmon index -t {input.sequences} -i {output.index_dir} --threads {threads} > {log} 2>&1

        # Check the exit code and create a placeholder file if indexing fails
        if [ $? -ne 0 ]; then
            echo "Salmon indexing failed. Please check the log for details." >> {log}
            exit 1
        fi
        """
