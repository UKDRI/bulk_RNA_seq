rule build_salmon_index_mouse:
    input:
        sequences="data/genome/mouse/mouse_transcriptome.fa.gz"  # Assuming the file is gzipped
    output:
        index_dir=directory("data/genome/mouse/salmon/salmon_index_mouse")
    log:
        "data/genome/mouse/salmon/salmon_index_mouse/genome_index.log"
    conda:
        "../envs/salmon.yaml"
    threads: 10  # Adjusted threads for stability
    shell:
        """
        # Create the output directory
        mkdir -p {output.index_dir}

        # Run Salmon indexing using gzipped file
        zcat {input.sequences} | salmon index -t /dev/stdin -i {output.index_dir} --threads {threads} > {log} 2>&1

        # Check the exit code and create a placeholder file if indexing fails
        if [ $? -ne 0 ]; then
            echo "Salmon indexing failed. Please check the log for details." >> {log}
            exit 1
        fi
        """
