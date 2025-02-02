rule infer_strand:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        annotation=lambda wildcards: f"data/genome/{selected_genome}/{selected_genome}.bed"
    output:
        txt="../test-dataset/data/qc/rseqc/{sample}_infer_experiment.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit immediately if a command fails

        # Ensure output directory exists
        mkdir -p $(dirname {output.txt})

        # Run strand inference
        infer_experiment.py -i {input.bam} -r {input.annotation} > {output.txt} 2>&1

        # Check if output was created successfully
        if [ ! -s {output.txt} ]; then
            echo "Error: infer_experiment.py did not generate output" >&2
            exit 1
        fi
        """
