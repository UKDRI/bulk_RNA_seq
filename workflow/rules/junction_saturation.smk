rule junction_saturation:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        annotation=lambda wildcards: f"data/genome/{selected_genome}/{selected_genome}.bed"
    output:
        pdf="../test-dataset/data/qc/rseqc/{sample}.junctionSaturation_plot.pdf",
        r="../test-dataset/data/qc/rseqc/{sample}.junctionSaturation_plot.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on first error

        # Ensure output directory exists
        mkdir -p $(dirname {output.pdf})

        # Run junction_saturation.py
        junction_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample}

        # Move generated files to expected names
        mv $(dirname {output.pdf})/{wildcards.sample}.junctionSaturation_plot.pdf {output.pdf}
        mv $(dirname {output.pdf})/{wildcards.sample}.junctionSaturation_plot.r {output.r}

        # Verify outputs exist
        if [ ! -s {output.pdf} ] || [ ! -s {output.r} ]; then
            echo "Error: junction_saturation.py did not generate expected outputs." >&2
            exit 1
        fi
        """
