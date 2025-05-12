rule junction_saturation:
    input:
        bam="../results/aligned/{sample}.bam",
        annotation=lambda wildcards: f"../results/genome/{selected_genome}/{selected_genome}.bed"
    output:
        pdf="../results/qc/rseqc/{sample}.junctionSaturation_plot.pdf",
        r="../results/qc/rseqc/{sample}.junctionSaturation_plot.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on first error

        # Ensure output directory exists
        mkdir -p $(dirname {output.pdf})

        # Run junction_saturation.py
        junction_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample}

        # Verify outputs exist
        if [ ! -s {output.pdf} ] || [ ! -s {output.r} ]; then
            echo "Error: junction_saturation.py did not generate expected outputs." >&2
            exit 1
        fi
        """
