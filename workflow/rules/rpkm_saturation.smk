rule rpkm_saturation:
    input:
        bam="../results/aligned/{sample}.bam",
        annotation=lambda wildcards: f"../results/genome/{selected_genome}/{selected_genome}.bed"
    output:
        eRPKM="../results/qc/rseqc/{sample}.eRPKM.xls",
        rawCount="../results/qc/rseqc/{sample}.rawCount.xls",
        saturation_pdf="../results/qc/rseqc/{sample}.saturation.pdf",
        saturation_r="../results/qc/rseqc/{sample}.saturation.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on first error

        # Ensure output directory exists
        mkdir -p $(dirname {output.saturation_pdf})

        # Run RPKM_saturation.py
        RPKM_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.saturation_pdf})/{wildcards.sample}

        # Verify outputs exist
        if [ ! -s {output.eRPKM} ] || [ ! -s {output.rawCount} ] || [ ! -s {output.saturation_pdf} ] || [ ! -s {output.saturation_r} ]; then
            echo "Error: RPKM_saturation.py did not generate expected outputs." >&2
            exit 1
        fi
        """
