rule rpkm_saturation:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        annotation=lambda wildcards: f"data/genome/{selected_genome}/{selected_genome}.bed"
    output:
        eRPKM="../test-dataset/data/qc/rseqc/{sample}.eRPKM.xls",
        rawCount="../test-dataset/data/qc/rseqc/{sample}.rawCount.xls",
        saturation_pdf="../test-dataset/data/qc/rseqc/{sample}.saturation.pdf",
        saturation_r="../test-dataset/data/qc/rseqc/{sample}.saturation.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on first error

        # Ensure output directory exists
        mkdir -p $(dirname {output.saturation_pdf})

        # Run RPKM_saturation.py
        RPKM_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.saturation_pdf})/{wildcards.sample}

        # Move generated files to expected names
        mv $(dirname {output.saturation_pdf})/{wildcards.sample}.eRPKM.xls {output.eRPKM}
        mv $(dirname {output.saturation_pdf})/{wildcards.sample}.rawCount.xls {output.rawCount}
        mv $(dirname {output.saturation_pdf})/{wildcards.sample}.saturation.pdf {output.saturation_pdf}
        mv $(dirname {output.saturation_pdf})/{wildcards.sample}.saturation.r {output.saturation_r}

        # Verify outputs exist
        if [ ! -s {output.eRPKM} ] || [ ! -s {output.rawCount} ] || [ ! -s {output.saturation_pdf} ] || [ ! -s {output.saturation_r} ]; then
            echo "Error: RPKM_saturation.py did not generate expected outputs." >&2
            exit 1
        fi
        """
