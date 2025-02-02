rule rpkm_saturation:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        annotation="data/genome/mouse/mouse.bed"
    output:
        eRPKM="../test-dataset/data/qc_CM/rseqc/{sample}.eRPKM.xls",
        rawCount="../test-dataset/data/qc_CM/rseqc/{sample}.rawCount.xls",
        saturation_pdf="../test-dataset/data/qc_CM/rseqc/{sample}.saturation.pdf",
        saturation_r="../test-dataset/data/qc_CM/rseqc/{sample}.saturation.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.saturation_pdf})

        # Run RPKM_saturation.py
        RPKM_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.saturation_pdf})/{wildcards.sample}
        """