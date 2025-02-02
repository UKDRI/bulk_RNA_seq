rule junction_saturation:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        annotation="data/genome/mouse/mouse.bed"
    output:
        pdf="../test-dataset/data/qc_CM/rseqc/{sample}.junctionSaturation_plot.pdf",
        r="../test-dataset/data/qc_CM/rseqc/{sample}.junctionSaturation_plot.r"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run junction_saturation.py
        junction_saturation.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample}
        """