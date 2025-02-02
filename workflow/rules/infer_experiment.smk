rule infer_strand:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        annotation="data/genome/mouse/mouse.bed"
    output:
        txt="../test-dataset/data/qc_CM/rseqc/{sample}_infer_experiment.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        infer_experiment.py -i {input.bam} -r {input.annotation} > {output.txt}
        """
