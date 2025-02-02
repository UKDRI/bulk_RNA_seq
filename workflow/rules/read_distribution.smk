rule read_distribution:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        annotation="data/genome/mouse/mouse.bed"
    output:
        txt="../test-dataset/data/qc_CM/rseqc/{sample}_read_distribution.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.annotation} > {output.txt}
        """
