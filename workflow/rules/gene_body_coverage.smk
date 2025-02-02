rule gene_body_coverage:
    input:
        bam="../test-dataset/data/aligned_CM/{sample}.bam",
        bai="../test-dataset/data/aligned_CM/{sample}.bam.bai",
        annotation="data/genome/mouse/mouse.bed"
    output:
        pdf="../test-dataset/data/qc_CM/rseqc/{sample}.geneBodyCoverage.curves.pdf",
        txt="../test-dataset/data/qc_CM/rseqc/{sample}.geneBodyCoverage.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run geneBody_coverage.py
        geneBody_coverage.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample}
        geneBody_coverage.py -i {input.bam} -r {input.annotation} -o $(dirname {output.txt})/{wildcards.sample}

        """
