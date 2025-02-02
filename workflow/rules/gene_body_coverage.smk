rule gene_body_coverage:
    input:
        bam="../test-dataset/data/aligned/{sample}.bam",
        bai="../test-dataset/data/aligned/{sample}.bam.bai",
        annotation=lambda wildcards: f"data/genome/{selected_genome}/{selected_genome}.bed"
    output:
        pdf="../test-dataset/data/qc/rseqc/{sample}.geneBodyCoverage.curves.pdf",
        txt="../test-dataset/data/qc/rseqc/{sample}.geneBodyCoverage.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit immediately if a command fails
        
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.pdf})

        # Run geneBody_coverage.py once and generate both PDF and TXT outputs
        geneBody_coverage.py -i {input.bam} -r {input.annotation} -o $(dirname {output.pdf})/{wildcards.sample}

        # Ensure the output file is correctly named
        mv $(dirname {output.pdf})/{wildcards.sample}.geneBodyCoverage.txt {output.txt}
        mv $(dirname {output.pdf})/{wildcards.sample}.geneBodyCoverage.curves.pdf {output.pdf}
        """
