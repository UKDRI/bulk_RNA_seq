rule generate_align_pct:
    input:
        bam_files=expand("../test-dataset/data/aligned/{sample}.bam", sample=samples)
    output:
        align_pct="../test-dataset/data/qc/samtools/align_pct.xlsx"
    conda:
        "../envs/samtools.yaml"
    threads: 15
    shell:
        """
        python workflow/scripts/generate_align_pct.py \
        --bam_files {input.bam_files} \
        --output {output.align_pct} \
        --threads {threads}
        """
