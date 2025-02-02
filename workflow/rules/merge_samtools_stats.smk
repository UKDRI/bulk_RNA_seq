rule merge_samtools_stats:
    input:
        stats_files=expand("../test-dataset/data/qc/samtools/{sample}_samtools_stats.txt", sample=samples)
    output:
        merged_stats="../test-dataset/data/qc/samtools/merged_samtools_stats.xlsx"
    conda:
        "../envs/merge_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.merged_stats})
        python workflow/scripts/merge_samtools_stats.py --input {input.stats_files} --output {output.merged_stats}
        """
