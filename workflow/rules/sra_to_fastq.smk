rule download_sra:
    output:
        "data/sra/{sample}.sra"
    conda:
        "../envs/data_retrieval.yaml"
    shell:
        """
        prefetch {wildcards.sample} -O data/sra/
        """

rule sra_to_fastq:
    input:
        "data/sra/{sample}.sra"
    output:
        "data/raw/{sample}_1.fastq.gz",
        "data/raw/{sample}_2.fastq.gz"
    conda:
        "../envs/data_retrieval.yaml"  # Conda environment with parallel-fastq-dump installed
    threads: config["threads"]  # Ensure this is defined in config.yaml
    shell:
        """
        parallel-fastq-dump --sra-id {wildcards.sample} \
            --threads {threads} \
            --outdir data/raw/ \
            --split-files \
            --gzip
        """
