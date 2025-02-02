rule download_test_data:
    conda:
        "../envs/download.yaml"  # Conda environment for downloading
    output:
        "test-data/raw/{sample}_{read}.fastq.gz"
    params:
        url=lambda wildcards: samplesheet_df.loc[
            (samplesheet_df['group'] == wildcards.sample.split('_rep')[0]) &
            (samplesheet_df['replicate'] == int(wildcards.sample.split('_rep')[1]))
        ][f'fastq_{1 if wildcards.read == "1" else 2}'].values[0]
    shell:
        """
        mkdir -p $(dirname {output}) && wget -O {output} {params.url}
        """
