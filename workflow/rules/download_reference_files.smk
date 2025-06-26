# Download the primary assembly (DNA) FASTA file 
rule download_genome:
    output:
        expand(config["genome_fasta"], selected_genome=config["selected_genome"])
    log:
        expand("logs/download_genome_{selected_genome}.log", selected_genome=config["selected_genome"])
    conda:
        "../envs/fetch.yaml"
    params:
        release_version = config[config["selected_genome"]]["release_version"],
        species_name = config[config["selected_genome"]]["species_name"],
        species_name_upper = config[config["selected_genome"]]["species_name_upper"],
        genome_version = config["selected_genome"]
    threads: 1
    resources:
        mem_mb = 500
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output})
        wget --quiet ftp://ftp.ensembl.org/pub/release-{params.release_version}/fasta/{params.species_name}/dna/{params.species_name_upper}.{params.genome_version}.dna.primary_assembly.fa.gz \
            -O {output}.gz > {log} 2>&1
        gunzip {output}.gz 2>> {log}
        """

# Download the GTF annotation file
rule download_annotation:
    output:
        expand(config["genes_gtf"], selected_genome=config["selected_genome"])
    log:
        expand("logs/download_annotation_{selected_genome}.log", selected_genome=config["selected_genome"])
    conda:
        "../envs/fetch.yaml"
    params:
        release_version = config[config["selected_genome"]]["release_version"],
        species_name = config[config["selected_genome"]]["species_name"],
        species_name_upper = config[config["selected_genome"]]["species_name_upper"],
        genome_version = config["selected_genome"]
    threads: 1
    resources:
        mem_mb = 500
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output})
        wget --quiet ftp://ftp.ensembl.org/pub/release-{params.release_version}/gtf/{params.species_name}/{params.species_name_upper}.{params.genome_version}.{params.release_version}.gtf.gz \
            -O {output}.gz > {log} 2>&1
        gunzip {output}.gz 2>> {log}
        """

# Download the cDNA FASTA file (transcriptome)
rule download_transcriptome:
    output:
        expand(config["cdna_fasta"], selected_genome=config["selected_genome"])
    log:
        expand("logs/download_transcriptome_{selected_genome}.log", selected_genome=config["selected_genome"])
    conda:
        "../envs/fetch.yaml"
    params:
        release_version = config[config["selected_genome"]]["release_version"],
        species_name = config[config["selected_genome"]]["species_name"],
        species_name_upper = config[config["selected_genome"]]["species_name_upper"],
        genome_version = config["selected_genome"]
    threads: 1
    resources:
        mem_mb = 500
        runtime_m = 10
    shell:
        """
        mkdir -p $(dirname {output})
        wget --quiet ftp://ftp.ensembl.org/pub/release-{params.release_version}/fasta/{params.species_name}/cdna/{params.species_name_upper}.{params.genome_version}.cdna.all.fa.gz \
            -O {output}.gz > {log} 2>&1
        gunzip {output}.gz 2>> {log}
        """

rule get_canonical_transcripts:
    input:
        expand(config["genes_gtf"], selected_genome=config["selected_genome"])
    output:
        "resources/{selected_genome}_canonical_transcripts.txt"
    localrule: True
    log:
        "logs/get_canonical_transcripts_{selected_genome}.log"
    conda:
        "../envs/fetch.yaml"
    threads: 1
    resources:
        mem_mb = 500
        runtime_m = 10
    shell:
        """
        grep Ensembl_canonical {input} | cut -d ';' -f3 | cut -d '"' -f2 | sort -u > {output} 2> {log}
        """

# (Optional) Build STAR index for alignment
rule build_star_index:
    input:
        genome_fa=expand(config["genome_fasta"], selected_genome=config["selected_genome"]),
        gtf_file=expand(config["genes_gtf"], selected_genome=config["selected_genome"])
    output:
        directory(expand("resources/star_index_{selected_genome}", selected_genome=config["selected_genome"]))
    log:
        "logs/build_star_index.log"
    conda:
        "../envs/align.yaml"
    threads: 20
    resources:
        mem_gb = 33
        runtime_h = 3
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf_file} \
             --runThreadN {threads} > {log} 2>&1
        """

# Build Salmon index from the unzipped transcriptome
rule build_salmon_index:
    input:
        transcriptome = expand(config["cdna_fasta"], selected_genome=config["selected_genome"])
    output:
        directory(expand(config["salmon_index"], selected_genome=config["selected_genome"]))
    log:
        "logs/build_salmon_index.log"
    conda:
        "../envs/quant.yaml"
    threads: 8
    resources:
        mem_gb = 1
        runtime_m = 10
    shell:
        """
        mkdir -p {output}
        salmon index \
            -t {input.transcriptome} \
            -i {output} \
            -p {threads} > {log} 2>&1
        """
