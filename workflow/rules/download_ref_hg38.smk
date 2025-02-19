##################################
# Download + Unzip + Index (hg38)
##################################

# Download the primary assembly (DNA) FASTA file
rule download_human_genome:
    output:
        "../results/genome/hg38/hg38_genome.fa.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
            -O {output}
        """

# Download the cDNA FASTA file (transcriptome)
rule download_human_transcriptome:
    output:
        "../results/genome/hg38/hg38_transcriptome.fa.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
            -O {output}
        """

# Download the GTF annotation file
rule download_human_annotation:
    output:
        "../results/genome/hg38/hg38_annotation.gtf.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz \
            -O {output}
        """

# Unzip the genome for STAR
rule unzip_human_genome:
    input:
        "../results/genome/hg38/hg38_genome.fa.gz"
    output:
        "../results/genome/hg38/hg38_genome.fa"
    shell:
        "gunzip -c {input} > {output}"

# Unzip the transcriptome for Salmon
rule unzip_human_transcriptome:
    input:
        "../results/genome/hg38/hg38_transcriptome.fa.gz"
    output:
        "../results/genome/hg38/hg38_transcriptome.fa"
    shell:
        "gunzip -c {input} > {output}"

# Unzip the annotation for STAR
rule unzip_human_annotation:
    input:
        "../results/genome/hg38/hg38_annotation.gtf.gz"
    output:
        "../results/genome/hg38/hg38_annotation.gtf"
    shell:
        "gunzip -c {input} > {output}"

# (Optional) Build STAR index for alignment
rule build_star_index_genome:
    input:
        genome_fa="../results/genome/hg38/hg38_genome.fa",
        gtf_file="../results/genome/hg38/hg38_annotation.gtf"
    output:
        directory("../results/genome/hg38/star_index_hg38")
    conda:
        "../envs/align.yaml"
    threads: 20
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf_file} \
             --runThreadN {threads} \
             --limitGenomeGenerateRAM 150000000000
        """

# Build Salmon index from the unzipped transcriptome
rule build_salmon_index:
    input:
        transcriptome="../results/genome/hg38/hg38_transcriptome.fa"
    output:
        directory("../results/genome/hg38/salmon/salmon_index_hg38")
    conda:
        "../envs/quant.yaml"
    threads: 8
    shell:
        """
        mkdir -p {output}
        salmon index \
            -t {input.transcriptome} \
            -i {output} \
            -p {threads}
        """
