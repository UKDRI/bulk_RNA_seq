# Download the primary assembly (DNA) FASTA file (Release 110)
rule download_mouse_genome:
    output:
        "../results/genome/mm39/mm39_genome.fa.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
            -O {output}
        """

# Download the GTF annotation file
rule download_mouse_annotation:
    output:
        "../results/genome/mm39/mm39_annotation.gtf.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz \
            -O {output}
        """

# Download the cDNA FASTA file (transcriptome)
rule download_mouse_transcriptome:
    output:
        "../results/genome/mm39/mm39_transcriptome.fa.gz"
    shell:
        """
        mkdir -p $(dirname {output})
        wget -N ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz \
            -O {output}
        """

# Unzip the genome for STAR
rule unzip_mouse_genome:
    input:
        "../results/genome/mm39/mm39_genome.fa.gz"
    output:
        "../results/genome/mm39/mm39_genome.fa"
    shell:
        "gunzip -c {input} > {output}"

# Unzip the annotation for STAR
rule unzip_mouse_annotation:
    input:
        "../results/genome/mm39/mm39_annotation.gtf.gz"
    output:
        "../results/genome/mm39/mm39_annotation.gtf"
    shell:
        "gunzip -c {input} > {output}"

# Unzip the transcriptome for Salmon
rule unzip_mouse_transcriptome:
    input:
        "../results/genome/mm39/mm39_transcriptome.fa.gz"
    output:
        "../results/genome/mm39/mm39_transcriptome.fa"
    shell:
        "gunzip -c {input} > {output}"

# (Optional) Build STAR index for alignment
rule build_star_index_mm39:
    input:
        genome_fa="../results/genome/mm39/mm39_genome.fa",
        gtf_file="../results/genome/mm39/mm39_annotation.gtf"
    output:
        directory("../results/genome/mm39/star_index_mm39")
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
             --runThreadN {threads}
        """

# Build Salmon index from the unzipped transcriptome
rule build_salmon_index_mm39:
    input:
        transcriptome = "../results/genome/mm39/mm39_transcriptome.fa"
    output:
        directory("../results/genome/mm39/salmon/salmon_index_mm39")
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
