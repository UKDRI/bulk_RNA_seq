# Download mouse genome FASTA file
rule download_mouse_genome:
    output:
        "data/genome/mouse_genome.fa.gz"
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O {output}
        """

# Download mouse GTF annotation file
rule download_mouse_annotation:
    output:
        "data/genome/mouse_annotation.gtf.gz"
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz -O {output}
        """

rule build_star_index:
    input:
        genome_fa="data/genome/mouse/mouse_genome.fa",  # Genome FASTA
        gtf_file="data/genome/mouse/mouse_annotation.gtf"  # GTF
    output:
        directory("data/genome/mouse/star_index_genome")
    conda:
        "../envs/align.yaml"  # Conda environment for STAR
    shell:
        """
        # Run STAR genome generation with provided genome FASTA
        STAR --runMode genomeGenerate --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf_file} --runThreadN 8
        """


# NEW: Download mouse transcriptome for Salmon quantification
rule download_mouse_transcriptome:
    output:
        "data/genome/mouse/mouse_transcriptome.fa.gz"
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz -O {output}
        """

# NEW: Build Salmon index from the transcriptome
rule build_salmon_index:
    input:
        transcriptome="data/genome/mouse/mouse_transcriptome.fa.gz"  # Gzipped transcriptome
    output:
        directory("data/genome/mouse/salmon_index")
    conda:
        "../envs/quant.yaml"  # Conda environment for Salmon
    shell:
        """
        # Run Salmon index generation with gzipped transcriptome file
        salmon index -t {input.transcriptome} -i {output}
        """
