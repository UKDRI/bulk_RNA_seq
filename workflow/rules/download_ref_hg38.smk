# Rule to download the latest human transcriptome (cDNA) from Ensembl
rule download_human_genome:
    output:
        "data/genome/hg38/hg38_genome.fa.gz"
    shell:
        """
        # Create directory if it doesn't exist
        mkdir -p $(dirname {output})

        # Download the DNA FASTA file for the transcriptome
        wget -N ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O {output}
        """


rule download_human_transcriptome:
    output:
        "data/genome/hg38/hg38_transcriptome.fa.gz"
    shell:
        """
        # Create directory if it doesn't exist
        mkdir -p $(dirname {output})
        
        # Download the cDNA FASTA file for the transcriptome
        wget -N ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O {output}
        """

# Rule to download the human GTF annotation file (GRCh38) from Ensembl
rule download_human_annotation:
    output:
        "data/genome/hg38/hg38_annotation.gtf.gz"
    shell:
        """
        # Create directory if it doesn't exist
        mkdir -p $(dirname {output})
        
        # Download the GTF annotation file
        wget -N ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz -O {output}
        """

# Rule to unzip the human genome FASTA file
rule unzip_human_genome:
    input:
        "data/genome/hg38/hg38_genome.fa.gz"
    output:
        "data/genome/hg38/hg38_genome.fa"
    shell:
        """
        # Unzip the transcriptome FASTA file
        gunzip -c {input} > {output}
        """

# Rule to unzip the human transcriptome FASTA file
rule unzip_human_transcriptome:
    input:
        "data/genome/hg38/hg38_transcriptome.fa.gz"
    output:
        "data/genome/hg38/hg38_transcriptome.fa"
    shell:
        """
        # Unzip the transcriptome FASTA file
        gunzip -c {input} > {output}
        """

# Rule to unzip the human annotation file
rule unzip_human_annotation:
    input:
        "data/genome/hg38/hg38_annotation.gtf.gz"
    output:
        "data/genome/hg38/hg38_annotation.gtf"
    shell:
        """
        # Unzip the GTF annotation file
        gunzip -c {input} > {output}
        """

# Rule to build the STAR index for the human transcriptome
rule build_star_index_genome:
    input:
        genome_fa="data/genome/hg38/hg38_genome.fa",
        gtf_file="data/genome/hg38/hg38_annotation.gtf"
    output:
        directory("data/genome/hg38/star_index_hg38")
    conda:
        "../envs/align.yaml"  # Conda environment for STAR
    threads: 50
    shell:
        """
        # Ensure the output directory exists
        mkdir -p {output}
        
        # Run STAR to build an index using the genome FASTA and GTF annotation
        STAR --runMode genomeGenerate --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf_file} --runThreadN {threads} \
             --limitGenomeGenerateRAM 150000000000
        """
