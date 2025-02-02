# RNA-seq Alignment Pipeline

This repository contains a Snakemake-based RNA-seq alignment pipeline using GitHub Actions for CI/CD.

## Requirements
- Conda
- Snakemake
- GitHub account for CI/CD

## Setting up locally

1. Set up the conda environment:

```bash
conda env create -f environment.yaml
conda activate rna-seq-env
```
   
```bash

STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir data/genome/star_index \
--genomeFastaFiles data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile data/genome/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 100
```

```bash
snakemake --use-conda --cores 4
```

