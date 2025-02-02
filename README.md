# bulk_RNA_seq
Snakemake pipeline processing fastq files to the DESeq2, and GO annotations

## Pipeline Features

- **FastQ Quality Control** (`FastQC`, `MultiQC`)
- **Adapter Trimming** (`Atria`)
- **Read Alignment** (`STAR`)
- **Transcript Quantification** (`Salmon`)
- **Quality Control Metrics** (`RSeQC`, `samtools`)
- **Differential Expression Analysis** (`DESeq2`)
- **GO Enrichment Analysis** (`clusterProfiler`)

## 📥 Installation

### 1️⃣ Install Dependencies

```bash
conda install -c bioconda snakemake
```

### 2️⃣ Clone the Repository
```bash
git clone https://github.com/your_username/Bulk_RNA_seq.git
cd Bulk_RNA_seq
```

## 🚀 Running the Pipeline

### 1️⃣ Test Pipeline Execution (Dry Run)
```bash
snakemake --use-conda -np
```

### 2️⃣ Run the Full Pipeline
```bash
snakemake --use-conda --cores 30
```
💡 *Adjust `--cores` based on your system's available CPUs.*

### 3️⃣ Run a Specific Step
```bash
snakemake --use-conda --cores 30 salmon_quant_reads
```

### 4️⃣ Force Rerun of a Specific Step
```bash
snakemake --use-conda -R build_salmon_index --cores 30
```

### 5️⃣ Rerun If File Modification Time Has Changed (recommended)

If you have input files in correct format, you should skip generating them again by: 

```bash
snakemake --use-conda --rerun-triggers mtime --cores 30
```

