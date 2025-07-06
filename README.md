# FFPE_TDS
A Snakemake pipeline to detect markers on targeted DNA/RNA sequencing data. DNA Markers includes SNV, Indel, SV, MSI. RNA Markers includes Fusion, TPM


## DNA
- fastQC: quality control
- umi_tools: extract UMIs from reads
- bwa-mem: alignment
- sambamba, umi_tools: sorting, indexing, UMI-aware deduplicate
- freebayes: call variants (SNV, Indel etc.)
- mantaSV: call SV
- MSIsensor-pro: MSI
- Quantification (?)


## RNA
- fastQC: quality control
- STAR: alignment
- samtools: sort
- qualimap: quality control
- FusionCatcher/CICERO(slow, but manage to call ITD): gene fusions
- featurecounts: TPM

## Requirements
Snakemake (>=6.0)
Conda for environment management

# Installation
## Linux
1. Clone the repo

2. Install Snakemake
```
# Step 1: Create the environment and install snakemake
conda create -n FFPE_TDS -c bioconda -c conda-forge snakemake

# Step 2: Activate the environment
conda activate FFPE_TDS
```

3. Prepare reference genome hs37d5, create BWA index
```
cd /path/to/
curl -O https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz
bwa index /path/to/hs37d5.fa.gz

```

# Configuration
Edit the config.yaml file to specify your data paths and parameters:
```
# Input data paths
input:
  fastq_dir: "/path/to/your/fastq/files"
  target_bed: "/path/to/your/target/regions.bed"

# Reference genome
reference:
  genome: "/path/to/hs37d5.fa"

# UMI configuration
umi:
  barcode_pattern: "NNNNNNNN"  # Adjust based on your UMI structure
```

# Usage
## Basic Execution
```
# Dry run to check the workflow
snakemake --dry-run

# Run with 8 cores
snakemake --cores 8 --use-conda

# Run with specific samples
snakemake --cores 8 --use-conda results/reports/SAMPLE1_quantification_report.tsv
```

