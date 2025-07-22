## FFPE_TDS
A Snakemake pipeline to detect markers on targeted DNA/RNA sequencing data. DNA Markers includes SNV, Indel, SV, MSI. RNA Markers includes Fusion, TPM


### DNA
- fastQC: quality control
- umi_tools: extract UMIs from reads  (working on using GATK picard UmiAwareMarkDuplicatesWithMateCigar)
- bwa-mem: alignment
- sambamba, umi_tools: sorting, indexing, UMI-aware deduplicate
- freebayes: call variants (SNV, Indel etc.)
- mantaSV: call SV
- MSIsensor-pro: MSI
- Quantification (?)


### RNA
- fastQC: quality control
- STAR: alignment
- samtools: sort
- qualimap: quality control
- FusionCatcher/CICERO(slow, but manage to call ITD): gene fusions
- featurecounts: TPM

## Requirements
Snakemake
Conda for environment management

##  Installation
### Linux
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

### Quick Start

1. **Configure your data:**
Edit config.yaml to set paths to your FASTQ files, reference genome, and BED file.

2. **Setup environment:**
```
./setup.sh
```

3. **Run the pipeline:**
```
./run_pipeline.sh -c 8  # Use 8 cores

snakemake --cores 18 --sdm conda # Use 18 cores

```

### Test on https://www.ncbi.nlm.nih.gov/sra/SRX958739[accn]

1. **Download with fastq-dump --split-files**
2. **get 36 gene locations: https://pmc.ncbi.nlm.nih.gov/articles/PMC4676270/**
3. **UMI deduplication not working, comment out**

### Test RNAseq on https://www.ncbi.nlm.nih.gov/sra/SRX7223888[accn]