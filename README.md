# FFPE_TDS
A Snakemake pipeline to detect markers on targeted DNA/RNA sequencing data. DNA Markers includes SNV, Indel, SV, MSI. RNA Markers includes Fusion, TPM


## DNA
- fastQC: quality control
- bwa-mem: alignment
- sambamba: sort
- freebayes: call variants (SNV, Indel etc.)
- mantaSV: call SV
- MSIsensor-pro: MSI


## RNA
- fastQC: quality control
- STAR: alignment
- samtools: sort
- qualimap: quality control
- FusionCatcher/CICERO(slow, but manage to call ITD): gene fusions
- featurecounts: TPM
