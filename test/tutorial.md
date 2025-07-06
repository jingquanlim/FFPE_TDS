# UMI-Aware Targeted DNA-Sequencing Pipeline

A comprehensive Snakemake pipeline for processing targeted DNA-sequencing data with UMI (Unique Molecular Identifier) support. This pipeline accurately detects and quantifies Single Nucleotide Variants (SNVs), Insertions/Deletions (Indels), Structural Variants (SVs), and Microsatellite Instability (MSI).

## Features

- **UMI-aware processing** using umi_tools for bias reduction
- **Multi-threaded execution** for improved performance
- **Comprehensive variant detection**: SNVs, INDELs, SVs, and MSI
- **Quality control** with FastQC reporting
- **Automated quantification** and summary reporting
- **Conda environment management** for reproducible execution

## Pipeline Overview

The pipeline consists of the following main steps:

1. **Quality Control** - FastQC analysis of raw reads
2. **UMI Extraction** - Extract UMIs from reads using umi_tools
3. **Alignment** - BWA-MEM alignment to hs37d5 reference
4. **Post-processing** - Sorting, indexing, and UMI-based deduplication
5. **Variant Calling**:
   - SNV/INDEL calling with FreeBayes
   - Structural variant calling with Manta
   - MSI analysis with MSIsensor-pro
6. **Quantification** - Generate summary reports and statistics


## Tutorial with test data

This tutorial demonstrates how to run the pipeline using the specified test dataset.

### 1. Setup Environment

```bash
# Create a working directory
mkdir umi_pipeline_test
cd umi_pipeline_test

# Copy pipeline files
cp -r /path/to/pipeline/* .
```

### 2. Configure for Test Data

Edit `config.yaml`:

```yaml
input:
  fastq_dir: "/data/ocklab2/jingquan/binlab2/Targeted_TWIST_NHL_TDSB9/input_files"
  target_bed: "/data/ocklab2/jingquan/binlab2/Targeted_TWIST_NHL_TDSB9/TWIST_NHL_rev3.bed"

reference:
  genome: "/path/to/hs37d5/hs37d5.fa"  # Update this path

umi:
  barcode_pattern: "NNNNNNNN"  # Adjust based on actual UMI structure
```

### 3. Run the Pipeline

```bash
# Test the workflow
snakemake --dry-run

# Run with conda environments
snakemake --cores 8 --use-conda
```

### 4. Expected Outputs

The pipeline will generate:

```
results/
├── qc/
│   └── fastqc/           # FastQC reports
├── umi_extract/          # UMI-extracted reads
├── alignment/            # Aligned and processed BAM files
├── variants/
│   ├── snv_indel/        # SNV/INDEL VCF files
│   ├── sv/               # Structural variant calls
│   └── msi/              # MSI analysis results
└── reports/
    ├── *_quantification_report.tsv  # Per-sample reports
    └── pipeline_summary.html        # Overall summary
```

## Output Files

### Key Output Files

- **`results/reports/pipeline_summary.html`** - Comprehensive HTML summary report
- **`results/reports/*_quantification_report.tsv`** - Per-sample quantification data
- **`results/alignment/*_processed.bam`** - Final processed BAM files
- **`results/variants/snv_indel/*_variants.vcf`** - SNV/INDEL calls
- **`results/variants/sv/*_sv.vcf`** - Structural variant calls
- **`results/variants/msi/*_msi.txt`** - MSI analysis results

### Report Contents

The quantification reports include:
- Total counts of each variant type per sample
- Allele frequencies for SNVs and INDELs
- Structural variant characteristics
- MSI status and scores

## Troubleshooting

### Common Issues

1. **Reference genome not found**:
   - Ensure the hs37d5 reference path is correct in `config.yaml`
   - Check that BWA index files exist

2. **UMI extraction fails**:
   - Verify the UMI pattern matches your library preparation
   - Check FASTQ file formats and headers

3. **Memory issues**:
   - Increase memory allocation in `config.yaml`
   - Use `--resources mem_mb=16000` for high-memory rules

4. **Conda environment issues**:
   - Clear conda cache: `conda clean --all`
   - Use `--conda-cleanup-pkgs cache` with snakemake

### Performance Optimization

- Use `--cores` to specify maximum CPU cores
- Adjust thread counts in `config.yaml` based on your system
- Consider using `mamba` instead of `conda` for faster dependency resolution

## File Structure

```
umi-targeted-sequencing-pipeline/
├── Snakefile                    # Main pipeline workflow
├── config.yaml                  # Configuration file
├── envs/                        # Conda environment files
│   ├── qc.yaml
│   ├── umi_tools.yaml
│   ├── alignment.yaml
│   ├── variant_calling.yaml
│   ├── msi.yaml
│   └── reporting.yaml
├── scripts/                     # Analysis scripts
│   ├── quantify_variants.py
│   └── generate_summary.py
└── README.md                    # This file
```

## Citation

If you use this pipeline in your research, please cite the relevant tools:

- **Snakemake**: Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine.
- **umi_tools**: Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers.
- **BWA**: Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform.
- **FreeBayes**: Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing.
- **Manta**: Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels.
- **MSIsensor-pro**: Jia, P., et al. (2020). MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability.

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review Snakemake documentation: https://snakemake.readthedocs.io/
3. Open an issue in the project repository

## License

This pipeline is released under the MIT License.
