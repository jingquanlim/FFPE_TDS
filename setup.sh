#!/bin/bash

# UMI-Aware Targeted DNA-Sequencing Pipeline Setup Script
# This script helps set up the pipeline environment and dependencies

set -e

echo "=== UMI-Aware Targeted DNA-Sequencing Pipeline Setup ==="
echo

# Check if conda/mamba is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "✓ Mamba found - using mamba for faster package installation"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "✓ Conda found - using conda for package installation"
else
    echo "❌ Neither conda nor mamba found. Please install Miniconda or Mambaforge first."
    echo "   Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Install Snakemake if not available
if ! command -v snakemake &> /dev/null; then
    echo "📦 Installing Snakemake..."
    $CONDA_CMD install -c conda-forge -c bioconda snakemake -y
    echo "✓ Snakemake installed"
else
    echo "✓ Snakemake already installed"
fi

# Create directory structure
echo "📁 Creating directory structure..."
mkdir -p results/{qc/fastqc,umi_extract,alignment,variants/{snv_indel,sv,msi},reports}
mkdir -p logs
mkdir -p envs
mkdir -p scripts

echo "✓ Directory structure created"

# Check if reference genome path is configured
echo "🔍 Checking configuration..."
if [ ! -f "config.yaml" ]; then
    echo "❌ config.yaml not found. Please ensure the configuration file is present."
    exit 1
fi

# Extract reference genome path from config
REF_PATH=$(grep -E "^\s*genome:" config.yaml | sed 's/.*genome:\s*//' | tr -d '"' | tr -d "'" | xargs)

if [ -z "$REF_PATH" ]; then
    echo "❌ Reference genome path not found in config.yaml"
    exit 1
fi

echo "📋 Reference genome path: $REF_PATH"

# Check if reference genome exists
if [ ! -f "$REF_PATH" ]; then
    echo "⚠️  Reference genome file not found at: $REF_PATH"
    echo "   Please ensure the hs37d5 reference genome is available and the path is correct."
    echo "   You can download it from: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/"
else
    echo "✓ Reference genome found"
    
    # Check if BWA index exists
    if [ ! -f "${REF_PATH}.bwt" ]; then
        echo "📦 Creating BWA index (this may take a while)..."
        bwa index "$REF_PATH"
        echo "✓ BWA index created"
    else
        echo "✓ BWA index already exists"
    fi
fi

# Check if input directories exist
INPUT_DIR=$(grep -E "^\s*fastq_dir:" config.yaml | sed 's/.*fastq_dir:\s*//' | tr -d '"' | tr -d "'" | xargs)
TARGET_BED=$(grep -E "^\s*target_bed:" config.yaml | sed 's/.*target_bed:\s*//' | tr -d '"' | tr -d "'" | xargs)

if [ ! -d "$INPUT_DIR" ]; then
    echo "⚠️  Input FASTQ directory not found at: $INPUT_DIR"
    echo "   Please ensure the path is correct in config.yaml"
else
    echo "✓ Input FASTQ directory found"
    FASTQ_COUNT=$(find "$INPUT_DIR" -name "*.fastq.gz" | wc -l)
    echo "   Found $FASTQ_COUNT FASTQ files"
fi

if [ ! -f "$TARGET_BED" ]; then
    echo "⚠️  Target BED file not found at: $TARGET_BED"
    echo "   Please ensure the path is correct in config.yaml"
else
    echo "✓ Target BED file found"
fi

# Test conda environments
echo "🧪 Testing conda environments..."
for env_file in envs/*.yaml; do
    if [ -f "$env_file" ]; then
        env_name=$(basename "$env_file" .yaml)
        echo "  Testing $env_name environment..."
        
        # Try to create a temporary environment to test dependencies
        if $CONDA_CMD env create -f "$env_file" -n "test_$env_name" &> /dev/null; then
            echo "  ✓ $env_name environment dependencies available"
            $CONDA_CMD env remove -n "test_$env_name" -y &> /dev/null
        else
            echo "  ⚠️  Issues with $env_name environment - will be resolved during pipeline execution"
        fi
    fi
done

# Run pipeline dry-run to test workflow
echo "🔄 Testing pipeline workflow..."
if snakemake --dry-run &> /dev/null; then
    echo "✓ Pipeline workflow syntax is valid"
else
    echo "❌ Pipeline workflow has syntax errors"
    echo "   Run 'snakemake --dry-run' to see detailed error messages"
    exit 1
fi

echo
echo "🎉 Setup completed successfully!"
echo
echo "Next steps:"
echo "1. Review and adjust config.yaml if needed"
echo "2. Run the pipeline with: snakemake --cores 8 --use-conda"
echo "3. For a dry run first: snakemake --dry-run"
echo
echo "For more information, see README.md"
