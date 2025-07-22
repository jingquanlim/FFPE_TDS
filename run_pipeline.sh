#!/bin/bash

# UMI-Aware Targeted DNA-Sequencing Pipeline Execution Script
# This script provides convenient ways to run the pipeline

set -e

# Default parameters
CORES=8
DRY_RUN=true
CLUSTER=false
CLUSTER_CMD=""
CONDA_PREFIX=""
UNLOCK=false

# Usage information
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "Options:"
    echo "  -c, --cores CORES        Number of CPU cores to use (default: 8)"
    echo "  -d, --dry-run           Perform a dry run without executing"
    echo "  -u, --unlock            Unlock the working directory"
    echo "  --cluster COMMAND       Run on cluster with specified command"
    echo "  --conda-prefix PATH     Specify conda prefix directory"
    echo "  -h, --help              Show this help message"
    echo
    echo "Examples:"
    echo "  $0                      # Run with default 8 cores"
    echo "  $0 -c 16                # Run with 16 cores"
    echo "  $0 -d                   # Dry run"
    echo "  $0 --cluster 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}'"
    echo "  $0 --unlock             # Unlock working directory"
    echo
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        --cluster)
            CLUSTER=true
            CLUSTER_CMD="$2"
            shift 2
            ;;
        --conda-prefix)
            CONDA_PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check if Snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "❌ Snakemake not found. Please install Snakemake first."
    echo "   conda install -c conda-forge -c bioconda snakemake"
    exit 1
fi

# Check if config file exists
if [ ! -f "config.yaml" ]; then
    echo "❌ config.yaml not found. Please ensure the configuration file is present."
    exit 1
fi

# Check if Snakefile exists
if [ ! -f "Snakefile" ]; then
    echo "❌ Snakefile not found. Please ensure the Snakefile is present."
    exit 1
fi

echo "=== UMI-Aware Targeted DNA-Sequencing Pipeline ==="
echo "Configuration:"
echo "  Cores: $CORES"
echo "  Dry run: $DRY_RUN"
echo "  Cluster: $CLUSTER"
echo "  Unlock: $UNLOCK"
echo

# Handle unlock
if [ "$UNLOCK" = true ]; then
    echo "🔓 Unlocking working directory..."
    snakemake --unlock
    echo "✓ Working directory unlocked"
    exit 0
fi

# Build snakemake command
SNAKEMAKE_CMD="snakemake --cores $CORES --use-conda"

# Add conda prefix if specified
if [ -n "$CONDA_PREFIX" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --conda-prefix $CONDA_PREFIX"
fi

# Add cluster command if specified
if [ "$CLUSTER" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cluster '$CLUSTER_CMD'"
fi

# Add dry run flag if specified
if [ "$DRY_RUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run"
    echo "🔄 Performing dry run..."
else
    echo "🚀 Starting pipeline execution..."
fi

# Add additional useful flags
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --printshellcmds logs/snakemake_stats.txt"

echo "Command: $SNAKEMAKE_CMD"
echo

# Create logs directory
mkdir -p logs

# Execute the pipeline
eval $SNAKEMAKE_CMD

if [ "$DRY_RUN" = false ]; then
    echo
    echo "🎉 Pipeline execution completed!"
    echo
    echo "📊 Key output files:"
    echo "  - Overall summary: results/reports/pipeline_summary.html"
    echo "  - Sample reports: results/reports/*_quantification_report.tsv"
    echo "  - Processed BAMs: results/alignment/*_processed.bam"
    echo "  - Variant calls: results/variants/"
    echo "  - Execution stats: logs/snakemake_stats.txt"
    echo
    echo "📖 Open results/reports/pipeline_summary.html in a web browser to view the summary report."
else
    echo "✅ Dry run completed successfully!"
fi
