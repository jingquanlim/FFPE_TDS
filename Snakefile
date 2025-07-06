#!/usr/bin/env python3

import os
import glob
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Define samples based on input FASTQ files
def get_samples():
    """Extract sample names from input FASTQ files."""
    input_dir = config["input"]["fastq_dir"]
    r1_files = glob.glob(os.path.join(input_dir, "*_R1_*.fastq.gz"))
    samples = []
    for r1 in r1_files:
        basename = os.path.basename(r1)
        # Extract sample name (assuming format: SAMPLE_R1_001.fastq.gz)
        sample_name = basename.split("_R1_")[0]
        samples.append(sample_name)
    return samples

SAMPLES = get_samples()

# Define final output files
rule all:
    input:
        # Quality control reports
        expand("results/qc/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("results/qc/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        
        # Final processed BAM files
        expand("results/alignment/{sample}_processed.bam", sample=SAMPLES),
        expand("results/alignment/{sample}_processed.bam.bai", sample=SAMPLES),
        
        # Variant calling results
        expand("results/variants/snv_indel/{sample}_variants.vcf", sample=SAMPLES),
        expand("results/variants/sv/{sample}_sv.vcf", sample=SAMPLES),
        expand("results/variants/msi/{sample}_msi.txt", sample=SAMPLES),
        
        # Quantification reports
        expand("results/reports/{sample}_quantification_report.tsv", sample=SAMPLES),
        
        # Summary report
        "results/reports/pipeline_summary.html"

# Quality control with FastQC
rule fastqc:
    input:
        r1=config["input"]["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
        r2=config["input"]["fastq_dir"] + "/{sample}_R2_001.fastq.gz"
    output:
        r1_html="results/qc/fastqc/{sample}_R1_fastqc.html",
        r1_zip="results/qc/fastqc/{sample}_R1_fastqc.zip",
        r2_html="results/qc/fastqc/{sample}_R2_fastqc.html",
        r2_zip="results/qc/fastqc/{sample}_R2_fastqc.zip"
    params:
        outdir="results/qc/fastqc"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} -o {params.outdir} -t {threads}
        """

# UMI extraction
rule umi_extract:
    input:
        r1=config["input"]["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
        r2=config["input"]["fastq_dir"] + "/{sample}_R2_001.fastq.gz"
    output:
        r1="results/umi_extract/{sample}_R1_umi.fastq.gz",
        r2="results/umi_extract/{sample}_R2_umi.fastq.gz"
    params:
        bc_pattern=config["umi"]["barcode_pattern"],
        bc_pattern2=config["umi"]["barcode_pattern2"] if config["umi"]["barcode_pattern2"] else ""
    threads: 1
    conda:
        "envs/umi_tools.yaml"
    shell:
        """
        mkdir -p results/umi_extract
        if [ -n "{params.bc_pattern2}" ]; then
            umi_tools extract --bc-pattern={params.bc_pattern} \
                             --bc-pattern2={params.bc_pattern2} \
                             --stdin {input.r1} \
                             --stdout {output.r1} \
                             --read2-in {input.r2} \
                             --read2-out {output.r2}
        else
            umi_tools extract --bc-pattern={params.bc_pattern} \
                             --stdin {input.r1} \
                             --stdout {output.r1} \
                             --read2-in {input.r2} \
                             --read2-out {output.r2}
        fi
        """

# BWA alignment
rule bwa_align:
    input:
        r1="results/umi_extract/{sample}_R1_umi.fastq.gz",
        r2="results/umi_extract/{sample}_R2_umi.fastq.gz",
        ref=config["reference"]["genome"]
    output:
        bam="results/alignment/{sample}_aligned.bam"
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    threads: 8
    conda:
        "envs/alignment.yaml"
    shell:
        """
        mkdir -p results/alignment
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam}
        """

# Sort and index BAM
rule sort_index_bam:
    input:
        bam="results/alignment/{sample}_aligned.bam"
    output:
        bam="results/alignment/{sample}_sorted.bam",
        bai="results/alignment/{sample}_sorted.bam.bai"
    threads: 4
    conda:
        "envs/alignment.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

# UMI deduplication
rule umi_dedup:
    input:
        bam="results/alignment/{sample}_sorted.bam",
        bai="results/alignment/{sample}_sorted.bam.bai"
    output:
        bam="results/alignment/{sample}_processed.bam",
        bai="results/alignment/{sample}_processed.bam.bai"
    threads: 2
    conda:
        "envs/umi_tools.yaml"
    shell:
        """
        umi_tools dedup --stdin={input.bam} --stdout={output.bam}
        samtools index {output.bam}
        """

# SNV and Indel calling with FreeBayes
rule freebayes_call:
    input:
        bam="results/alignment/{sample}_processed.bam",
        bai="results/alignment/{sample}_processed.bam.bai",
        ref=config["reference"]["genome"],
        bed=config["input"]["target_bed"]
    output:
        vcf="results/variants/snv_indel/{sample}_variants.vcf"
    threads: 4
    conda:
        "envs/variant_calling.yaml"
    shell:
        """
        mkdir -p results/variants/snv_indel
        freebayes -f {input.ref} -t {input.bed} {input.bam} > {output.vcf}
        """

# Structural variant calling with Manta
rule manta_call:
    input:
        bam="results/alignment/{sample}_processed.bam",
        bai="results/alignment/{sample}_processed.bam.bai",
        ref=config["reference"]["genome"]
    output:
        vcf="results/variants/sv/{sample}_sv.vcf",
        dir=directory("results/variants/sv/{sample}_manta_work")
    threads: 8
    conda:
        "envs/variant_calling.yaml"
    shell:
        """
        mkdir -p results/variants/sv
        
        # Configure Manta
        configManta.py --bam {input.bam} --referenceFasta {input.ref} \
                       --runDir {output.dir}
        
        # Run Manta
        {output.dir}/runWorkflow.py -m local -j {threads}
        
        # Copy results
        cp {output.dir}/results/variants/diploidSV.vcf.gz {output.vcf}.gz
        gunzip {output.vcf}.gz
        """

# MSI calling with MSIsensor-pro
rule msi_call:
    input:
        bam="results/alignment/{sample}_processed.bam",
        bai="results/alignment/{sample}_processed.bam.bai",
        ref=config["reference"]["genome"]
    output:
        msi="results/variants/msi/{sample}_msi.txt"
    threads: 4
    conda:
        "envs/msi.yaml"
    shell:
        """
        mkdir -p results/variants/msi
        
        # Create MSI list if not exists
        if [ ! -f results/variants/msi/msi_list.txt ]; then
            msisensor-pro scan -d {input.ref} -o results/variants/msi/msi_list.txt
        fi
        
        # Run MSI analysis
        msisensor-pro msi -d results/variants/msi/msi_list.txt \
                          -t {input.bam} \
                          -o results/variants/msi/{wildcards.sample} \
                          -b {threads}
        
        # Rename output file
        mv results/variants/msi/{wildcards.sample} {output.msi}
        """

# Generate quantification report
rule quantify_variants:
    input:
        snv_vcf="results/variants/snv_indel/{sample}_variants.vcf",
        sv_vcf="results/variants/sv/{sample}_sv.vcf",
        msi_txt="results/variants/msi/{sample}_msi.txt"
    output:
        report="results/reports/{sample}_quantification_report.tsv"
    conda:
        "envs/reporting.yaml"
    script:
        "scripts/quantify_variants.py"

# Generate summary report
rule generate_summary:
    input:
        reports=expand("results/reports/{sample}_quantification_report.tsv", sample=SAMPLES)
    output:
        summary="results/reports/pipeline_summary.html"
    conda:
        "envs/reporting.yaml"
    script:
        "scripts/generate_summary.py"
