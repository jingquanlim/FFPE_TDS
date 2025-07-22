#!/usr/bin/env python3

import os
import glob
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Define samples based on input FASTQ files

def get_rna_samples():
    """Extract sample names from input RNA FASTQ files."""
    input_dir = config["rnaseq"]["input"]["fastq_dir"]
    r1_files = glob.glob(os.path.join(input_dir, "*_R1_*.fastq.gz"))
    samples = []
    for r1 in r1_files:
        basename = os.path.basename(r1)
        # Extract sample name (assuming format: SAMPLE_R1_001.fastq.gz)
        sample_name = basename.split("_R1_")[0]
        samples.append(sample_name)
    return samples

RNA_SAMPLES = get_rna_samples()

# Define final output files
rule all:
    input:
        # Quality control reports
        expand("results/rnaseq/fastqc/{sample}_R1_001_fastqc.html", sample=RNA_SAMPLES),
        expand("results/rnaseq/fastqc/{sample}_R2_001_fastqc.html", sample=RNA_SAMPLES),
        
        # RNA-seq outputs
        expand("results/rnaseq/quantification/{sample}.genes.results", sample=RNA_SAMPLES),

        expand("results/rnaseq/quantification/{sample}.isoforms.results", sample=RNA_SAMPLES),
        expand("results/rnaseq/fusion/{sample}_cicero.fusions.tsv", sample=RNA_SAMPLES)


# Quality control with FastQC
rule fastqc:
    input:
        r1=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
        r2=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R2_001.fastq.gz"
    output:
        r1_html="results/rnaseq/fastqc/{sample}_R1_001_fastqc.html",
        r1_zip="results/rnaseq/fastqc/{sample}_R1_001_fastqc.zip",
        r2_html="results/rnaseq/fastqc/{sample}_R2_001_fastqc.html",
        r2_zip="results/rnaseq/fastqc/{sample}_R2_001_fastqc.zip"
    params:
        outdir="results/rnaseq/fastqc"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} -o {params.outdir} -t {threads}
        """


# MeRNA for rRNA removal
# rule merna:
#     input:
#         r1=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
#         r2=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R2_001.fastq.gz"
#     output:
#         r1="results/rnaseq/merna/{sample}_R1_norna.fastq.gz",
#         r2="results/rnaseq/merna/{sample}_R2_norna.fastq.gz"
#     threads: 8
#     conda:
#         "envs/rnaseq.yaml"
#     shell:
#         """
#         mkdir -p results/rnaseq/merna
#         merna -1 {input.r1} -2 {input.r2} -o results/rnaseq/merna/{wildcards.sample} --threads {threads}
#         mv results/rnaseq/merna/{wildcards.sample}_1.fastq.gz {output.r1}
#         mv results/rnaseq/merna/{wildcards.sample}_2.fastq.gz {output.r2}
#         """

# BBSplit for contaminant removal
# rule bbsplit:
#     input:
#         r1="results/rnaseq/merna/{sample}_R1_norna.fastq.gz",
#         r2="results/rnaseq/merna/{sample}_R2_norna.fastq.gz",
#         human_ref=config["rnaseq"]["bbsplit"]["human_ref"],
#         mouse_ref=config["rnaseq"]["bbsplit"]["mouse_ref"]
#     output:
#         r1="results/rnaseq/bbsplit/{sample}_R1_clean.fastq.gz",
#         r2="results/rnaseq/bbsplit/{sample}_R2_clean.fastq.gz"
#     threads: 8
#     conda:
#         "envs/rnaseq.yaml"
#     shell:
#         """
#         mkdir -p results/rnaseq/bbsplit
#         bbsplit.sh -Xmx40g threads={threads} ref_human={input.human_ref} ref_mouse={input.mouse_ref} \
#                    in1={input.r1} in2={input.r2} \
#                    basename=results/rnaseq/bbsplit/{wildcards.sample}_%.fastq.gz
#         mv results/rnaseq/bbsplit/{wildcards.sample}_human.fastq.gz {output.r1}
#         mv results/rnaseq/bbsplit/{wildcards.sample}_human_2.fastq.gz {output.r2}
#         """

# STAR alignment
rule star_align_rnaseq:
    input:
        r1=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
        r2=config["rnaseq"]["input"]["fastq_dir"] + "/{sample}_R2_001.fastq.gz",
        stardir=config["rnaseq"]["STAR"]["ref"]
    output:
        bam="results/rnaseq/alignment/{sample}_Aligned.toTranscriptome.out.bam",
        chimeric_junction="results/rnaseq/alignment/{sample}_Chimeric.out.junction"
    threads: 12
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/rnaseq/alignment
        STAR --genomeDir {input.stardir} \
             --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
             --chimOutType Junctions \
             --chimJunctionOverhangMin 10 \
             --chimScoreMin 1 \
             --chimScoreDropMax 30 \
             --chimScoreJunctionNonGTAG 0 \
             --chimScoreSeparation 1 \
             --alignSJstitchMismatchNmax 5 -1 5 5 \
             --chimSegmentReadGapMax 3 \
             --outFileNamePrefix results/rnaseq/alignment/{wildcards.sample}_ \
             --runThreadN {threads}
        # Create an empty chimeric file if STAR didn't find any chimeric reads
        [[ -f {output.chimeric_junction} ]] || touch {output.chimeric_junction}
        """

# RSEM quantification
rule rsem_quant:
    input:
        bam="results/rnaseq/alignment/{sample}_Aligned.toTranscriptome.out.bam",
        
    params:
        rsem_ref=config["rnaseq"]["rsem"]["ref"]
    output:
        genes="results/rnaseq/quantification/{sample}.genes.results",
        isoforms="results/rnaseq/quantification/{sample}.isoforms.results"
    threads: 8
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/rnaseq/quantification
        rsem-calculate-expression --bam --no-bam-output -p {threads} \
                                  --paired-end {input.bam} \
                                  {params.rsem_ref} \
                                  results/rnaseq/quantification/{wildcards.sample}
        """

# CICERO for fusion detection
rule cicero:
    input:
        chimeric_junction="results/rnaseq/alignment/{sample}_Chimeric.out.junction",
        fusion_annot=config["rnaseq"]["cicero"]["fusion_annot"]
    output:
        fusions="results/rnaseq/fusion/{sample}_cicero.fusions.tsv"
    threads: 1
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        mkdir -p results/rnaseq/fusion
        cicero.py -t {input.chimeric_junction} -o {output.fusions} -a {input.fusion_annot}
        """
