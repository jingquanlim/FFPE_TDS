rule fastqc:
    input:
        fastq="results/fastq/{sample}_{read}.fastq.gz",
    output:
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip",
    
    params:
        extra="--quiet",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "results/fastqc/{sample}_{read}.log",
    threads: 1
    wrapper:
        "v6.0.0/bio/fastqc"

rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_{read}_fastqc.{ext}",
            sample=SAMPLES,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    
    params:
        extra="--verbose --dirs",
    message:
        """--- Generating MultiQC report for seq data."""
    log:
        "results/multiqc/multiqc.log",
    wrapper:
        "v6.0.0/bio/multiqc"