SAMPLES = ["A", "B", "C"]

rule all:
    input:
        expand("fastq/{sample}.trimmed.fastq.gz", sample=SAMPLES),
        "multiqc_report.html"

# Rule to process samples
rule trim_fastq:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastq/{sample}.trimmed.fastq.gz"
    shell:
        "fastp -i {input} -o {output}"

# Aggregate QC report
rule multiqc:
    input:
        expand("fastq/{sample}.trimmed.fastq.gz", sample=SAMPLES)
    output:
        "multiqc_report.html"
    shell:
        "multiqc . -o ."