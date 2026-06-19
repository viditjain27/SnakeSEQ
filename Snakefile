# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SnakeSEQ — Automated RNA-seq Differential Expression Analysis Pipeline     ║
# ║  Author  : Vidit Jain                                                        ║
# ║  Dataset : GSE52778 — Airway Smooth Muscle Cells (Human, Paired-end)        ║
# ║  Run     : snakemake --use-conda --cores 4                                   ║
# ║  Dry-run : snakemake --use-conda -n --cores 4                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

import pandas as pd
import os

# ── Load configuration and sample sheet ───────────────────────────────────────
configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples"], sep="\t", index_col="sample")
SAMPLES    = samples_df.index.tolist()

# ══════════════════════════════════════════════════════════════════════════════
# rule all — Defines FINAL outputs. Snakemake works BACKWARDS from these.
# ══════════════════════════════════════════════════════════════════════════════
rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/deseq2/volcano_plot.png",
        "results/deseq2/heatmap.png",
        "results/deseq2/pca_plot.png",
        "results/deseq2/results_table.csv"


# ══════════════════════════════════════════════════════════════════════════════
# RULE 1: Download raw FASTQ from NCBI SRA
# ══════════════════════════════════════════════════════════════════════════════
rule fasterq_dump:
    output:
        r1 = "results/fastq/raw/{sample}_1.fastq.gz",
        r2 = "results/fastq/raw/{sample}_2.fastq.gz"
    params:
        outdir = "results/fastq/raw",
        sradir = "results/sra"
    conda: "envs/download.yaml"
    threads: 4
    log: "logs/fasterq_dump/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} {params.sradir}
        prefetch {wildcards.sample} -O {params.sradir} 2> {log}
        fasterq-dump {params.sradir}/{wildcards.sample} \
            --split-3 \
            --outdir {params.outdir} \
            --threads {threads} 2>> {log}
        gzip {params.outdir}/{wildcards.sample}_1.fastq
        gzip {params.outdir}/{wildcards.sample}_2.fastq
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 2: Quality control — raw reads
# ══════════════════════════════════════════════════════════════════════════════
rule fastqc_raw:
    input:
        r1 = "results/fastq/raw/{sample}_1.fastq.gz",
        r2 = "results/fastq/raw/{sample}_2.fastq.gz"
    output:
        html1 = "results/qc/raw/{sample}_1_fastqc.html",
        html2 = "results/qc/raw/{sample}_2_fastqc.html",
        zip1  = "results/qc/raw/{sample}_1_fastqc.zip",
        zip2  = "results/qc/raw/{sample}_2_fastqc.zip"
    conda: "envs/qc.yaml"
    threads: 2
    log: "logs/fastqc_raw/{sample}.log"
    shell:
        """
        mkdir -p results/qc/raw
        fastqc -t {threads} -o results/qc/raw {input.r1} {input.r2} 2> {log}
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 3: Adapter trimming with fastp
# ══════════════════════════════════════════════════════════════════════════════
rule trim_fastp:
    input:
        r1 = "results/fastq/raw/{sample}_1.fastq.gz",
        r2 = "results/fastq/raw/{sample}_2.fastq.gz"
    output:
        r1   = "results/fastq/trimmed/{sample}_1.trimmed.fastq.gz",
        r2   = "results/fastq/trimmed/{sample}_2.trimmed.fastq.gz",
        html = "results/qc/fastp/{sample}_fastp.html",
        json = "results/qc/fastp/{sample}_fastp.json"
    conda: "envs/qc.yaml"
    threads: 4
    log: "logs/fastp/{sample}.log"
    shell:
        """
        mkdir -p results/fastq/trimmed results/qc/fastp
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --length_required 36 \
            --qualified_quality_phred 20 \
            --detect_adapter_for_pe \
            2> {log}
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 4: Quality control — trimmed reads
# ══════════════════════════════════════════════════════════════════════════════
rule fastqc_trimmed:
    input:
        r1 = "results/fastq/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = "results/fastq/trimmed/{sample}_2.trimmed.fastq.gz"
    output:
        html1 = "results/qc/trimmed/{sample}_1.trimmed_fastqc.html",
        html2 = "results/qc/trimmed/{sample}_2.trimmed_fastqc.html",
        zip1  = "results/qc/trimmed/{sample}_1.trimmed_fastqc.zip",
        zip2  = "results/qc/trimmed/{sample}_2.trimmed_fastqc.zip"
    conda: "envs/qc.yaml"
    threads: 2
    log: "logs/fastqc_trimmed/{sample}.log"
    shell:
        """
        mkdir -p results/qc/trimmed
        fastqc -t {threads} -o results/qc/trimmed {input.r1} {input.r2} 2> {log}
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 5: Splice-aware alignment with HISAT2
# ══════════════════════════════════════════════════════════════════════════════
rule hisat2_align:
    input:
        r1 = "results/fastq/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = "results/fastq/trimmed/{sample}_2.trimmed.fastq.gz"
    output:
        bam = temp("results/aligned/{sample}.unsorted.bam")
    params:
        index = config["genome_index"]
    conda: "envs/align.yaml"
    threads: config["threads"]
    log: "logs/hisat2/{sample}.log"
    shell:
        """
        mkdir -p results/aligned
        hisat2 \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --dta \
            --new-summary \
            2> {log} \
        | samtools view -bS -o {output.bam} -
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 6: Sort BAM by coordinate
# ══════════════════════════════════════════════════════════════════════════════
rule samtools_sort:
    input:  "results/aligned/{sample}.unsorted.bam"
    output: "results/aligned/{sample}.sorted.bam"
    conda: "envs/align.yaml"
    threads: 4
    log: "logs/samtools_sort/{sample}.log"
    shell: "samtools sort -@ {threads} -o {output} {input} 2> {log}"


# ══════════════════════════════════════════════════════════════════════════════
# RULE 7: Index sorted BAM
# ══════════════════════════════════════════════════════════════════════════════
rule samtools_index:
    input:  "results/aligned/{sample}.sorted.bam"
    output: "results/aligned/{sample}.sorted.bam.bai"
    conda: "envs/align.yaml"
    log: "logs/samtools_index/{sample}.log"
    shell: "samtools index {input} 2> {log}"


# ══════════════════════════════════════════════════════════════════════════════
# RULE 8: Count reads per gene with featureCounts
# ══════════════════════════════════════════════════════════════════════════════
rule featurecounts:
    input:
        bam = "results/aligned/{sample}.sorted.bam",
        bai = "results/aligned/{sample}.sorted.bam.bai",
        gtf = config["gtf"]
    output:
        counts  = "results/counts/{sample}.counts.txt",
        summary = "results/counts/{sample}.counts.txt.summary"
    params:
        strand = config["strandedness"]
    conda: "envs/counts.yaml"
    threads: 4
    log: "logs/featurecounts/{sample}.log"
    shell:
        """
        mkdir -p results/counts
        featureCounts \
            -T {threads} \
            -p --countReadPairs \
            -s {params.strand} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam} \
            2> {log}
        """


# ══════════════════════════════════════════════════════════════════════════════
# RULE 9: Merge all per-sample counts into one matrix
# ══════════════════════════════════════════════════════════════════════════════
rule merge_counts:
    input:
        counts = expand("results/counts/{sample}.counts.txt", sample=SAMPLES)
    output:
        matrix = "results/counts/count_matrix.csv"
    conda: "envs/deseq2.yaml"
    log: "logs/merge_counts.log"
    script: "scripts/merge_counts.py"


# ══════════════════════════════════════════════════════════════════════════════
# RULE 10: DEA with PyDESeq2 + all visualizations
# ══════════════════════════════════════════════════════════════════════════════
rule deseq2_analysis:
    input:
        matrix  = "results/counts/count_matrix.csv",
        samples = config["samples"]
    output:
        results  = "results/deseq2/results_table.csv",
        norm_mat = "results/deseq2/normalized_counts.csv",
        volcano  = "results/deseq2/volcano_plot.png",
        heatmap  = "results/deseq2/heatmap.png",
        pca      = "results/deseq2/pca_plot.png"
    params:
        numerator   = config["contrast"]["numerator"],
        denominator = config["contrast"]["denominator"],
        fdr         = 0.05,
        lfc         = 1.0
    conda: "envs/deseq2.yaml"
    log: "logs/deseq2.log"
    script: "scripts/run_deseq2.py"


# ══════════════════════════════════════════════════════════════════════════════
# RULE 11: MultiQC — aggregate all QC reports
# ══════════════════════════════════════════════════════════════════════════════
rule multiqc:
    input:
        expand("results/qc/raw/{sample}_1_fastqc.zip",             sample=SAMPLES),
        expand("results/qc/raw/{sample}_2_fastqc.zip",             sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_1.trimmed_fastqc.zip", sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_2.trimmed_fastqc.zip", sample=SAMPLES),
        expand("results/qc/fastp/{sample}_fastp.json",             sample=SAMPLES),
        expand("logs/hisat2/{sample}.log",                         sample=SAMPLES),
        expand("results/counts/{sample}.counts.txt.summary",       sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    conda: "envs/qc.yaml"
    log: "logs/multiqc.log"
    shell:
        """
        mkdir -p results/multiqc
        multiqc results/qc/ results/qc/fastp/ logs/hisat2/ results/counts/ \
            -o results/multiqc --force 2> {log}
        """
