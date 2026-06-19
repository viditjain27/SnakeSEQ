"""
merge_counts.py
────────────────────────────────────────────────────────────────────────
Called by Snakemake rule: merge_counts

Purpose:
    featureCounts produces ONE file per sample, each with columns:
        Geneid  Chr  Start  End  Strand  Length  <sample_bam_path>
    The last column (raw count) has a long messy name (the BAM path).

    This script reads all per-sample files, extracts just the
    Geneid + count column from each, and merges them into a single
    matrix:
                  SRR1039508  SRR1039509  SRR1039512 ...
        GeneID1        120         98          145
        GeneID2          5          3            8
        ...

    This count_matrix.csv is the direct input to PyDESeq2.
────────────────────────────────────────────────────────────────────────
"""

import pandas as pd

# Snakemake injects 'snakemake' object automatically when using script: directive
input_files = snakemake.input.counts
output_file = snakemake.output.matrix
log_file = snakemake.log[0]

with open(log_file, "w") as log:
    log.write(f"Merging {len(input_files)} count files...\n")

    merged = None

    for f in input_files:
        # Extract sample name from filename: results/counts/SRR1039508.counts.txt -> SRR1039508
        sample_name = f.split("/")[-1].replace(".counts.txt", "")

        # featureCounts output: first line is a comment (# Program:featureCounts ...)
        # second line is the real header
        df = pd.read_csv(f, sep="\t", comment="#")

        # Columns are: Geneid, Chr, Start, End, Strand, Length, <bam_path>
        # The count column is always the LAST column
        count_col = df.columns[-1]

        sample_df = df[["Geneid", count_col]].copy()
        sample_df.columns = ["Geneid", sample_name]
        sample_df = sample_df.set_index("Geneid")

        log.write(f"  {sample_name}: {sample_df.shape[0]} genes, "
                  f"total reads = {sample_df[sample_name].sum()}\n")

        if merged is None:
            merged = sample_df
        else:
            merged = merged.join(sample_df, how="outer")

    # Fill any missing values with 0 (shouldn't happen, but safe)
    merged = merged.fillna(0).astype(int)

    log.write(f"\nFinal count matrix shape: {merged.shape[0]} genes x {merged.shape[1]} samples\n")

    merged.to_csv(output_file)
    log.write(f"Written to {output_file}\n")
