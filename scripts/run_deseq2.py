"""
run_deseq2.py
────────────────────────────────────────────────────────────────────────
Called by Snakemake rule: deseq2_analysis

Purpose:
    1. Load the merged count matrix + sample metadata
    2. Run PyDESeq2 (negative binomial GLM + Wald test) to find
       differentially expressed genes between conditions
    3. Generate 3 plots:
         - Volcano plot   (log2FC vs -log10 padj)
         - Heatmap        (top 50 variable genes, z-scored)
         - PCA plot       (sample clustering)
    4. Save results table + normalized counts as CSV

PyDESeq2 replicates R's DESeq2 algorithm:
    - Estimates size factors (library size normalization)
    - Estimates per-gene dispersion (variance modeling)
    - Fits negative binomial GLM per gene
    - Wald test for log2FoldChange != 0
    - Benjamini-Hochberg FDR correction -> padj
────────────────────────────────────────────────────────────────────────
"""

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # no display needed, just save files
import matplotlib.pyplot as plt
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

# ── Snakemake I/O ──────────────────────────────────────────────────────
count_matrix_file = snakemake.input.matrix
samples_file       = snakemake.input.samples

out_results  = snakemake.output.results
out_norm     = snakemake.output.norm_mat
out_volcano  = snakemake.output.volcano
out_heatmap  = snakemake.output.heatmap
out_pca      = snakemake.output.pca

numerator    = snakemake.params.numerator
denominator  = snakemake.params.denominator
fdr_cutoff   = snakemake.params.fdr
lfc_cutoff   = snakemake.params.lfc

log_file = snakemake.log[0]
log = open(log_file, "w")

def logp(msg):
    log.write(msg + "\n")
    log.flush()

# ════════════════════════════════════════════════════════════════════
# STEP 1 — Load data
# ════════════════════════════════════════════════════════════════════
logp("Loading count matrix and sample metadata...")

counts_df = pd.read_csv(count_matrix_file, index_col=0)   # genes x samples
metadata  = pd.read_csv(samples_file, sep="\t", index_col="sample")

# PyDESeq2 expects: rows = samples, columns = genes
counts_df = counts_df.T
counts_df = counts_df.loc[metadata.index]  # ensure matching order

# Remove genes with zero counts across ALL samples
counts_df = counts_df.loc[:, counts_df.sum(axis=0) > 0]

logp(f"Count matrix: {counts_df.shape[0]} samples x {counts_df.shape[1]} genes (after removing all-zero genes)")
logp(f"Conditions: {metadata['condition'].value_counts().to_dict()}")

# ════════════════════════════════════════════════════════════════════
# STEP 2 — Run PyDESeq2
# ════════════════════════════════════════════════════════════════════
logp("\nFitting DESeq2 model (size factors + dispersion + GLM)...")

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
)
dds.deseq2()

logp("Running Wald test for contrast: "
     f"condition {numerator} vs {denominator}")

stat_res = DeseqStats(
    dds,
    contrast=["condition", numerator, denominator],
)
stat_res.summary()

results_df = stat_res.results_df.copy()
results_df = results_df.sort_values("padj")

n_sig = ((results_df["padj"] < fdr_cutoff) &
         (results_df["log2FoldChange"].abs() > lfc_cutoff)).sum()
logp(f"\nSignificant DE genes (padj < {fdr_cutoff} & |log2FC| > {lfc_cutoff}): {n_sig}")

results_df.to_csv(out_results)
logp(f"Results table written to {out_results}")

# Normalized counts (size-factor adjusted) for heatmap/PCA
normalized_counts = pd.DataFrame(
    dds.layers["normed_counts"],
    index=counts_df.index,
    columns=counts_df.columns,
)
normalized_counts.T.to_csv(out_norm)
logp(f"Normalized counts written to {out_norm}")

# ════════════════════════════════════════════════════════════════════
# STEP 3 — Volcano Plot
# ════════════════════════════════════════════════════════════════════
logp("\nGenerating volcano plot...")

plot_df = results_df.copy()
plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].replace(0, 1e-300))
plot_df["significant"] = (
    (plot_df["padj"] < fdr_cutoff) & (plot_df["log2FoldChange"].abs() > lfc_cutoff)
)

plt.figure(figsize=(9, 7))
colors = plot_df["significant"].map({True: "#E63946", False: "#A8A8A8"})
plt.scatter(
    plot_df["log2FoldChange"], plot_df["neg_log10_padj"],
    c=colors, s=12, alpha=0.6, edgecolors="none"
)
plt.axvline(lfc_cutoff, color="grey", linestyle="--", linewidth=0.8)
plt.axvline(-lfc_cutoff, color="grey", linestyle="--", linewidth=0.8)
plt.axhline(-np.log10(fdr_cutoff), color="grey", linestyle="--", linewidth=0.8)

# Label top 10 most significant genes
top_genes = plot_df[plot_df["significant"]].nsmallest(10, "padj")
for gene_id, row in top_genes.iterrows():
    plt.annotate(
        gene_id,
        (row["log2FoldChange"], row["neg_log10_padj"]),
        fontsize=7, alpha=0.8,
        xytext=(3, 3), textcoords="offset points"
    )

plt.xlabel("log2 Fold Change", fontsize=12)
plt.ylabel("-log10(adjusted p-value)", fontsize=12)
plt.title(f"Volcano Plot: {numerator} vs {denominator}\n"
          f"{n_sig} significant genes (padj<{fdr_cutoff}, |log2FC|>{lfc_cutoff})",
          fontsize=12)
plt.tight_layout()
plt.savefig(out_volcano, dpi=150)
plt.close()
logp(f"Volcano plot saved to {out_volcano}")

# ════════════════════════════════════════════════════════════════════
# STEP 4 — Heatmap of top 50 variable genes
# ════════════════════════════════════════════════════════════════════
logp("\nGenerating heatmap of top 50 most variable genes...")

# log-transform normalized counts
log_norm = np.log2(normalized_counts + 1).T  # Transpose to genes x samples

# Pick top 50 genes by variance across samples
top_var_genes = log_norm.var(axis=1).sort_values(ascending=False).head(50).index
heatmap_data = log_norm.loc[top_var_genes]

# z-score per gene (row-wise) so each gene is on comparable scale
heatmap_z = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(
    heatmap_data.std(axis=1).replace(0, 1), axis=0
)

# Column colors by condition
condition_colors = metadata["condition"].map(
    {numerator: "#E63946", denominator: "#457B9D"}
)

plt.figure(figsize=(10, 12))
g = sns.clustermap(
    heatmap_z,
    cmap="vlag",
    col_colors=condition_colors,
    figsize=(10, 12),
    yticklabels=True,
    xticklabels=True,
    cbar_kws={"label": "Row z-score"},
)
g.fig.suptitle("Top 50 Most Variable Genes (z-scored log2 normalized counts)",
               y=1.02, fontsize=12)
g.savefig(out_heatmap, dpi=150, bbox_inches="tight")
plt.close()
logp(f"Heatmap saved to {out_heatmap}")

# ════════════════════════════════════════════════════════════════════
# STEP 5 — PCA plot
# ════════════════════════════════════════════════════════════════════
logp("\nGenerating PCA plot...")

# Use top 500 most variable genes for PCA (standard practice)
top_pca_genes = log_norm.var(axis=1).sort_values(ascending=False).head(500).index
pca_input = log_norm.loc[top_pca_genes].T  # samples x genes

pca = PCA(n_components=2)
pca_coords = pca.fit_transform(pca_input)
var_explained = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(
    pca_coords, columns=["PC1", "PC2"], index=pca_input.index
)
pca_df["condition"] = metadata.loc[pca_df.index, "condition"]

plt.figure(figsize=(8, 7))
for cond, color in zip([numerator, denominator], ["#E63946", "#457B9D"]):
    subset = pca_df[pca_df["condition"] == cond]
    plt.scatter(subset["PC1"], subset["PC2"], label=cond, s=100,
                color=color, edgecolors="black", linewidths=0.8)

for sample_id, row in pca_df.iterrows():
    plt.annotate(sample_id, (row["PC1"], row["PC2"]),
                  fontsize=8, xytext=(5, 5), textcoords="offset points")

plt.xlabel(f"PC1 ({var_explained[0]:.1f}% variance)", fontsize=12)
plt.ylabel(f"PC2 ({var_explained[1]:.1f}% variance)", fontsize=12)
plt.title("PCA of Samples (top 500 variable genes)", fontsize=12)
plt.legend(title="Condition")
plt.tight_layout()
plt.savefig(out_pca, dpi=150)
plt.close()
logp(f"PCA plot saved to {out_pca}")

logp("\n✓ DESeq2 analysis complete.")
log.close()
