#!/bin/bash
# ════════════════════════════════════════════════════════════════════
# setup_references.sh
# Run this ONCE before the first pipeline run.
# Downloads: HISAT2 pre-built hg38 genome index + Ensembl GTF annotation
# Usage: bash setup_references.sh
# ════════════════════════════════════════════════════════════════════

set -e  # exit on any error

echo "═══════════════════════════════════════════════════════════"
echo " SnakeSEQ — Downloading reference files"
echo " This will download ~4-5 GB. Ensure good internet connection."
echo "═══════════════════════════════════════════════════════════"

mkdir -p resources/hg38_index
cd resources/hg38_index

if [ ! -f genome.1.ht2 ]; then
    echo "[1/2] Downloading HISAT2 hg38 genome index..."
    curl -L -o hg38_genome.tar.gz \
        https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
    tar -xzf hg38_genome.tar.gz
    rm hg38_genome.tar.gz
    echo "      Done. Index files: $(ls genome.*.ht2 | wc -l) parts"
else
    echo "[1/2] HISAT2 index already present — skipping."
fi

cd ../..

if [ ! -f resources/Homo_sapiens.GRCh38.111.gtf ]; then
    echo "[2/2] Downloading Ensembl GRCh38 GTF annotation..."
    curl -L -o resources/Homo_sapiens.GRCh38.111.gtf.gz \
        https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
    gunzip resources/Homo_sapiens.GRCh38.111.gtf.gz
    echo "      Done."
else
    echo "[2/2] GTF annotation already present — skipping."
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo " ✓ All reference files ready."
echo " Next step: snakemake --use-conda -n --cores 4   (dry run)"
echo "═══════════════════════════════════════════════════════════"
