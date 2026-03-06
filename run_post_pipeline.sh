#!/bin/bash
# Script post processing: combine forward/reverse results and score nodes

set -euo pipefail

source functions_motif_network.bash

CELLTYPE=$1
OUTDIR="results/${CELLTYPE}"
ANNOTATED_GENES_FILE="/data/scRNAseq/seurat/zf_annot/wt_expression_fish.txt"

# 1. Merge networks and gene expression
merge_network_strands \
  "$OUTDIR/forward/network_final.tsv" \
  "$OUTDIR/reverse/network_final.tsv" \
  "$OUTDIR/network_final.tsv"

merge_expression_files \
  "$OUTDIR/forward/gene_expression.tsv" \
  "$OUTDIR/reverse/gene_expression.tsv" \
  "$OUTDIR/forward/TF_expression.tsv" \
  "$OUTDIR/reverse/TF_expression.tsv" \
  "$OUTDIR/gene_expression.tsv" \
  "$OUTDIR/TF_expression.tsv" \
  "$OUTDIR/expression.tsv"

# 2. Calculate node degree
calculate_degrees \
  "$OUTDIR/network_final.tsv" \
  "$OUTDIR/degree_TFs.txt" \
  "$OUTDIR/degree_targets.txt" \
  "$OUTDIR/degrees.txt" \
  "$OUTDIR/degrees_sorted.txt"

# 3. Join gene expression and node degree
join_expression_and_degrees \
  "$OUTDIR/degrees_sorted.txt" \
  "$OUTDIR/expression.tsv" \
  "$OUTDIR/expression_sorted.tsv" \
  "$OUTDIR/merged.tsv"

# 4. Scoring y filtrado de red
score_network \
  "$OUTDIR/network_final.tsv" \
  "$OUTDIR/merged.tsv" \
  "$OUTDIR/network_scored.txt" \
#  "$OUTDIR/filtered_network.txt"

echo "Post-processing completed for $CELLTYPE"
