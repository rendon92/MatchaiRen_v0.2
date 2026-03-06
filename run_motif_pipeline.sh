
#!/bin/bash

# run_motif_pipeline.sh
# ----------------------
# Wrapper script to execute some modular functions for HOMER analysis. Can support forward, reverse and all strands


set -euo pipefail
#set -xe
#set -euxo pipefail

source functions_motif_network.bash 

strand=$1
celltype=$2
ORTH_MOUSE=$3
ORTH_HUMAN=$4

if [[ -z "$strand" || -z "$celltype" ]]; then
  echo "Uso: $0 <forward|reverse|all> <celltype>"
  exit 1
fi

run_strand() {
  local strand_name=$1
  local SEQ_FILE="results/${celltype}/docr_sequences_${strand_name}_lin.tsv"
  local PEAK_LIST="results/${celltype}/peaks_within_degs.tsv"
  local MOTIF_DB="data/motifs/jaspar2022.motifs"
  local OUTDIR="results/${celltype}/${strand_name}"

  mkdir -p "$OUTDIR"

  # Linearize peaks
  linearize_peak_sequences "$PEAK_LIST" "$SEQ_FILE" "$OUTDIR/peaks_degs_homer_linearized.tsv"

  # Process each peak


mapfile -t peaks < "$OUTDIR/peaks_degs_homer_linearized.tsv"
total=${#peaks[@]}
count=0
bar_width=40

echo "→ Processing HOMER for strand [$strand_name] with $total entries"

# Initial time
last_update=$(date +%s)

for peak_entry in "${peaks[@]}"; do
  echo "$peak_entry" | tr '\t' '\n' > "$OUTDIR/helper.fa"
  coord=$(echo "$peak_entry" | awk '{print $1}')

  generate_background_fasta "$OUTDIR/helper.fa" "$OUTDIR/bg.fa"
  run_homer_on_sequence "$OUTDIR/helper.fa" "$OUTDIR/bg.fa" "$MOTIF_DB" "$OUTDIR/homer_${coord}.tsv"
  filter_homer_results "$OUTDIR/homer_${coord}.tsv" "$coord" "$OUTDIR/motifs_filtered.tsv"
  awk -F':::|\t' '{print $2 "\t" $4}' "$OUTDIR/motifs_filtered.tsv" > "$OUTDIR/putative_interactions"

  count=$((count + 1))
  percent=$((count * 100 / total))
  filled=$((percent * bar_width / 100))
  empty=$((bar_width - filled))
  bar=$(printf "%${filled}s" | tr ' ' '#')$(printf "%${empty}s" | tr ' ' '-')

  now=$(date +%s)
  if (( now - last_update >= 5 || count == total )); then
    printf "[%s] [%s] %3d%% (%d/%d)\n" "$strand_name" "$bar" "$percent" "$count" "$total"
    last_update=$now
  fi
done

echo -e " Processing of HOMER finished for strand $strand_name"


    if [[ ! -s "$OUTDIR/putative_interactions" ]]; then
      echo "[✘] ERROR: putative_interactions was not correctly generated."
      exit 1
    fi

    # Analize orthologs using parse_genes.R
    parse_orthologs "$OUTDIR/putative_interactions" "$ORTH_MOUSE" "$ORTH_HUMAN" "$OUTDIR/human_mouse_query"


  # Orthologs and network

  cut -f4 "results/${celltype}/positive_genes.tsv" > "$OUTDIR/pos_genes"
  cut -f4 "results/${celltype}/negative_genes.tsv" > "$OUTDIR/neg_genes"

  grep -wf "$OUTDIR/pos_genes" "$OUTDIR/human_mouse_query" | sed 's/>//' > "$OUTDIR/de_TFs_peaks"


    merge_tfs_and_degs \
      "results/${celltype}/DEGs_peaks.tsv" \
      "$OUTDIR/de_TFs_peaks" \
      "$OUTDIR/tf_peaks_degs_network.tsv" \
      "$OUTDIR/tf_degs_network.tsv" \
      "$OUTDIR/tf_degs_network_interaction.tsv"

    classify_interactions \
      "$OUTDIR/pos_genes" \
      "$OUTDIR/neg_genes" \
      "$OUTDIR/tf_degs_network_interaction.tsv" \
      "$OUTDIR/network_activated.tsv" \
      "$OUTDIR/network_repressed.tsv" \
      "$OUTDIR/network_final.tsv"

    collect_tf_expression \
      "$OUTDIR/network_final.tsv" \
      "results/${celltype}/positive_genes.tsv" \
      "results/${celltype}/negative_genes.tsv" \
      "$OUTDIR/gene_expression.tsv" \
      "$OUTDIR/TF.tsv" \
      "$OUTDIR/TF_expression.tsv"
    echo "Cleaning intermediate files in $OUTDIR"
    rm -f "$OUTDIR"/homer_*.tsv "$OUTDIR"/helper.fa "$OUTDIR"/bg.fa
  echo "Pipeline finished for strand: $strand_name"
}



if [[ "$strand" == "all" ]]; then
  run_strand forward
  run_strand reverse
else
  run_strand "$strand"
fi

