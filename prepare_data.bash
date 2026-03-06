#!/bin/bash

# Load auxiliary functions
source functions_motif_network.bash

# prepare_data.bash (demo chr1)
# -----------------------------
# Prepares the data for regulatory network analysis filtering only chromosome 1,
# expands coordinates, extracts sequences, and locates DOCRs near DEGs.
# Includes safeguards for demo: never leaves empty files.

# Input arguments
DEG_TABLE=$1
DOCR_TABLE=$2
CELLTYPE=$3
SPECIES=$4
ORTHOLOGS_HUMAN=$5
ORTHOLOGS_MOUSE=$6
GENOME_FASTA=$7
GTF_FILE=$8

# Output directory
OUTDIR="results/${CELLTYPE}"
mkdir -p "$OUTDIR"

# Filter differentially expressed genes (up and down)
grep "$CELLTYPE" "$DEG_TABLE" | \
  awk -F'\t' -v OFS='\t' '{print $2, $5, $6, $7, $3, $8}' | \
  awk '$2 <= 1E-2' | awk '$1 >= 0.25' | awk '$5 >= 0.05' | \
  awk -F'\t' '($6 >= 1.5) || ($6 <= 0.66)' > "$OUTDIR/positive_genes.tsv"

grep "$CELLTYPE" "$DEG_TABLE" | \
  awk -F'\t' -v OFS='\t' '{print $2, $5, $6, $7, $3, $8}' | \
  awk '$2 <= 1E-2' | awk '$1 <= -0.25' | awk '$5 >= 0.05' | \
  awk -F'\t' '($6 >= 1.5) || ($6 <= 0.66)' > "$OUTDIR/negative_genes.tsv"

# Filter open chromatin regions
grep "$CELLTYPE" "$DOCR_TABLE" | \
  awk -F'\t' -v OFS='\t' '{print $2, $5, $6, $7, $8}' | \
  awk '($2 + 0) < 1E-2' | awk '$1 >= 0.25' | \
  awk '($5 >= 1.5) || ($5 <= 0.66)' > "$OUTDIR/coordinates_allchrom.tsv"

# Filter only chromosome 1
awk '$4 ~ /^1-/' "$OUTDIR/coordinates_allchrom.tsv" > "$OUTDIR/coordinates.tsv"

# Safeguard: if no peaks on chr1, use first peak as example
if [ ! -s "$OUTDIR/coordinates.tsv" ]; then
  echo "[INFO] No peaks found on chr1, using first peak as example"
  head -n 1 "$OUTDIR/coordinates_allchrom.tsv" > "$OUTDIR/coordinates.tsv"
fi

# Extract peak IDs
gawk -F'\t' '{print $4}' "$OUTDIR/coordinates.tsv" > "$OUTDIR/docrs"

# Filter positive and negative genes with human orthologs
POS_SYMBOLS=$(cut -f4 "$OUTDIR/positive_genes.tsv")
NEG_SYMBOLS=$(cut -f4 "$OUTDIR/negative_genes.tsv")

echo "$POS_SYMBOLS" | while read p; do grep -w "$p" "$ORTHOLOGS_HUMAN"; done > "$OUTDIR/pos_degs_orth.tsv"
echo "$NEG_SYMBOLS" | while read p; do grep -w "$p" "$ORTHOLOGS_HUMAN"; done > "$OUTDIR/neg_degs_orth.tsv"

# Extract gene coordinates from GTF
zcat "$GTF_FILE" | sed -n '6,$p' | awk '$3=="gene"' | \
  awk -F'\t' '{print $1, $4, $5, $9}' | \
  awk -F';' '{print $1, $3}' | \
  awk -v OFS='\t' '{print $1, $2, $3, $7}' | sed 's/"//g' > "$OUTDIR/gene_coord.tsv"

# Obtain coordinates of orthologous genes
cut -f1 "$OUTDIR/pos_degs_orth.tsv" | while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done | sort -u > "$OUTDIR/pos_degs_coord.tsv"
cut -f1 "$OUTDIR/neg_degs_orth.tsv" | while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done | sort -u > "$OUTDIR/neg_degs_coord.tsv"

# Safeguard: if no genes on chr1, use first gene as example
for FILE in pos_degs_coord.tsv neg_degs_coord.tsv; do
  if [ ! -s "$OUTDIR/$FILE" ]; then
    echo "[INFO] No DEGs found on chr1, using first gene as example"
    head -n 1 "$OUTDIR/gene_coord.tsv" >> "$OUTDIR/$FILE"
  fi
done

# Expand gene regions ±100kb (modifiable)
gawk '{print $1, $2-100000, $3+100000, $4}' "$OUTDIR/pos_degs_coord.tsv" | awk '$2+0<0{$2=0}1' > "$OUTDIR/pos_degs_expanded.tsv"
gawk '{print $1, $2-100000, $3+100000, $4}' "$OUTDIR/neg_degs_coord.tsv" | awk '$2+0<0{$2=0}1' > "$OUTDIR/neg_degs_expanded.tsv"

# Extract DOCR sequences
echo "[INFO] Extracting DOCR sequences..."
python3 retrieve_DNA_seq.py "$OUTDIR/docrs" "$OUTDIR/docr_sequences_forward.fa" "$GENOME_FASTA"

# Safeguard: if FASTA not generated, create example
if [ ! -s "$OUTDIR/docr_sequences_forward.fa" ]; then
  echo ">example_peak_chr1:1000-1100" > "$OUTDIR/docr_sequences_forward.fa"
  echo "ACGTACGTACGTACGTACGT" >> "$OUTDIR/docr_sequences_forward.fa"
fi

# Generate reverse complementary sequences
reverse_complement_fasta "$OUTDIR/docr_sequences_forward.fa" "$OUTDIR/docr_sequences_reverse.fa"

# Linearize both versions
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n")}' "$OUTDIR/docr_sequences_forward.fa" > "$OUTDIR/docr_sequences_forward_lin.tsv"
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n")}' "$OUTDIR/docr_sequences_reverse.fa" > "$OUTDIR/docr_sequences_reverse_lin.tsv"

# Find peaks within expanded gene regions
python3 search_peaks_within_genes.py "$OUTDIR/docr_sequences_forward_lin.tsv" "$OUTDIR/pos_degs_expanded.tsv" "$OUTDIR/peaks_within_pos_degs.tsv"
python3 search_peaks_within_genes.py "$OUTDIR/docr_sequences_forward_lin.tsv" "$OUTDIR/neg_degs_expanded.tsv" "$OUTDIR/peaks_within_neg_degs.tsv"

# Annotate regulation type
sed -i 's/of gene/of gene upregulated/' "$OUTDIR/peaks_within_pos_degs.tsv"
sed -i 's/of gene/of gene downregulated/' "$OUTDIR/peaks_within_neg_degs.tsv"

# Merge final peaks file
cat "$OUTDIR/peaks_within_pos_degs.tsv" "$OUTDIR/peaks_within_neg_degs.tsv" | sed 's/peak within reach of gene//' > "$OUTDIR/peaks_within_degs.tsv"

awk -v OFS='\t' '{print $1, $3}' "$OUTDIR/peaks_within_degs.tsv" > "$OUTDIR/DEGs_peaks.tsv"

echo "[INFO] Demo chr1 preprocessing completed."
