#!/bin/bash

# Cargar funciones auxiliares
source functions_motif_network.bash

# prepare_data.bash
# ------------------
# Prepara los datos para análisis de red reguladora: filtra DEGs y DOCRs,
# expande coordenadas, extrae secuencias, y localiza regiones reguladoras cercanas a genes DEGs.

# Entrada
DEG_TABLE=$1
DOCR_TABLE=$2
CELLTYPE=$3
SPECIES=$4
ORTHOLOGS_HUMAN=$5
ORTHOLOGS_MOUSE=$6
GENOME_FASTA=$7
GTF_FILE=$8

# Output dir
OUTDIR="results/${CELLTYPE}"
mkdir -p "$OUTDIR"

# Filtrar genes diferencialmente expresados (up y down)
grep "$CELLTYPE" "$DEG_TABLE" | \
  awk -F'	' -v OFS='\t' '{print $2, $5, $6, $7, $3, $8}' | \
  awk '$2 <= 1E-2' | awk '$1 >= 0.75' | awk '$5 >= 0.05' | \
  awk -F'\t' '($6 >= 1.5) || ($6 <= 0.66)' > "$OUTDIR/positive_genes.tsv"

grep "$CELLTYPE" "$DEG_TABLE" | \
  awk -F'	' -v OFS='\t' '{print $2, $5, $6, $7, $3, $8}' | \
  awk '$2 <= 1E-2' | awk '$1 <= -0.75' | awk '$5 >= 0.05' | \
  awk -F'\t' '($6 >= 1.5) || ($6 <= 0.66)' > "$OUTDIR/negative_genes.tsv"

# Filtrar regiones de cromatina abierta
grep "$CELLTYPE" "$DOCR_TABLE" | \
  awk -F'\t' -v OFS='\t' '{print $2, $5, $6, $7, $8}' | \
  awk '($2 + 0)  <  1E-2' | awk '$1 >= 0.25' | \
  awk '($5 >= 1.5) || ($5 <= 0.66)' | \
  awk '$5 ~ /^1-/' > "$OUTDIR/coordinates.tsv"

# Extraer IDs de picos
gawk -F'\t' '{print $4}' "$OUTDIR/coordinates.tsv" > "$OUTDIR/docrs"

# Filtrar genes positivos y negativos con ortología a humano
POS_SYMBOLS=$(cut -f4 "$OUTDIR/positive_genes.tsv")
NEG_SYMBOLS=$(cut -f4 "$OUTDIR/negative_genes.tsv")

echo "$POS_SYMBOLS" | while read p; do grep -w "$p" "$ORTHOLOGS_HUMAN"; done > "$OUTDIR/pos_degs_orth.tsv"
echo "$NEG_SYMBOLS" | while read p; do grep -w "$p" "$ORTHOLOGS_HUMAN"; done > "$OUTDIR/neg_degs_orth.tsv"

# Extraer coordenadas génicas del GTF
zcat "$GTF_FILE" | sed -n '6,$p' | awk '$3=="gene"' | \
  awk -F'\t' '{print $1, $4, $5, $9}' | \
  awk -F';' '{print $1, $3}' | \
  awk -v OFS='\t' '{print $1, $2, $3, $7}' | sed 's/"//g' > "$OUTDIR/gene_coord.tsv"

# Obtener coordenadas de genes ortólogos
#while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done < "$OUTDIR/pos_degs_orth.tsv" | sort -u > "$OUTDIR/pos_degs_coord.tsv"
#while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done < "$OUTDIR/neg_degs_orth.tsv" | sort -u > "$OUTDIR/neg_degs_coord.tsv"

cut -f1 "$OUTDIR/pos_degs_orth.tsv" | while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done | sort -u > "$OUTDIR/pos_degs_coord.tsv"
cut -f1 "$OUTDIR/neg_degs_orth.tsv" | while read p; do grep -w "$p" "$OUTDIR/gene_coord.tsv"; done | sort -u > "$OUTDIR/neg_degs_coord.tsv"

awk '$1=="1"' "$OUTDIR/pos_degs_coord.tsv" > "$OUTDIR/pos_degs_coord_chr1.tsv"
mv "$OUTDIR/pos_degs_coord_chr1.tsv" "$OUTDIR/pos_degs_coord.tsv"

awk '$1=="1"' "$OUTDIR/neg_degs_coord.tsv" > "$OUTDIR/neg_degs_coord_chr1.tsv"
mv "$OUTDIR/neg_degs_coord_chr1.tsv" "$OUTDIR/neg_degs_coord.tsv"


# Expandir regiones génicas ±100kb
gawk '{print $1, $2-100000, $3+100000, $4}' "$OUTDIR/pos_degs_coord.tsv" | \
  awk '$2+0<0{$2=0}1' > "$OUTDIR/pos_degs_expanded.tsv"

gawk '{print $1, $2-100000, $3+100000, $4}' "$OUTDIR/neg_degs_coord.tsv" | \
  awk '$2+0<0{$2=0}1' > "$OUTDIR/neg_degs_expanded.tsv"

# Extraer secuencias de DOCRs
echo "[INFO] Extrayendo secuencias de DOCRs..."
python3 retrieve_DNA_seq.py "$OUTDIR/docrs" "$OUTDIR/docr_sequences_forward.fa" "$GENOME_FASTA"

# Generar secuencias reverse complementarias
reverse_complement_fasta "$OUTDIR/docr_sequences_forward.fa" "$OUTDIR/docr_sequences_reverse.fa"

# Linearizar ambas versiones
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n")}' "$OUTDIR/docr_sequences_forward.fa" > "$OUTDIR/docr_sequences_forward_lin.tsv"
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n")}' "$OUTDIR/docr_sequences_reverse.fa" > "$OUTDIR/docr_sequences_reverse_lin.tsv"

# Buscar peaks dentro de regiones expandidas
python3 search_peaks_within_genes.py "$OUTDIR/docr_sequences_forward_lin.tsv" "$OUTDIR/pos_degs_expanded.tsv" "$OUTDIR/peaks_within_pos_degs.tsv"
python3 search_peaks_within_genes.py "$OUTDIR/docr_sequences_forward_lin.tsv" "$OUTDIR/neg_degs_expanded.tsv" "$OUTDIR/peaks_within_neg_degs.tsv"

# Anotar tipo de regulación
sed -i 's/of gene/of gene upregulated/' "$OUTDIR/peaks_within_pos_degs.tsv"
sed -i 's/of gene/of gene downregulated/' "$OUTDIR/peaks_within_neg_degs.tsv"

# Unificar archivo final
cat "$OUTDIR/peaks_within_pos_degs.tsv" "$OUTDIR/peaks_within_neg_degs.tsv" | sed 's/peak within reach of gene//' > "$OUTDIR/peaks_within_degs.tsv"

awk -v OFS='\t' '{print $1, $3}' "$OUTDIR/peaks_within_degs.tsv" > "$OUTDIR/DEGs_peaks.tsv" 

# Fin del preprocesamiento

