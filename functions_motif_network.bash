#!/bin/bash

# functions_motif_network.bash
# -----------------------------
# Contains the functions necessary to process enriched motifs, TFs and genes regulated in the context of each inferred GRN.


# Linearize DNA sequences associated to chromatin peaks near DEGs
function linearize_peak_sequences() {
  local peak_list=$1
  local all_seqs=$2
  local output=$3

  awk '{print $1}' "$peak_list" | sort -u | while read peak_id; do
    grep "$peak_id" "$all_seqs"
  done > "$output"
}


# Extract reverse complement sequence
function reverse_complement_fasta() {
  local input_fasta=$1
  local output_fasta=$2

  awk 'BEGIN{ORS=""} 
    /^>/ {
      if (seq) { print rc(seq) "\n" };
      print $0 "\n"; seq=""; next
    } 
    {seq=seq $0} 
    END {print rc(seq) "\n"} 
    function rc(s,    i, out) {
      out=""
      for(i=length(s); i>0; i--) {
        c=substr(s,i,1)
        out=out complement[c]
      }
      return out
    } 
    BEGIN{
      complement["A"]="T"; complement["T"]="A";
      complement["C"]="G"; complement["G"]="C";
      complement["a"]="t"; complement["t"]="a";
      complement["c"]="g"; complement["g"]="c"
    }' "$input_fasta" > "$output_fasta"
}


# Executes HOMER for an individual sequence using a background sequence
function run_homer_on_sequence() {
  local seq_file=$1
  local background_fasta=$2
  local motif_db=$3
  local output_file=$4

  homer2 known -i "$seq_file" -b "$background_fasta" -p 8 -m "$motif_db" > "$output_file" 2>/dev/null
}


# Create a background sequence with the same length and nucleotide composition as the query DNA peak
function generate_background_fasta() {
  local query_fasta=$1
  local output_fasta=$2
  python3 create_background_fasta_same_freqs.py "$output_fasta" "$query_fasta"
}


# Filter enriched motifs by score and p-value, associate to DNA peaks
function filter_homer_results() {
  local homer_output=$1
  local peak_coord=$2
  local output_filtered=$3

  sed '/^Motif/d' "$homer_output" | sed 's/%//g' | \
    awk -v peak="$peak_coord" -F'\t' -v OFS='\t' '$6==1 && $9<=10 {print $1, $2, peak}' >> "$output_filtered"
}


# Execute R script to parse orthologs
function parse_orthologs() {
  local interactions_file=$1
  local mouse_orth=$2
  local human_orth=$3
  local human_mouse_query=$4
  Rscript parse_genes.R "$interactions_file" "$mouse_orth" "$human_orth" "$human_mouse_query"
}


# Execute R script to merge TFs and DEGs in a regulatory network
function merge_tfs_and_degs() {

  local DEGs_peaks_path=$1
  local de_TFs_peaks_path=$2
  local tf_peaks_degs_outpath=$3
  local tf_degs_network_outpath=$4
  local tf_degs_network_interaction=$5

  Rscript merge_tfs_degs.R "$DEGs_peaks_path" "$de_TFs_peaks_path" "$tf_peaks_degs_outpath" "$tf_degs_network_outpath"
#  awk '$1 = toupper($1)' "$tf_degs_network_outpath" | uniq | sort > "$tf_degs_network_TF_uppercase"
  awk '{print toupper($1) "\tinteracts\t" $2}' "$tf_degs_network_outpath" > "$tf_degs_network_interaction"
}

# Classify interactions as activators or repressors
function classify_interactions() {
  local deg_pos=$1
  local deg_neg=$2
  local network_in=$3
  local output_activates=$4
  local output_represses=$5
  local network_out=$6

  grep -wf "$deg_pos" "$network_in" | sed 's/interacts/activates/' > "$output_activates"
  grep -wf "$deg_neg" "$network_in" | sed 's/interacts/represses/' > "$output_represses"
  cat "$output_activates" "$output_represses" > "$network_out"
}


# Extract expression of regulatory TFs
function collect_tf_expression() {
  local network_file=$1
  local pos_deg_table=$2
  local neg_deg_table=$3
  local gene_expression=$4
  local tfs=$5
  local tf_expression=$6

  
  cat "$pos_deg_table" "$neg_deg_table" | awk -F'\t' -v OFS='\t' '{print $4, $1, $2, $5, $6}' > "$gene_expression"
  awk '{print $1}' "$network_file" | sort -u > "$tfs"
  while read tf; do
    grep -i "$tf" "$gene_expression"
  done < "$tfs" | awk '{ $1 = toupper($1); print }' > "$tf_expression"
}

# Merge forward and reverse networks
function merge_network_strands() {
  local forward=$1
  local reverse=$2
  local output=$3

  cat "$forward" "$reverse" | sort | uniq | sort > "$output"
}

# Combine expression files
function merge_expression_files() {
  local fw_gene_expr=$1
  local rv_gene_expr=$2
  local fw_tf_expr=$3
  local rv_tf_expr=$4
  local gene_expr_output=$5
  local tf_expr_output=$6
  local expr_output=$7

  cat "$fw_gene_expr" "$rv_gene_expr" | sort | uniq | sort > "$gene_expr_output"
  cat "$fw_tf_expr" "$rv_tf_expr" | sort | uniq | sort > "$tf_expr_output"
  #cat "$gene_expr_output" "$tf_expr_output" > "$expr_output"
  (cat "$gene_expr_output"; awk '{$1=$1}1' OFS='\t' "$tf_expr_output") > "$expr_output"
}

# Calculate TF and target degrees
function calculate_degrees() {
  local network_file=$1
  local out_TFs=$2
  local out_targets=$3
  local out_all=$4
  local out_sorted=$5

  awk '{print $1}' "$network_file" | sort | uniq -c | sort -nr | sed 's/^\s\+//g' > "$out_TFs"
  awk '{print $3}' "$network_file" | sort | uniq -c | sort -nr | sed 's/^\s\+//g' > "$out_targets"
  cat "$out_TFs" "$out_targets" | sort -nr > "$out_all"
  awk '{print $2"\t"$1}' "$out_all" | sort -k1,1 > "$out_sorted"
}

# Join expression and degrees
function join_expression_and_degrees() {
  local degrees_sorted=$1
  local expression=$2
  local expression_sorted=$3
  local output=$4

  sort -k1,1 "$expression" > "$expression_sorted"
  join -t $'\t' -a 2 -e "NA" -o 2.1 1.2 2.2 2.3 2.4 2.5 -1 1 -2 1 "$degrees_sorted" "$expression_sorted" > "$output"
}


# R script for Network scoring and filtering
function score_network() {
  local network_file=$1
  local merged_file=$2
  local scored_out=$3
#  local filtered_out=$4

  Rscript network_scoring.R "$network_file" "$merged_file" "$scored_out"  #"$filtered_out"
}




