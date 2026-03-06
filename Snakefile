import os

configfile: "config.yaml"

SPECIES = config["species"]
SPINFO = config["species_info"][SPECIES]
CELLTYPES = config["celltypes"]
STRANDS = config["strands"]

rule all:
  input:
    expand("results/{celltype}/{strand}/network_final.tsv", strand=STRANDS, celltype=CELLTYPES),
    expand("results/{celltype}/{strand}/TF_expression.tsv", strand=STRANDS, celltype=CELLTYPES),
    expand("results/{celltype}/network_scored.txt", celltype=CELLTYPES),

rule prepare_data:
  input:
    deg="data/degs/all_degs.tsv",
    docr="data/docrs/all_docrs.tsv",
    gtf=SPINFO["annotation_gtf"],
    orth_human=SPINFO["orthologs_to_human"],
    orth_mouse=SPINFO["orthologs_to_mouse"],
    fasta=SPINFO["genome_fasta"]
  output:
    peaks="results/{celltype}/peaks_within_degs.tsv",
    forward="results/{celltype}/docr_sequences_forward_lin.tsv",
    reverse_seq="results/{celltype}/docr_sequences_reverse_lin.tsv",
    pos_genes="results/{celltype}/positive_genes.tsv",
    neg_genes="results/{celltype}/negative_genes.tsv"
  shell:
    """
    echo 'Running prepare_data for {wildcards.celltype}'
    bash prepare_data.bash {input.deg} {input.docr} {wildcards.celltype} {SPECIES} {input.orth_human} {input.orth_mouse} {input.fasta} {input.gtf}
    """

rule run_pipeline:
  input:
    peaks="results/{celltype}/peaks_within_degs.tsv",
    seqs="results/{celltype}/docr_sequences_{strand}_lin.tsv",
    pos_genes="results/{celltype}/positive_genes.tsv",
    orth_mouse=SPINFO["orthologs_to_mouse"],
    orth_human=SPINFO["orthologs_to_human"]
  output:
    network="results/{celltype}/{strand}/network_final.tsv",
    tf_expr="results/{celltype}/{strand}/TF_expression.tsv",
    gene_expr="results/{celltype}/{strand}/gene_expression.tsv"
  shell:
   """
    echo 'Running run_pipeline for {wildcards.celltype} strand {wildcards.strand}'
    bash run_motif_pipeline.sh {wildcards.strand} {wildcards.celltype} {input.orth_mouse} {input.orth_human}
   """

rule post_processing:
    input:
        fwd_network = "results/{celltype}/forward/network_final.tsv",
        rev_network = "results/{celltype}/reverse/network_final.tsv",
        fwd_gene_expr = "results/{celltype}/forward/gene_expression.tsv",
        rev_gene_expr = "results/{celltype}/reverse/gene_expression.tsv",
        fwd_tf_expr   = "results/{celltype}/forward/TF_expression.tsv",
        rev_tf_expr   = "results/{celltype}/reverse/TF_expression.tsv"
    output:
        scored   = "results/{celltype}/network_scored.txt",
    shell:
        """
        echo 'Executing post-processing for {wildcards.celltype}'
        bash run_post_pipeline.sh {wildcards.celltype}
        """
