# network_scoring.R
score_network <- function(network_file, merged_file, scored_out) {

  network <- read.table(network_file, quote="\"", comment.char="")
  network_table <- read.delim(merged_file, header=FALSE)
  colnames(network_table) <- c("gene_name", "Degree", "log2FC", "padj", "pct", "pct.dif")
  network_table <- na.omit(network_table)

  # Normalizaciones
  network_table$Degree_norm_0_100 <- (
    (network_table$Degree - min(network_table$Degree)) /
    (max(network_table$Degree) - min(network_table$Degree))
  ) * 100

  network_table$pct.dif_norm_100 <- (
    (network_table$pct.dif - min(network_table$pct.dif)) /
    (max(network_table$pct.dif) - min(network_table$pct.dif))
  ) * 100

  network_table$pct_norm_100 <- (
    (network_table$pct - min(network_table$pct)) /
    (max(network_table$pct - min(network_table$pct)))
  ) * 100

  network_table$padj_log <- -log10(network_table$padj + .Machine$double.eps)

  network_table$padj_log_norm_100 <- (
    (network_table$padj_log - min(network_table$padj_log, na.rm = TRUE)) /
    (max(network_table$padj_log, na.rm = TRUE) - min(network_table$padj_log, na.rm = TRUE))
  ) * 100

  network_table$log2FC_abs <- abs(network_table$log2FC)

  network_table$log2FC_abs_norm_100 <- (
    (network_table$log2FC_abs - min(network_table$log2FC_abs, na.rm = TRUE)) /
    (max(network_table$log2FC_abs, na.rm = TRUE) - min(network_table$log2FC_abs, na.rm = TRUE))
  ) * 100

  network_table$score_norm <- network_table$log2FC_abs_norm_100*0.5/4 +
                              network_table$padj_log_norm_100*0.5/4 +
                              network_table$Degree_norm_0_100*0.5 +
                              network_table$pct.dif_norm_100*0.5/4 +
                              network_table$pct_norm_100*0.5/4

  TFs <- unique(network$V1)
#  network_table$TF <- ifelse(network_table$gene_name %in% TFs, "YES", "NO")

  network_table$TF <- NA
  network_table$TF[network_table$gene_name %in% TFs] <- "YES"
  network_table$TF[is.na(network_table$TF)] <- "NO"





  write.table(network_table, scored_out, sep = "\t", quote = FALSE, row.names = FALSE)

#  disc_nodes <- subset(network_table, Degree < degree_threshold)$gene_name
#  filtered_network <- network[!(network$V3 %in% disc_nodes), ]
#  write.table(filtered_network, filtered_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Parse command line arguments when run as script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  score_network(
    network_file = args[1],
    merged_file = args[2],
    scored_out = args[3]
#    filtered_out = args[4]
  )
}
