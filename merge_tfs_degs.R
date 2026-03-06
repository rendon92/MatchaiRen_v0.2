args <- commandArgs(trailingOnly = TRUE)

DEGs_peaks_path <- args[1]
de_TFs_peaks_path <- args[2]
tf_peaks_degs_outpath <- args[3]
tf_degs_network_outpath <- args[4]

DEGs_peaks <- read.delim(DEGs_peaks_path, header=FALSE)
de_TFs_peaks <- read.delim(de_TFs_peaks_path, header=FALSE)
network <- merge(de_TFs_peaks, DEGs_peaks, by = "V1")
write.table(network, tf_peaks_degs_outpath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
network <- network[, c(2,3)]
write.table(network, tf_degs_network_outpath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
