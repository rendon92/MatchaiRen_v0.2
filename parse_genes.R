# parse_genes.R

args <- commandArgs(trailingOnly = TRUE)

putative_path <- args[1]
mouse_orth_path <- args[2]
human_orth_path <- args[3]
out_path <- args[4]

# Cargar interacciones candidatas (ya debe existir el fichero)
putative_interactions <- read.delim(putative_path, header=FALSE)

# checking zebra-human correspondence for JASPAR motifs
human_orth <- read.table(human_orth_path, quote="\"", comment.char="")
colnames(human_orth) <- c("V2", "V1")
human_query <- merge(putative_interactions, human_orth, by = "V1")

# checking zebra-mouse correspondence
mouse_orth <- read.table(mouse_orth_path, quote="\"", comment.char="")
colnames(mouse_orth) <- c("V2", "V1")
mouse_query <- merge(putative_interactions, mouse_orth, by = "V1")

# Combinar y exportar
human_mouse_query <- as.data.frame(rbind(human_query, mouse_query))[, c(2,3)]
write.table(human_mouse_query, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
