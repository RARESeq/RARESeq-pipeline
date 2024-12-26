args = commandArgs(T)
input_file = args[1]
output_file = paste0(args[1], ".with_symbol")
protein_coding_genes = args[2]

gene_mapping = read.delim(protein_coding_genes, header = F)
gene_mapping$gene_id_simplified = unlist(lapply(as.character(gene_mapping$V1), function(x) strsplit(x, "[.]")[[1]][1]))

data = read.delim(input_file)
data$gene_id_simplified = unlist(lapply(as.character(data[,1]), function(x) strsplit(x, "[.]")[[1]][1]))
data = data[data$gene_id_simplified %in% gene_mapping$gene_id_simplified,]
data$gene_symbol = gene_mapping[match(data$gene_id_simplified, gene_mapping$gene_id_simplified),]$V2
data$gene_symbol = ifelse(is.na(data$gene_symbol), as.character(data[,1]), as.character(data$gene_symbol))

write.table(data, output_file, sep = "\t", row.names = F)

