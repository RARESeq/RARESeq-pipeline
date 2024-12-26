library(reshape2)
args = commandArgs(T)
input_dir = args[1]
output_dir = args[2]
dir.create(output_dir, recursive = T, showWarning = F)

files = list.files(input_dir, pattern = ".genes.results.with_symbol$")

all_data = data.frame()
all_conditions = data.frame()
for(file in files)
{
	input_path = file.path(input_dir, file)
	print(input_path)
	data = read.delim(input_path)
	data = data[!duplicated(data$gene_symbol),]
	s = strsplit(file, ".genes")[[1]][1]
	data$ID = s
	
	data = data[,c('ID', "gene_symbol", "TPM")]

	if(nrow(all_data) == 0)
	{
		all_data = data
	}else{
		all_data = rbind(all_data, data)
	}
}
new_data = dcast(all_data, gene_symbol~ID, value.var = "TPM")
new_data[,-1] = apply(new_data[,-1], 2, function(x) x * 1e6 / sum(x))
write.table(new_data, file.path(output_dir, "expression_matrix_protein_coding.txt"), sep = "\t", row.names = F, quote = F)

