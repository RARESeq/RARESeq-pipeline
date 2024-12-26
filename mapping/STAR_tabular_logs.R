args = commandArgs(T)
input_dir = args[1]
output_dir = dirname(input_dir)
dir.create(output_dir, recursive = T, showWarning = F)

files = sort(list.files(input_dir, pattern = "*Log.final.out"))
files = files[grepl("Sample_", files)]
all_data = data.frame()
for(file in files)
{
	
	id = strsplit(basename(file), ".Log.", fixed = T)[[1]][1]
			
	data = data.frame()
	con  <- file(file.path(input_dir, file), open = "r")
	#skip first 5 lines
	readLines(con, n = 5, warn = FALSE)
	while (length(l <- readLines(con, n = 1, warn = FALSE)) > 0) {
		l = trimws(l)
		spl = strsplit(l, "\t")[[1]]
		if(length(spl) == 1){
			tmp = cbind(spl, NA)
		}else{
			spl[1] = gsub(" |", "", spl[1], fixed = T)
			spl[2] = gsub("%", "", spl[2], fixed = T)
			tmp = cbind(spl[1], spl[2])
		}
		if(nrow(data) == 0)
		{
			data = tmp
		}else{
			data = rbind(data, tmp)
		}
	} 
	close(con)
	colnames(data) = c("Descr", id)
	if(nrow(all_data) == 0)
	{
		all_data = data
	}else{
		all_data = merge(all_data, data, by = "Descr")
	}
}

all_data = as.data.frame(all_data)
all_data = all_data[match(data[,1], as.data.frame(all_data)$Descr),]
write.table(all_data, file.path(output_dir, "tabular_logs.txt"), row.names = F, sep = "\t", quote = F)

s = do.call(cbind, strsplit(c("__", colnames(all_data[-1])), "__"))
colnames(s) = colnames(all_data)
all_data = rbind(s, all_data)
all_data = apply(all_data, 2, function(x) ifelse(is.na(x), "", x))
write.table(all_data, file.path(output_dir, "tabular_logs_pretty.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
