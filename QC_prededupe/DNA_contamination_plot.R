library(ggplot2)
library(reshape2)
library(RColorBrewer)

args = commandArgs(T)

input_dir = args[1]
output_dir = input_dir

folders = list.files(input_dir)
all_data = NULL
for(folder in folders)
{
	inp = file.path(input_dir, folder)
	files = list.files(inp, pattern = "*.contamination_stats.txt$")
	if(length(files) == 0)
	{
		next
	}
	data = read.delim(file.path(inp, files[1]), header = F)
	data$ID = gsub("Sample_", "", strsplit(files, "__")[[1]][1])
	data$Fraction = data$V2 / sum(data$V2) * 100
	all_data = rbind(all_data, data)
}

all_data = all_data[all_data$V1 == "NotSplit",]
pdf(file.path(output_dir, "DNA_contamination.pdf"), height = 6 + 0.01 * length(unique(all_data$ID)), width = 4, family = "Helvetica", useDingbats = F)
g <- ggplot(all_data, aes(y = Fraction, x = ID)) + 
	geom_bar(stat = "identity", width = .7) +
	theme_bw() + 
	theme(aspect.ratio = 1 + 0.005 * length(unique(all_data$ID))) + 
	ylim(0, 100) + 
	theme(panel.grid = element_blank()) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
	ylab("DNA contamination\n(non-spliced reads) (%)") + 
	xlab("") +
	coord_flip()
plot(g)
dev.off()

