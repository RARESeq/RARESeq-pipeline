library(ggplot2)
library(RColorBrewer)
library(reshape2)   

args = commandArgs(T) 
input_dir = args[1]
output_dir = input_dir

folders = list.files(input_dir)
all_stats = NULL
for(folder in folders)
{
	inp = file.path(input_dir, folder)
	metrics_file = list.files(inp, pattern = "*.sorted.bam.metrics.tsv$")
	if(length(metrics_file) == 0)
	{
		next
	}
	metrics = read.delim(file.path(inp, metrics_file[1]))
	colnames(metrics)[2] = strsplit(colnames(metrics)[2], "__")[[1]][1]
	if(is.null(all_stats))
	{
		all_stats = metrics
	}else{
		all_stats = merge(all_stats, metrics, by = "Sample")
	}
}
colnames(all_stats)[1] = "Type"
write.table(all_stats, file.path(output_dir, "combined_RNASeQC_metrics.tsv"), sep = "\t", row.names = F)


melted = melt(all_stats[all_stats$Type %in% c("End 1 Sense Rate", "End 2 Sense Rate"),], id.vars = "Type")
melted$variable = gsub("Sample_", "", as.character(melted$variable))
melted$Type = gsub(" Rate", "", as.character(melted$Type))
pdf(file.path(output_dir, "strandness.pdf"), width = 6, height = 5 + 0.05 * length(unique(melted$variable)), family = "Helvetica", useDingbats = F)
g <- ggplot(melted, aes(y = value, x = variable, fill = Type)) + 
		geom_bar(stat = "identity", position = "dodge", width = 0.7) +
		scale_fill_manual(values = brewer.pal(8, "Set1")) + 
		theme_bw() + 
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 2 + 0.001 * length(unique(melted$variable))) + 
		#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
		xlab("") + 
		ylab("Sense rate") +
		coord_flip()
plot(g)
dev.off()


melted = melt(all_stats[all_stats$Type %in% c("rRNA Rate"),], id.vars = "Type")
melted$variable = gsub("Sample_", "", as.character(melted$variable))
melted$value = melted$value * 100
pdf(file.path(output_dir, "rRNA.pdf"), width = 6, height = 5 + 0.05 * length(unique(melted$variable)), family = "Helvetica", useDingbats = F)
g <- ggplot(melted, aes(y = value, x = variable)) + 
		geom_bar(stat = "identity", position = "dodge", width = 0.7) +
		scale_fill_manual(values = brewer.pal(8, "Set1")) + 
		theme_bw() + 
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 2 + 0.001 * length(unique(melted$variable))) + 
		#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
		xlab("") + 
		ylab("rRNA rate (%)") +
		coord_flip()
plot(g)
dev.off()

melted = melt(all_stats[all_stats$Type %in% c("Exonic Rate", "Intronic Rate", "Intergenic Rate", "Ambiguous Alignment Rate"),], id.vars = "Type")
melted$variable = gsub("Sample_", "", as.character(melted$variable))
melted$value = melted$value * 100
melted$Type = gsub(" Rate", "", as.character(melted$Type))
pdf(file.path(output_dir, "read_regions.pdf"), width = 6, height = 5 + 0.05 * length(unique(melted$variable)), family = "Helvetica", useDingbats = F)
g <- ggplot(melted, aes(y = value, x = variable, fill = Type)) + 
		geom_bar(stat = "identity", position = "dodge", width = 0.85) +
		scale_fill_manual(values = brewer.pal(8, "Set1")) + 
		theme_bw() + 
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 2 + 0.001 * length(unique(melted$variable))) + 
		#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
		xlab("") + 
		ylab("Reads (%)") +
		coord_flip()
plot(g)
dev.off()
  