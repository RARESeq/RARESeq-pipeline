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
	files = list.files(inp, pattern = "*.molecule_size.freq$")
	if(length(files) == 0)
	{
		next
	}
	data = read.delim(file.path(inp, files[1]))
	data$ID = gsub("Sample_", "", strsplit(files, "__")[[1]][1])
	data$Density = data$Count / sum(data$Count)
	all_data = rbind(all_data, data)
}


pdf(file.path(output_dir, "molecule_length_histograms.pdf"), width = 4, height = 4, family = "Helvetica", useDingbats = F)
splits = split(all_data, all_data$ID)
lapply(splits, function(spl){
g <- ggplot(spl, aes(x = Molecule_size, y = Density)) +
	geom_bar(stat = "identity") + 
	theme_bw() + 
	theme(aspect.ratio = 0.5) + 
	theme(panel.grid = element_blank()) + 
	ggtitle(gsub(".", "-", as.character(spl$ID[1]), fixed = T)) + 
	ylab("Density") + 
	xlab("Molecule length") 
plot(g)
})
dev.off()

all_data = NULL
for(folder in folders)
{
	inp = file.path(input_dir, folder)
	files = list.files(inp, pattern = "*.molecule_size.summary$")
	if(length(files) == 0)
	{
		next
	}
	data = read.delim(file.path(inp, files[1]))
	data$ID = gsub("Sample_", "", strsplit(files, "__")[[1]][1])
	data$Total = data$Overlapping + data$Non_Overlapping 
	data$Overlapping.perc = (data$Overlapping - data$Inconsistent) / data$Total
	data$Non_Overlapping.perc = data$Non_Overlapping / data$Total
	data$Inconsistent.perc = data$Inconsistent / data$Total
	all_data = rbind(all_data, data)
}

melted_data = melt(all_data[,c("ID", "Overlapping.perc", "Non_Overlapping.perc", "Inconsistent.perc")], id.vars = "ID")
melted_data$variable = gsub(".perc", "", as.character(melted_data$variable))
melted_data$variable = gsub("_", " ", as.character(melted_data$variable))
melted_data$variable = gsub("Non Overlapping", "Non-overlapping", as.character(melted_data$variable))

pdf(file.path(output_dir, "molecule_length_barplots.pdf"), width = 7, height = 3.5, family = "Helvetica", useDingbats = F)
g <- ggplot(melted_data, aes(x = ID, y = value, fill = variable)) +
	geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
	theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
		legend.position = "none") + 
	ylab("Molecules (%)") + 
	xlab("") + 
	scale_fill_manual(values = brewer.pal(8, "Set1"), name = "Overlap type") + 
	theme_bw() + 
	theme(aspect.ratio = 0.5) + 
	theme(panel.grid = element_blank()) + 
	theme(axis.text.x= element_text(angle = 45, hjust = 1, vjust = 1))
plot(g)
dev.off()
