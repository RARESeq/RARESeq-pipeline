library(reshape2)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(T)
input_dir = args[1]
output_dir = args[1]
dir.create(output_dir, recursive = T, showWarning = F)

raw_data = read.delim(file.path(input_dir, "tabular_logs_pretty.txt"), header = F)

data = t(raw_data[,-1,drop = F])
colnames(data) <- raw_data[,1]

nvars = sum(raw_data[,1] == "")

data = data.frame(data[,1:nvars], (apply(data[,(nvars + 1):ncol(data),drop = F], 2, function(x) as.numeric(as.character(x)))), check.names = F)	

data[,"Number of unmapped reads"] = data[,"Number of input reads"] - data[, "Uniquely mapped reads number"] - data[, "Number of reads mapped to multiple loci"]
- data[,"Number of reads mapped to too many loci"]

colnames(data)[1] <- "ID"
melted = melt(data, id.vars = c("ID"))
melted

data = melted
melted = data
plot_variables = c("Number of input reads", "Uniquely mapped reads number", "Number of reads mapped to multiple loci", "Number of unmapped reads")
to_plot = subset(melted, melted$variable %in% plot_variables, droplevels = T)

to_plot$ID = gsub("Sample_", "", as.character(to_plot$ID))
pdf(file.path(output_dir, "read_stats_line.pdf"), width = 6, height = 4 + 0.1 * length(unique(to_plot$ID)), family = "Helvetica", useDingbats = F)
g <- ggplot(to_plot, aes(x = ID, y = value, fill = variable)) +
	geom_bar(stat = "identity", position=position_dodge(), width = 0.75) + 
	theme_bw() + 
	#scale_y_continuous() + #breaks=pretty_breaks(n=6)) +  
	theme(
		aspect.ratio = 1.5 + 0.1 * length(unique(to_plot$ID)) ,
		legend.position = "bottom",
		#text = element_text(size = 15),
		axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5),
		panel.grid = element_blank()
	) + 
	xlab("") + 
	ylab("Count") + 
	labs(fill = "Read type") + 
	coord_flip() + 
	scale_fill_manual(values = brewer.pal(8, "Set1")) + 
	guides(fill = guide_legend(ncol = 1, position = "bottom", direction = "vertical")) #+ 
	#facet(to_plot, conditions) 
	
plot(g)
dev.off()

plot_variables = c("Uniquely mapped reads %",  "% of reads unmapped: too short", "% of reads mapped to multiple loci")
to_plot = subset(melted, melted$variable %in% plot_variables, droplevels = T)
to_plot$variable = factor(as.character(to_plot$variable), levels = plot_variables)
to_plot$ID = gsub("Sample_", "", as.character(to_plot$ID))


pdf(file.path(output_dir, "read_stats_perc_line.pdf"), width = 6, height = 4 + 0.1 * length(unique(to_plot$ID)), family = "Helvetica", useDingbats = F)
g <- ggplot(to_plot, aes(x = ID, y = value, fill = variable)) +
	geom_bar(stat = "identity", position=position_dodge(), width = 0.75) + 
	theme_bw() + 
	ylim(0, 100) + 
	#scale_y_continuous() + #breaks=pretty_breaks(n=10)) +  
	theme(
		aspect.ratio = 1.5 + 0.1 * length(unique(to_plot$ID)) ,
		legend.position = "bottom",
		#text = element_text(size = 20),
		#axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
		panel.grid = element_blank()
		#axis.text.y = element_text(hjust = 0.5, vjust = 1)
	) + 
	xlab("") + 
	ylab("Percentage reads") + 
	labs(fill = "Read type") + 
	coord_flip() + 
	scale_fill_manual(values = brewer.pal(8, "Set1")) + 
	guides(fill = guide_legend(ncol = 1, position = "bottom", direction = "vertical", title = "")) #+ 
	#facet(to_plot, conditions)
	
plot(g)
dev.off()