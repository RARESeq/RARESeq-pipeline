library(jsonlite)
library(reshape2)
library(ggplot2)

args = commandArgs(T)

input_dir = args[1]
output_dir = dirname(args[1])

files = sort(list.files(input_dir, pattern = "*.json"))

all_data = data.frame()

for(file in files){
	s = strsplit(strsplit(basename(file), "_R1.", fixed = T)[[1]][1], "__")[[1]]
	sample_name = strsplit(basename(file), "_R1.", fixed = T)[[1]][1]
	con  <- file(file.path(input_dir, file), open = "r")
	rd <- readLines(con, warn = "F")
	df <- fromJSON(rd)
	close(con)
	data = data.frame(
		#"Total reads" = df$filtering_result$total_reads,
		"Good reads" = df$filtering_result$passed_filter_reads / (df$filtering_result$passed_filter_reads + df$filtering_result$low_quality_reads),
		#"Bad reads" = df$filtering_result$bad_reads / df$filtering_result$total_reads,
		"Reads with low quality" = df$filtering_result$low_quality_reads / (df$filtering_result$passed_filter_reads + df$filtering_result$low_quality_reads)
		)
	data$ID = s[1]
	data$Sample = sample_name
	
	all_data = rbind(all_data, data)
}
write.table(all_data, file.path(output_dir, "fastp_qc_data.txt"), row.names = F, sep = "\t")
melted = melt(all_data, id.vars = c("ID", "Sample"))
melted$variable = gsub(".", " ", melted$variable, fixed = T)
melted$Sample = gsub("Sample_", "", as.character(melted$Sample), fixed = T)
	

pdf(file.path(output_dir, "fastp_barplots.pdf"), width = 7, height = 6, family = "Helvetica", useDingbats = F)

g <- ggplot(melted, aes(x = Sample, y = value, fill = variable)) +
	geom_bar( stat = "identity", position=position_dodge()) +
	#geom_bar( stat = "identity", width = 1) + #position=position_dodge()) +
	scale_fill_brewer(palette = "Set1") + 
	theme_bw() + 
	ylim(0, 1) + 
	ylab("Percentage") + 
	theme_bw() + 
	theme(
		aspect.ratio = 1.5,
		#text = element_text(size = 25),
		#axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		panel.grid = element_blank(),
		#panel.grid.major.x = element_line(color = "#44444444"),
		#legend.position = "none",
		#panel.grid=element_blank(),
		#axis.ticks = element_blank(),
		#axis.text.x = element_blank()) + 
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + 
	#coord_polar("y") +
	coord_flip() + 
	ylab("Fraction") + 
	guides(fill = guide_legend(nrow = 2, title = "Read type"))
	#facet_wrap(~Sample)
	#facet(melted, conditions)
plot(g)
dev.off()
#save(g, file = file.path(output_dir, "plot.RData")) 