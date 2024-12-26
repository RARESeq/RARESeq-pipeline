library(ggplot2)
library(RColorBrewer)
library(reshape2)

args = commandArgs(T) 

input_dir = args[1]
dedupe_dir = args[2]
output_dir = input_dir
all_stats = NULL

pre = grepl("pre-ontarget-filtering", input_dir)

colors = rep(
	c(brewer.pal(9, "Set1")[-c(3)], 
	  brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], 
	  brewer.pal(8, "Dark2")[-c(5,7,8)],
	  brewer.pal(9, "Set1")[c(3)], 
	  brewer.pal(12, "Paired")[c(2, 4, 6, 8, 10, 11, 12)], 
	  brewer.pal(8, "Set2"),
	  brewer.pal(8, "Set3"),
	  brewer.pal(8, "Pastel1"),
	  brewer.pal(8, "Pastel2"),
	  brewer.pal(8, "Accent")), 50) 
colors = gsub("#FFFF33", "gold", colors)

folders = list.files(input_dir)
for(folder in folders)
{
	inp = file.path(input_dir, folder)
	metrics_file = list.files(inp, pattern = "*.dualindex-deduped.sorted.nr_reads$")
	if(length(metrics_file) == 0)
	{
		next
	}

	deduped_reads = read.delim(file.path(inp, metrics_file[1]), header = F)[1,1]
	non_deduped_reads = read.delim(file.path(inp, gsub(".sorted.dualindex-deduped", "", metrics_file[1])), header = F)[1,1]
	data = read.delim(file.path(inp, gsub(".dualindex-deduped.sorted.nr_reads", ".dualindex-deduped.sorted.exons_cov", metrics_file[1])), header = F)
	fraction = sum(data[data[,1] != 0,2]) / sum(data[,2])
	sample = strsplit(metrics_file[1], "__")[[1]][1]
	sample = strsplit(sample, ".dualindex-deduped.")[[1]][1]
	sample = gsub(".sorted", "", sample)		
		
	df = data.frame(Sample = sample, Non_deduped_molecules = as.integer(non_deduped_reads), Deduped_reads_molecules = as.integer(deduped_reads), Fraction_exonic_bases_covered = fraction)
	all_stats = rbind(all_stats, df)
}

samples = sort(unique(all_stats$Sample))
for(sample in samples)
{
	inp = file.path(dedupe_dir, sample)
	
	metrics_file = paste0(sample, ".sorted.dualindex-deduped.sorted.nr_reads")
	deduped_reads = read.delim(file.path(inp, metrics_file), header = F)[1,1]
	
	if(pre)
	{
		non_deduped_reads = read.delim(file.path(inp, gsub(".sorted.dualindex-deduped", "", metrics_file)), header = F)[1,1]				
	}else{		
		non_deduped_reads = read.delim(file.path(inp, gsub(".sorted.dualindex-deduped", ".sorted.filt", metrics_file)), header = F)[1,1]	
	}
	
	data = read.delim(file.path(inp, gsub(".dualindex-deduped.sorted.nr_reads", ".dualindex-deduped.sorted.exons_cov", metrics_file)), header = F)
	fraction = sum(data[data[,1] != 0,2]) / sum(data[,2])
		
	df = data.frame(Sample = sample, Non_deduped_molecules = as.integer(non_deduped_reads), Deduped_reads_molecules = as.integer(deduped_reads), Fraction_exonic_bases_covered = fraction)
	all_stats = rbind(all_stats, df)
}

write.table(all_stats, file.path(output_dir, "combined_downsampling_data.tsv"), sep = "\t", row.names = F)

melted = melt(all_stats, id.vars = c("Sample", "Non_deduped_molecules"))
melted$Type = ifelse(melted$variable == "Deduped_reads_molecules", "Number of unique reads", "Fraction of exonic bases covered")
melted$Type = factor(as.character(melted$Type), levels = rev(levels(as.factor(melted$Type))))
melted$Sample = gsub("Sample_", "", as.character(melted$Sample))
melted$Sample = gsub("_cfrna", "", as.character(melted$Sample))
melted$Sample = gsub("_tumor", "", as.character(melted$Sample))
melted$Sample = gsub("_normal", "", as.character(melted$Sample))

pdf(file.path(output_dir, "downsampling_plot.pdf"), height = 5.5 + 0.1 * length(unique(melted$Sample)), width = 8, family = "Helvetica", useDingbats = F) 
g <- ggplot(melted, aes(x = Non_deduped_molecules, y = value)) + 
	geom_point(aes(color = Sample), size = 2) + 	
	geom_line(aes(group = Sample, color = Sample), lwd = 1)+ 
	theme_bw() + 
	theme(panel.grid = element_blank()) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
	theme(aspect.ratio = 1) + 
	theme(strip.background = element_blank(), strip.placement = "outside") + 
	facet_wrap(~Type, scales = "free_y", nrow = 1, strip.position = "left") +
                #labeller = as_labeller(c(A = "Currents (A)", V = "Voltage (V)") ) )  +
    ylab(NULL) + 
	xlab("Number of input reads") + 
	scale_color_manual(values= colors) +
	guides(color = guide_legend(ncol = 3, title = "")) + 
	theme(legend.position = "bottom") + 
	theme(strip.text.x = element_text(size = 12)) +
	theme(strip.text.y  = element_text(size = 12)) +
	theme(axis.title.x = element_text(size = 12)) 

plot(g)
dev.off()

pdf(file.path(output_dir, "downsampling_plot_interpolated.pdf"), height = 5.5 + 0.1 * length(unique(melted$Sample)), width = 8, family = "Helvetica", useDingbats = F) 
splits = split(melted, paste(melted$Type, melted$Sample))
 df = do.call(rbind, lapply(splits, function(spl){
	spl = spl[order(spl$Non_deduped_molecules),]
	MonotonicSpline <- splinefun(x=spl$Non_deduped_molecules, y = spl$value, method = "monoH.FC")
	x = seq(min(spl$Non_deduped_molecules), max(spl$Non_deduped_molecules), length.out = 1000)
	monotonicFit <- MonotonicSpline(x, extrapol = "linear")	 
	data.frame(Non_deduped_molecules = x, value = monotonicFit, Sample = spl$Sample[1], Type = spl$Type[1])
	}))

g <- ggplot(melted, aes(x = Non_deduped_molecules, y = value)) + 
	geom_point(aes(color = Sample), size = 2) + 
	#geom_smooth(aes(group = Sample, color = Sample), se=FALSE,method = "glm",
    #          family = gaussian()) +  
    geom_line(aes(group = Sample, color = Sample), data=df, lwd = 1)+
	#geom_spline(aes(group = Sample, color = Sample))+ 
	#geom_line(aes(group = Sample, color = Sample))+ 
	theme_bw() + 
	theme(panel.grid = element_blank()) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
	theme(aspect.ratio = 1) + 
	theme(strip.background = element_blank(), strip.placement = "outside") + 
	facet_wrap(~Type, scales = "free_y", nrow = 1, strip.position = "left") +
                #labeller = as_labeller(c(A = "Currents (A)", V = "Voltage (V)") ) )  +
    ylab(NULL) +
	xlab("Number of input reads") + 
	scale_color_manual(values= colors) +
	guides(color = guide_legend(ncol = 3, title = "")) + 
	theme(legend.position = "bottom") + 
	theme(strip.text.x = element_text(size = 12)) +
	theme(strip.text.y  = element_text(size = 12)) +
	theme(axis.title.x = element_text(size = 12)) 

plot(g)

dev.off()
