#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")

DNase_data = fread(file = "work/DNase_seq.example.matrix", sep='\t', header=T, stringsAsFactors=FALSE)

gene_id     = "CD69"
enhancer_id = "chr12:9764781-9765131"

promoter_activity = DNase_data[ id == gene_id ]
enhancer_activity = DNase_data[ id == enhancer_id]

dot_size = 1.2

promoter_activity[,id:=NULL]
promoter_activity = as.numeric(promoter_activity)

enhancer_activity[,id:=NULL]
enhancer_activity = as.numeric(enhancer_activity)

plot_data = data.frame(promoter = promoter_activity, enhancer = enhancer_activity)
cat("Max promoter:",max(promoter_activity),"\n")
cat("Max enhancer:",max(enhancer_activity),"\n")

plot = ggplot(data = plot_data, aes(x=promoter, y=enhancer)) + theme_bw()
plot = plot + stat_smooth(method = "lm", se = T, color="black", linewidth=0, linetype = "dotted", fullrange=T, level=0.95)
plot = plot + geom_point(size=dot_size, stroke=0)

plot = plot + theme(    legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(1.1, 'lines'), legend.text=element_text(size=6, family="Helvetica", color='black', hjust=0),
                        legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=0.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
                        axis.title.x = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

plot = plot + labs(x = 'Promoter DNase-seq')
plot = plot + labs(y = 'Enhancer DNase-seq')
plot = plot + ggtitle(sprintf("%s enhancer at %s", gene_id, enhancer_id))

plot = plot + scale_y_continuous(breaks=c(0, 1, 2, 3), labels=c("0.0", "1.0", "2.0", "3.0"))
plot = plot + scale_x_continuous(limits = c(0, 5.4), breaks=c(0, 1, 2, 3, 4, 5), labels=c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0"))

plot = plot + theme(plot.margin = margin(t = 7, r = 7, b = 3, l = 14.5, unit = 'pt'))

ggsave(file="Figure_S3K.pdf", plot = plot, width = 100, height = 40, units = "mm", dpi = 1200)
ggsave(file="Figure_S3K.svg", plot = plot, width = 100, height = 40, units = "mm", dpi = 1200)
