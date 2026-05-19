#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")
options(scipen = 999)

input_data = fread(file = "work/Figure_S3A.smooth", sep='\t', header=T, stringsAsFactors=FALSE)

Mean_data   = data.table(distance = input_data$distance, signal = input_data$mean,   group = "mean")
Lower_data  = data.table(distance = input_data$distance, signal = input_data$lower,  group = "lower")
Higher_data = data.table(distance = input_data$distance, signal = input_data$higher, group = "higher")
plot_data   = rbindlist(list(Mean_data,Lower_data,Higher_data), use.names=TRUE, fill=TRUE)

plot_data$group = factor(plot_data$group, levels = c("mean","lower","higher"))
plot_data = plot_data[ distance < 1000000 ]

plot = ggplot(plot_data, aes(y=signal, x=distance)) + theme_bw()

plot = plot + geom_line(aes(color = group, linewidth = group))
plot = plot + scale_color_manual(values=c("#000000", "#c2c2c2", "#c2c2c2"))
plot = plot + scale_linewidth_manual(values = c(0.5, 0.2, 0.2))

plot = plot + theme(	legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(	legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.text=element_text(size=6, family="Helvetica", color='black', hjust=1),
			legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=0.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
                        axis.title.x = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())
plot = plot + guides(color = guide_legend(override.aes = list(size = 1.5, stroke = 0)))

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

plot = plot + labs(y = 'GLS Coefficient (DNase-seq in 89 biosamples)')
plot = plot + labs(x = 'TSS-element distance')
plot = plot + ggtitle("GLS Coefficient genome-wide")

plot = plot + scale_y_continuous(limits = c(-1,2.5), breaks=c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5))
plot = plot + scale_x_continuous(limits = c(0,1000000), breaks=c(0, 500000, 1000000), labels = c("0","500,000","1,000,000"))

plot = plot + theme(plot.margin = margin(t = 7, r = 13, b = 3, l = 12, unit = 'pt'))

ggsave(file="Figure_S3A.pdf", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3A.svg", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3A.png", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)

warnings()
