#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")
options(scipen = 999)

broad_broad = fread(file = "work/Classes.broad_broad.density", sep='\t', header=T, stringsAsFactors=FALSE)
broad_broad$group = "g:broad\ne:broad"

broad_specific = fread(file = "work/Classes.broad_specific.density", sep='\t', header=T, stringsAsFactors=FALSE)
broad_specific$group = "g:broad\ne:specific"

variable_broad = fread(file = "work/Classes.variable_broad.density", sep='\t', header=T, stringsAsFactors=FALSE)
variable_broad$group = "g:variable\ne:broad"

variable_specific = fread(file = "work/Classes.variable_specific.density", sep='\t', header=T, stringsAsFactors=FALSE)
variable_specific$group = "g:variable\ne:specific"

plot_data = rbindlist(list(broad_broad,broad_specific,variable_broad,variable_specific),use.names=TRUE)

plot_data$group = factor(plot_data$group, levels=c("g:broad\ne:broad","g:broad\ne:specific","g:variable\ne:broad","g:variable\ne:specific"))

plot_data = plot_data[ distance <= 1000000]

plot = ggplot(plot_data, aes(y=density, x=distance, fill=group)) + theme_bw()
plot = plot + geom_line(aes(color = group), linewidth = 0.4)
plot = plot + scale_color_manual(values=c("#808080","darkseagreen2","#ff7f00","#4E79A7"))

plot = plot + theme(	legend.position = "inside", legend.position.inside = c(0.47,0.73), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(1.3, 'lines'), legend.text=element_text(size=6, family="Helvetica", color='black', hjust=1),
                        legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=1.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
                        axis.title.x = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())
plot = plot + guides(color = guide_legend(override.aes = list(size = 1.5, stroke = 0)))

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

plot = plot + labs(y = 'Density')
plot = plot + labs(x = 'TSS-element distance')
plot = plot + scale_y_continuous(limits = c(0, 0.05),  breaks=c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05))
plot = plot + scale_x_continuous(limits = c(0, 1000000), breaks=c(0, 500000, 1000000), labels = c("0","500,000","1,000,000"))

plot = plot + theme(plot.margin = margin(t = 21.00, r = 11.8, b = 3, l = 12, unit = 'pt'))

ggsave(file="Figure_S3G.pdf", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3G.svg", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3B.png", plot = plot, width = 70, height = 65, units = "mm", dpi = 1200)
