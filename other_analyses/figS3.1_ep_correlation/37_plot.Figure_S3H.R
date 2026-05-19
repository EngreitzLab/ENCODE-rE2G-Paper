#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library("extrafont")
options(scipen = 999)

broad_broad = fread(file = "work/Classes.broad_broad.PRC.data", sep='\t', header=T, stringsAsFactors=FALSE)
broad_broad$group = "g:broad\ne:broad"

broad_specific = fread(file = "work/Classes.broad_specific.PRC.data", sep='\t', header=T, stringsAsFactors=FALSE)
broad_specific$group = "g:broad\ne:specific"

variable_broad = fread(file = "work/Classes.variable_broad.PRC.data", sep='\t', header=T, stringsAsFactors=FALSE)
variable_broad$group = "g:variable\ne:broad"

variable_specific = fread(file = "work/Classes.variable_specific.PRC.data", sep='\t', header=T, stringsAsFactors=FALSE)
variable_specific$group = "g:variable\ne:specific"

plot_data = rbindlist(list(broad_broad,broad_specific,variable_broad,variable_specific),use.names=TRUE)

plot_data$group = factor(plot_data$group, levels=c("g:broad\ne:broad","g:broad\ne:specific","g:variable\ne:broad","g:variable\ne:specific"))

plot = ggplot(plot_data, aes(x=recall, y=precision, fill=group)) + theme_bw()
plot = plot + geom_line(aes(color = group), linewidth = 0.3)
plot = plot + scale_color_manual(values=c("#808080","darkseagreen2","#ff7f00","#4E79A7"))

plot = plot + theme(    legend.position = "inside", legend.position.inside = c(0.67,0.73), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
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

plot = plot + labs(x = "Recall")
plot = plot + labs(y = "Precision")
plot = plot + ggtitle("CRISPRi in K562 cells")

plot = plot + scale_x_continuous(limits = c(0, 1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
plot = plot + scale_y_continuous(limits = c(0, 1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))

plot = plot + geom_hline(yintercept=0, color='black', linewidth=0.3, linetype='dotted')
plot = plot + geom_vline(xintercept=0, color='black', linewidth=0.3, linetype='dotted')

plot = plot + theme(plot.margin = margin(t = 8.9, r = 5, b = 4.9, l = 20, unit = 'pt'))

ggsave(file="Figure_S3H.pdf", plot = plot, width = 75, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3H.svg", plot = plot, width = 75, height = 65, units = "mm", dpi = 1200)
