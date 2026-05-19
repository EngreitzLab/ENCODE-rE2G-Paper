#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")

broad_broad = fread(file = "work/Classes.broad_broad.data", sep='\t', header=T, stringsAsFactors=FALSE)
broad_broad$group = "g:broad\ne:broad"

broad_specific = fread(file = "work/Classes.broad_specific.data", sep='\t', header=T, stringsAsFactors=FALSE)
broad_specific$group = "g:broad\ne:specific"

variable_broad = fread(file = "work/Classes.variable_broad.data", sep='\t', header=T, stringsAsFactors=FALSE)
variable_broad$group = "g:variable\ne:broad"

variable_specific = fread(file = "work/Classes.variable_specific.data", sep='\t', header=T, stringsAsFactors=FALSE)
variable_specific$group = "g:variable\ne:specific"

plot_data = rbindlist(list(broad_broad,broad_specific,variable_broad,variable_specific),use.names=TRUE)
plot_data = plot_data[Regulated == "TRUE" ]

plot_data$group = factor(plot_data$group, levels=c("g:broad\ne:broad","g:broad\ne:specific","g:variable\ne:broad","g:variable\ne:specific"))

qqq = plot_data[, mean(glsCoefficient), by="group"]
colnames(qqq) = c("group","mean")
print(qqq)
cat("\n")
qqq = plot_data[, max(glsCoefficient), by="group"]
colnames(qqq) = c("group","max")
print(qqq)
cat("\n")
qqq = plot_data[, min(glsCoefficient), by="group"]
colnames(qqq) = c("group","min")
print(qqq)
cat("\n")

plot = ggplot(plot_data, aes(y=glsCoefficient, x=group, fill=group)) + theme_bw()
plot = plot + geom_violin(scale = "width", linewidth=0.2, bw=0.15, width=0.95) + coord_cartesian(ylim = c(-1, 2.5))
plot = plot + stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black", linewidth=0.1, position = position_dodge(0.9))
plot = plot + scale_fill_manual(values=c("#808080","darkseagreen2","#ff7f00","#4E79A7"))


plot = plot + theme(    legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(0.9, 'lines'), legend.text=element_text(size=5, family="Helvetica", color='black', hjust=1),
                        legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=0.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
#                        axis.title.x = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
			axis.title.x = element_blank(),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

plot = plot + labs(y = 'GLS Coefficient (DNase-seq in 89 biosamples)')
plot = plot + ggtitle("CRISPRi in K562 cells")

plot = plot + scale_y_continuous(breaks=c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5))

plot = plot + theme(plot.margin = margin(t = 7, r = 7, b = 4.9, l = 12, unit = 'pt'))

ggsave(file="Figure_S3I.pdf", plot = plot, width = 82, height = 64.1, units = "mm", dpi = 1200)
ggsave(file="Figure_S3I.svg", plot = plot, width = 82, height = 64.1, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3C.png", plot = plot, width = 82, height = 65, units = "mm", dpi = 1200)

warnings()
