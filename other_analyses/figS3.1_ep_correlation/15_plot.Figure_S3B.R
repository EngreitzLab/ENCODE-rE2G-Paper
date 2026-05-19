#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")
options(scipen = 999)

Pos_glsCoefficient = fread(file = "work/EPCrisprBenchmark_Pos_glsCoefficient.density", sep='\t', header=T, stringsAsFactors=FALSE)
Pos_glsCoefficient$group = "CRISPRi regulated\npairs (N=438)"

Neg_glsCoefficient = fread(file = "work/EPCrisprBenchmark_Neg_glsCoefficient.density", sep='\t', header=T, stringsAsFactors=FALSE)
Neg_glsCoefficient$group = "CRISPRi not regulated\npairs (N=8,550)"

Neg_matched = fread(file = "work/EPCrisprBenchmark_Match_glsCoefficient.density", sep='\t', header=T, stringsAsFactors=FALSE)
Neg_matched$group = "CRISPRi matched\nnegative pairs (N=300)"

GTEx_regulatory = fread(file = "work/GTEx_ge_LCL.eQTLs.positive.density", sep='\t', header=T, stringsAsFactors=FALSE)
GTEx_regulatory$group = "Lymphoblastoid eQTLs\nfrom GTEx (N=844)"

GTEx_random = fread(file = "work/GTEx_ge_LCL.eQTLs.match.density", sep='\t', header=T, stringsAsFactors=FALSE)
GTEx_random$group = "Non-significant eQTLs\n(N=3,366)"

plot_data = rbindlist(list(Pos_glsCoefficient,Neg_glsCoefficient,Neg_matched,GTEx_regulatory,GTEx_random),use.names=TRUE)

plot_data$group = factor(plot_data$group, levels=c("CRISPRi regulated\npairs (N=438)","CRISPRi not regulated\npairs (N=8,550)","CRISPRi matched\nnegative pairs (N=300)","Lymphoblastoid eQTLs\nfrom GTEx (N=844)","Non-significant eQTLs\n(N=3,366)"))

plot_data = plot_data[ distance <= 205000]

plot = ggplot(plot_data, aes(y=density, x=distance, fill=group)) + theme_bw()
plot = plot + geom_line(aes(color = group), linewidth = 0.4)
plot = plot + scale_color_manual(values=c("#4c79ff","#808080","#ff7f00","#008000","#ff66ff"))


plot = plot + theme(	legend.position = "inside", legend.position.inside = c(0.70,0.73), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(0.9, 'lines'), legend.text=element_text(size=5, family="Helvetica", color='black', hjust=1),
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
#plot = plot + ggtitle("GTEx eQTLs in lymphoblastoid cells")
plot = plot + scale_y_continuous(limits = c(0, 0.05),  breaks=c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05))
plot = plot + scale_x_continuous(limits = c(0, 205000), breaks=c(0, 100000, 200000), labels = c("0", "100,000", "200,000"))

plot = plot + theme(plot.margin = margin(t = 21.00, r = 7, b = 3, l = 12, unit = 'pt'))

ggsave(file="Figure_S3B.pdf", plot = plot, width = 70, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3B.svg", plot = plot, width = 70, height = 65, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3B.png", plot = plot, width = 70, height = 65, units = "mm", dpi = 1200)
