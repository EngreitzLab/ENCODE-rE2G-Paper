#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library("extrafont")
options(scipen = 999)

Encode_E2G = fread(file = "work/PRC.Encode_E2G.data", sep='\t', header=T, stringsAsFactors=FALSE)
Encode_E2G$group = "Encode-E2G complete"

Encode_E2G_minus = fread(file = "work/PRC.Encode_E2G_minus.data", sep='\t', header=T, stringsAsFactors=FALSE)
Encode_E2G_minus$group = "Encode-E2G without\nGLS Coefficient"

ABC = fread(file = "work/PRC.ABC.data", sep='\t', header=T, stringsAsFactors=FALSE)
ABC$group = "ABC model"

ABC_plus = fread(file = "work/PRC.ABC_plus.data", sep='\t', header=T, stringsAsFactors=FALSE)
ABC_plus$group = "ABC with GLS Coefficient"

Baseline = fread(file = "work/PRC.Baseline.data", sep='\t', header=T, stringsAsFactors=FALSE)
Baseline$group = "Baseline model"

Baseline_plus = fread(file = "work/PRC.Baseline_plus.data", sep='\t', header=T, stringsAsFactors=FALSE)
Baseline_plus$group = "Baseline with\nGLS Coefficient"

glsCoefficient = fread(file = "work/PRC.glsCoefficient.data", sep='\t', header=T, stringsAsFactors=FALSE)
glsCoefficient$group = "GLS Coefficient"

RNA_glsCoefficient = fread(file = "work/PRC.RNAglsCoefficient.data", sep='\t', header=T, stringsAsFactors=FALSE)
RNA_glsCoefficient$group = "RNA GLS Coefficient"

Pearson = fread(file = "work/PRC.pearsonCorrelation.data", sep='\t', header=T, stringsAsFactors=FALSE)
Pearson$group = "Pearson correlation"

plot_data = rbindlist(list(Encode_E2G,Encode_E2G_minus,ABC,ABC_plus,Baseline,Baseline_plus,glsCoefficient,RNA_glsCoefficient,Pearson),use.names=TRUE)
plot_data$group = factor(plot_data$group, levels = c("Encode-E2G complete","Encode-E2G without\nGLS Coefficient","ABC model","ABC with GLS Coefficient","Baseline model","Baseline with\nGLS Coefficient","GLS Coefficient","RNA GLS Coefficient","Pearson correlation"))

plot = ggplot(plot_data, aes(x=recall, y=precision, fill=group)) + theme_bw()
plot = plot + geom_line(aes(color = group), linewidth = 0.3)

plot = plot + theme(    legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(0.9, 'lines'), legend.text=element_text(size=6, family="Helvetica", color='black', hjust=0),
                        legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=0.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
                        axis.title.x = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

legend_labels = c("ENCODE-rE2G","ENCODE-rE2G\nwithout GLS Coefficient","ABC","ABC\nwith GLS Coefficient","Baseline","Baseline with\nGLS Coefficient","GLS Coefficient","RNA GLS Coefficient","Pearson correlation")
#legend_labels = c('Encode-rE2G'^Extended),"Encode-E2G\nwithout GLS Coefficient",bquote('ABC'^'A=DNase x H3K27ac, C=Intact Hi-C'),"ABC\nwith GLS Coefficient","Baseline","Baseline with\nGLS Coefficient","GLS Coefficient","RNA GLS Coefficient","Pearson correlation")

plot = plot + scale_color_manual( labels = legend_labels , values=c("#B07AA1","#e29dcf","#4E79A7","#114781","darkseagreen2","#79a079","#ff7f00","#808080","#000000"))

plot = plot + labs(x = "Recall")
plot = plot + labs(y = "Precision")
plot = plot + ggtitle("CRISPRi in K562 cells")


plot = plot + scale_x_continuous(limits = c(0, 1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
plot = plot + scale_y_continuous(limits = c(0, 1.0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))

plot = plot + geom_hline(yintercept=0, color='black', linewidth=0.3, linetype='dotted')
plot = plot + geom_vline(xintercept=0, color='black', linewidth=0.3, linetype='dotted')

plot = plot + theme(plot.margin = margin(t = 7, r = 132.0, b = 4.9, l = 20, unit = 'pt'))

ggsave(file="Figure_S3E.pdf", plot = plot, width = 152, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3E.svg", plot = plot, width = 152, height = 65, units = "mm", dpi = 1200)
