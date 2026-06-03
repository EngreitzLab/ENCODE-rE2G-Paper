#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")

Pos_glsCoefficient   = fread(file = "work/EPCrisprBenchmark_Pos_glsCoefficient.data", sep='\t', header=T, stringsAsFactors=FALSE)
Neg_glsCoefficient   = fread(file = "work/EPCrisprBenchmark_Neg_glsCoefficient.data", sep='\t', header=T, stringsAsFactors=FALSE)
Match_glsCoefficient = fread(file = "work/EPCrisprBenchmark_Match_glsCoefficient.data", sep='\t', header=T, stringsAsFactors=FALSE)

cat("Pos versus Neg:\n")
analysis = wilcox.test(Pos_glsCoefficient$glsCoefficient,Neg_glsCoefficient$glsCoefficient)
print(analysis)
print(analysis$p.value)

cat("Pos versus Match:\n")
analysis = wilcox.test(Pos_glsCoefficient$glsCoefficient,Match_glsCoefficient$glsCoefficient)
print(analysis)
print(analysis$p.value)
cat("\nDone wilcox.test\n\n")

###################################
###################################
###################################
Pos_glsCoefficient$group = "regulated"
Neg_glsCoefficient$group = "notregulated"
Match_glsCoefficient$group = "matched"

plot_data = rbindlist(list(Neg_glsCoefficient,Pos_glsCoefficient,Match_glsCoefficient),use.names=TRUE)

plot_data$group = factor(plot_data$group, levels = c("notregulated","regulated","matched"))
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
plot = plot + geom_violin(scale = "width", linewidth=0.2, bw=0.10, width=0.95) + coord_cartesian(ylim = c(-1, 2.5))
plot = plot + stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black", linewidth=0.1, position = position_dodge(0.9))
plot = plot + scale_fill_manual(values=c("#808080", "#4c79ff","#ff7f00ff"))

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
plot = plot + scale_x_discrete(labels=c('Not regulated\npairs (N=8,550)','Regulated\npairs (N=438)','Distance matched\nneg pairs (N=300)'))

plot = plot + theme(plot.margin = margin(t = 7, r = 7, b = 4.9, l = 12, unit = 'pt'))

shift = 0.11
plot = plot + annotate(geom="text",x=1, y=0.09391446 - shift, label="0.0939", size=2.5)
plot = plot + annotate(geom="text",x=2, y=0.32087345 - shift, label="0.3209", size=2.5)
plot = plot + annotate(geom="text",x=3, y=0.13114130 - shift, label="0.1311", size=2.5)

plot = plot + annotate(geom="segment",x=1.1, xend=1.9, y=2.1, yend=2.1, color = "black", linewidth = 0.2)
plot = plot + annotate(geom="text",x=1.5, y=2.25, label = "p=3.13e-55", size=2.8)

plot = plot + annotate(geom="segment",x=2.1, xend=2.9, y=2.1, yend=2.1, color = "black", linewidth = 0.2)
plot = plot + annotate(geom="text",x=2.5, y=2.25, label = "p=3.28e-14", size=2.8)

ggsave(file="Figure_S3C.pdf", plot = plot, width = 82, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3C.svg", plot = plot, width = 82, height = 65, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3C.png", plot = plot, width = 82, height = 65, units = "mm", dpi = 1200)

warnings()
