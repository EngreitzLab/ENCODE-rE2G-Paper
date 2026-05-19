#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(showtext)
font_add(family = "Helvetica", regular = "./fonts/Helvetica.ttf", bold = "./fonts/Helvetica-Bold.ttf")

EQTLs  = fread(file = "work/GTEx_ge_LCL.eQTLs.positive.data.gz", sep='\t', header=T, stringsAsFactors=FALSE)
Random = fread(file = "work/GTEx_ge_LCL.eQTLs.match.data.gz", sep='\t', header=T, stringsAsFactors=FALSE)

cat("EQTLs:  ",mean(EQTLs$glsCoefficient),"\n")
cat("Random: ",mean(Random$glsCoefficient),"\n")

analysis = wilcox.test(EQTLs$glsCoefficient,Random$glsCoefficient)
print(analysis)
print(analysis$p.value)

cat("\nDone wilcox.test\n\n")

##########################################
##########################################
##########################################
EQTLs$group  = "regulated"
Random$group = "notregulated"

plot_data = rbindlist(list(EQTLs,Random), use.names=TRUE, fill=TRUE)
###have to do this for ploting. coord_cartesian(ylim = c(-1, 2.5) does not limit the plot
plot_data = plot_data[glsCoefficient > 2.5, glsCoefficient:=2.5]
plot_data = plot_data[glsCoefficient < -1, glsCoefficient:=-1]

plot_data$group = factor(plot_data$group, levels = c("regulated","notregulated"))
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
plot = plot + scale_fill_manual(values=c("#008000","#ff66ff"))

plot = plot + theme(    legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.margin = margin(0, 0, 0, 0))
plot = plot + theme(    legend.title = element_blank(), legend.key.size = unit(0.9, 'lines'), legend.text=element_text(size=5, family="Helvetica", color='black', hjust=1),
                        legend.key = element_rect(color = "transparent", fill = "transparent"), legend.background=element_rect(color = "transparent", fill = "transparent"),)

plot = plot + theme(    plot.title   = element_text(vjust=0.5, hjust=0.0, size=7, family="Helvetica", color='black', margin = margin(t = 1, r = 0.5, b = 3, l = 0.5)),
                        axis.title.y = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 3.0, b = 1, l = 0.5)),
                        axis.title.x = element_blank(),
                        axis.text.y  = element_text(vjust=0.5, hjust=1.0, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 1, l = 0.5)),
                        axis.text.x  = element_text(vjust=0.5, hjust=0.5, size=6, family="Helvetica", color='black',  margin = margin(t = 1, r = 0.5, b = 0, l = 0.5)),
                        panel.border = element_blank(), panel.background = element_blank())

plot = plot + theme(axis.line  = element_line(colour = 'black', linewidth = 0.2))
plot = plot + theme(axis.ticks = element_line(colour = "black", linewidth = 0.2), axis.ticks.length=unit(2.5, "points"))

plot = plot + labs(y = 'GLS Coefficient (DNase-seq in 89 biosamples)')
plot = plot + ggtitle("eQTL in lymphoblastoid cells")

plot = plot + scale_y_continuous(breaks=c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5))
plot = plot + scale_x_discrete(labels=c('eQTLs\n(N=844)','Matched pairs\n(N=3,366)'))

plot = plot + theme(plot.margin = margin(t = 7, r = 13, b = 4.9, l = 12, unit = 'pt'))

shift = 0.11
plot = plot + annotate(geom="text",x=1, y=0.2898974 - shift, label="0.2899", size=2.5)
plot = plot + annotate(geom="text",x=2, y=0.3447724 - shift, label="0.3447", size=2.5)

plot = plot + annotate(geom="segment",x=1.1, xend=1.9, y=2.2, yend=2.2, color = "black", linewidth = 0.2)
plot = plot + annotate(geom="text",x=1.5, y=2.35, label = "p=0.2141", size=2.5)

ggsave(file="Figure_S3D.pdf", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)
ggsave(file="Figure_S3D.svg", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)

#showtext_auto()
#showtext_opts(dpi = 1200)
#ggsave(file="Figure_S3D.png", plot = plot, width = 52, height = 65, units = "mm", dpi = 1200)
warnings()
