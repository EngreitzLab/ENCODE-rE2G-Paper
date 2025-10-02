suppressPackageStartupMessages({
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)})

## data entry
# all file names
all.files = c("ENCSR154AVE.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR622RJN.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR374YBC.chr11.vs.ENCSR154AVE.chrX.tsv",
             "ENCSR577TXK.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR945ENZ.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR101TGA.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR137EEM.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR914PJE.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR967DXR.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR643OAB.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR826CTO.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR274JIR.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR141VJJ.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR982JMP.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR627GNV.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR087QUG.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR165BGV.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR901BRM.chr16.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR229PFU.chr9.vs.ENCSR622RJN.chr11.tsv", "ENCSR527IGO.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR167DFU.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR724SCF.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR632CQZ.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR815JAA.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR048BUO.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR560CCF.chr11.vs.ENCSR154AVE.chrX.tsv", "ENCSR997TXZ.chrX.vs.ENCSR622RJN.chr11.tsv",
             "ENCSR423UEH.chrX.vs.ENCSR622RJN.chr11.tsv", "ENCSR064EFK.chr11.vs.ENCSR154AVE.chrX.tsv",
             "ENCSR064EFK.chr11.vs.ENCSR154AVE.chrX.chr9.ROI.tsv", "ENCSR064EFK.chr11.vs.ENCSR154AVE.chrX.chr17.ROI.tsv")

all.files = paste0("/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles_noFilter/", all.files)

## add experimental data to cts

df.all = read.table(all.files[1], header=TRUE, sep="\t")
for (i in 2:length(all.files)){
  temp = read.table(all.files[i], header=TRUE, sep="\t")
  df.all = rbind(df.all, temp)
}

## output

file.name = paste0("/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles_noFilter/allFiles")
write.table(df.all, file=paste0(file.name, ".tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# box plot
# filter based on logFC SE and padj
#df = dplyr::filter(df, lfcSE<1)
#df = dplyr::filter(df, lfcSE<0.1 | padj<0.05)
df = dplyr::filter(df.all, lfcSE<0.15 | padj<0.01)
# df = df.all
#df = dplyr::filter(df.all, lfcSE<0.15)


# add binned distance
# filter to distance>1 kb
df = dplyr::filter(df, distFromDeletion>1000)
df$binned.dist = ""
df$binned.dist[df$distFromDeletion<=10000] = "1-10kb"
n1 = filter(df, binned.dist=="1-10kb") %>% nrow()
df$binned.dist[df$distFromDeletion<=10000] = paste0("1-10kb\nN=", n1)

df$binned.dist[df$distFromDeletion>10000 & df$distFromDeletion<100000] = "10-100kb"
n2 = filter(df, binned.dist=="10-100kb") %>% nrow()
df$binned.dist[df$distFromDeletion>10000 & df$distFromDeletion<100000] = paste0("10-100kb\nN=", n2)

df$binned.dist[df$distFromDeletion>100000 & df$distFromDeletion<2E6] = ">100kb"
n3 = filter(df, binned.dist==">100kb") %>% nrow()
df$binned.dist[df$distFromDeletion>100000 & df$distFromDeletion<2E6] = paste0(">100kb\nN=", n3)

# are there non-other chromosome values >1Mb? None (2E6 = value set for other chromosome)
test = filter(df, distFromDeletion>1E6, distFromDeletion!=2E6)
# message("N entries with distance>1Mb: ", nrow(test))

df$binned.dist[df$distFromDeletion==2E6] = "Other chromosome"
n4 = filter(df, binned.dist=="Other chromosome") %>% nrow()
df$binned.dist[df$distFromDeletion==2E6] = paste0("Ctrls (other\nchrom.\nN=", n4)

df$binned.dist = factor(df$binned.dist)

# boxplot with other chrom
#pdf(file=paste0(file.name, ".boxplot.otherChrom.v2.pdf"), width=5, height=7)
df$binned.dist = ordered(df$binned.dist, levels = c(paste0("1-10kb\nN=", n1), paste0("10-100kb\nN=", n2),
                                                    paste0(">100kb\nN=", n3), paste0("Ctrls (other\nchrom.\nN=", n4)))
fig = ggplot(df, aes(x = binned.dist,y = log2FoldChange)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1), size=0.5, alpha=0.5) +  geom_hline(yintercept=0, lty='solid') +
  geom_hline(yintercept=0, lty='dashed', size=0.5, color='red') +
  labs(x="Distance between\ndeletion and enhancer", y="DNase-seq signal fold-change\n(log2 Enhancer KO / Ctrl)") + 
  #ylim(c(-5, 5)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

#ggsave(file.path('/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles/', 'main_fig4f_ylim.eps'), fig, height=3, width=2.5)
ggsave(file.path('/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles_noFilter/', 'main_fig4f_noFilter.lfcSE0.15ORpadj0.01.pdf'), fig, height=3, width=2.5)

#dev.off()

hist = ggplot(df.all, aes(lfcSE)) + 
  geom_histogram(binwidth=0.05) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
ggsave(file.path('/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles_noFilter/', 'lfcSE_histogram.pdf'), hist, height=3, width=2.5)

scatter = ggplot(df.all, aes(lfcSE, baseMean)) +
  geom_point(alpha=0.5, size=0.5) + 
  xlim(c(0, 0.75)) + ylim(c(0, 3000)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
ggsave(file.path('/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/GM_DNase_Stam/outFiles_noFilter/', 'lfcSE_vs_baseMean.pdf'), scatter, height=3, width=2.5)

## Wilcoxon rank sum test
g1 = filter(df, binned.dist== paste0("1-10kb\nN=", n1))
g2 = filter(df, binned.dist==paste0("10-100kb\nN=", n2))
g3 = filter(df, binned.dist==paste0(">100kb\nN=", n3))
g4 = filter(df, binned.dist==paste0("Ctrls (other\nchrom.\nN=", n4))

g1.v.g3 = wilcox.test(g1$log2FoldChange, g3$log2FoldChange)
g1.v.g2 = wilcox.test(g1$log2FoldChange, g2$log2FoldChange)
g2.v.g3 = wilcox.test(g2$log2FoldChange, g3$log2FoldChange)

g1.v.g4 = wilcox.test(g1$log2FoldChange, g4$log2FoldChange)
g2.v.g4 = wilcox.test(g2$log2FoldChange, g4$log2FoldChange)
g3.v.g4 = wilcox.test(g3$log2FoldChange, g4$log2FoldChange)

print(g1.v.g3)
print(g1.v.g2)
print(g2.v.g3)

print(g1.v.g4)
print(g2.v.g4)
print(g3.v.g4)

## t-test
t1 = t.test(g1$log2FoldChange, mu=0)
t2 = t.test(g2$log2FoldChange, mu=0)
t3 = t.test(g3$log2FoldChange, mu=0)
t4 = t.test(g4$log2FoldChange, mu=0)

# print("Bin 1")
# print(t1)
# print("Bin 2")
# print(t2)
# print("Bin 3")
# print(t3)
# print("Bin 4")
# print(t4)


                                            
                                            
                                            