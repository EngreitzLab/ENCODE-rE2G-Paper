## Exploratory analyses for ENCODE Distal Regulation manuscript
## December 12, 2022
## Jesse Engreitz
## Scratch / notes

'sdev -p engreitz -m 15G -t 24:00:00
cd /home/groups/engreitz/Software/anaconda3
conda activate ./envs/EngreitzLab
cd /oak/stanford/groups/engreitz/Users/ejagoda/abc/global_properties
R'


library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(reshape2)

###work on the fig3 
#have to do on the server wayyy too big now
setwd("/oak/stanford/groups/engreitz/Users/ejagoda/encode_e2g_revision/")

k562 = read.table("encode_e2g_predictions_K562_ENCDO000AAD_ENCFF325RTP_DNaseOnly_full_predictions_just_chr21_1mil.txt",header=T)
k562_use = k562[k562$isSelfPromoter == FALSE & k562$class != "promoter",]
#colnames(k562_use) =c("enh_chr","enh_start","enh_end","enh_name","enh_type","gene_name")
#need to make it as a matrix
#gene_info = read.table("/Users/ejagoda/Documents/HGRM/EP_benchmarking/ABC/upded_7.7.22.RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.txt")

##gonna split in half because it's too big
x = acast(k562_use, k562_use$start ~ k562_use$TargetGeneTSS, value.var = "Score",fun.aggregate = sum)
m = as.matrix(x)

pdf("K562_chr21_top1mil_for_tensor_test_2.pdf")
heatmap(m,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562", heat.colors(n = 12))
dev.off()

png("K562_chr21_top1mil_for_tensor_test_2.png")
heatmap(m,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562", heat.colors(n = 12))
dev.off()

#try smaller 1st 5mb

testes = read.table("encode_e2g_predictions_testis_ENCDO836UAR_ENCFF882ETO_DNaseOnly_full_predictions_just_chr21_1mil.txt",header=T)
testes_use =testes[testes$isSelfPromoter == FALSE & testes$class != "promoter",]

testes_use_small = testes_use[as.numeric(paste0(testes_use$TargetGeneTSS)) <= 30000000 & as.numeric(paste0(testes_use$start)) <= 30000000,]
x = acast(testes_use_small, testes_use_small$start ~ testes_use_small$TargetGeneTSS, value.var = "Score",fun.aggregate = sum)
m = as.matrix(x)

pdf("testes_chr21_30mb.pdf")
heatmap(m,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main ="Testes", heat.colors(n = 12))
dev.off()




#for the small one
k562_use_small = k562_use[as.numeric(paste0(k562_use$TargetGeneTSS)) <= 5000000 & as.numeric(paste0(k562_use$start)) <= 5000000,]
x = acast(k562_use_small, k562_use_small$start ~ k562_use_small$TargetGeneTSS, value.var = "Score",fun.aggregate = sum)
m = as.matrix(x)

#115981 4778603
# ACAP3  WRAP73
pdf("K562_chr21_5mb.pdf")
heatmap(m,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562", heat.colors(n = 12))
dev.off()




library(gplots)

pdf("K562_chr21_for_tensor_w_line.pdf")
heatmap.2(m,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562",add.expr = abline(v=(37073183 + 5000000)))
dev.off()


m_use1 = m[,colnames(m) >= (37073183 - 2500000) & colnames(m) <= (37073183 + 2500000)]
m_use2 = m_use1[row.names(m_use1) >= (37073183 - 2500000) & row.names(m_use1) <= (37073183 + 2500000),]


m_use3 = m[,colnames(m) >= (37073183 - 5000000) & colnames(m) <= (37073183 + 5000000)]
m_use4 = m_use3[row.names(m_use3) >= (37073183 - 5000000) & row.names(m_use3) <= (37073183 + 5000000),]

k562_use_small = k562_use[as.numeric(paste0(k562_use$TargetGeneTSS)) >= (37073183 - 2500000) & as.numeric(paste0(k562_use$TargetGeneTSS))  <= (37073183 + 2500000),]





#x = acast(k562_use_small, k562_use_small$start ~ k562_use_small$TargetGeneTSS, value.var = "Full.Score",fun.aggregate = sum)
#m = as.matrix(x)
pdf("K562_chr21_for_tensor_5mb.pdf")
heatmap(m_use2,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562")
#legend(fill = heat.colors(12), legend = 1:12,x = "left" )
dev.off()

pdf("K562_chr21_for_tensor_10mb.pdf")
heatmap(m_use4,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "K562")
#legend(fill = heat.colors(12), legend = 1:12,x = "left" )
dev.off()


#testis small
testis = read.table("testis_chr21_only_not_filtered_DNaseOnly_features.tsv",header=T,sep = '\t')
testis_use = testis[testis$isSelfPromoter == "False",]

k562 = read.table("testis_chr21_only_not_filtered_DNaseOnly_features.tsv",header=T,sep = '\t')
k562_use = k562[k562$isSelfPromoter == "False",]
#colnames(k562_use) =c("enh_chr","enh_start","enh_end","enh_name","enh_type","gene_name")
#need to make it as a matrix
#gene_info = read.table("/Users/ejagoda/Documents/HGRM/EP_benchmarking/ABC/upded_7.7.22.RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.txt")

x = acast(k562_use, k562_use$start ~ k562_use$TargetGeneTSS, value.var = "Full.Score",fun.aggregate = sum)
tm = as.matrix(x)

tm_use3 = tm[,colnames(tm) >= (37073183 - 5000000) & colnames(tm) <= (37073183 + 5000000)]
tm_use4 = tm_use3[row.names(tm_use3) >= (37073183 - 5000000) & row.names(tm_use3) <= (37073183 + 5000000),]


pdf("testis_chr21_for_tensor.pdf")
heatmap(tm,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "Testis")
dev.off()

testis_use_small = testis_use[as.numeric(paste0(testis_use$TargetGeneTSS)) >= (37073183 - 5000000) & as.numeric(paste0(testis_use$TargetGeneTSS))  > (37073183 + 5000000),]

#x = acast(testis_use_small, testis_use_small$start ~ testis_use_small$TargetGeneTSS, value.var = "Full.Score",fun.aggregate = sum)
#m = as.matrix(x)
pdf("testis_chr21_for_tensor_10mb.pdf")
heatmap(tm_use4,dendrogram='none', Rowv=NA, Colv=NA,trace='none',xlab = "Genes",ylab = "Enhancers",main = "Testis")
dev.off()


##get total bp for the collapsed tab

library(stringr)
sums = c()
for (chr in 1:23){
  collapsed_tab = read.table(paste0("chr",chr,"_test_biosamples_enhancer_enhancers_collapse.txt"),header = T, sep = '\t')
  collapsed_tab$bp = 0
  for (i in 1:nrow(collapsed_tab)){
    enhancers = str_replace(pattern = paste0("chr",chr,":"), str_split(collapsed_tab[i,"enhancers_in_pair"],";")[[1]], "")
    pos = unlist(str_split(enhancers,"-"))
    bp = max(na.omit(as.numeric(paste0(pos))))-min(na.omit(as.numeric(paste0(pos))))
    collapsed_tab$bp[i] = bp
  }
  sums = c(sums,sum(collapsed_tab$bp))
}
sum(sums)
#22949408
#22949408/3000000000


#make summary table

loadEnhancersPerGene <- function(predictors, biosamples, predictor) {
  #dir <- predictors %>% filter(`pred_id`==predictor) %>% pull(`directory_rel_path`) %>% as.matrix() %>% as.character()
  #dir <- paste0("/oak/stanford/groups/engreitz/Projects/Benchmarking/",dir,"/genome_wide_predictions/")
  dir <- "/oak/stanford/groups/engreitz/Projects/Benchmarking/Predictors/LogRegClassifiers/processed/DNaseOnly_220907/genome_wide_predictions/"
  #biosampleList <- biosamples$`Biosample_id`
  biosampleList <- biosamples$`V2`
  
  enhancersPerGene <- list()
  for (biosample in biosampleList) {
    file <- paste0(dir, biosample, "_EnhancerPerGene.txt")
    if (!file.exists(file)) {
      print(paste0("Missing file for biosample: ", biosample, ": ", file))
    } else {
      x <- read.delim(file) %>% mutate(Biosample=biosample)
      colnames(x)[2] <- "nEnhancers"
      enhancersPerGene <- c(enhancersPerGene, list(x))
    }
  }	
  enhancersPerGene <- do.call(rbind, enhancersPerGene)
  return(enhancersPerGene)
}
enhancersPerGene <- loadEnhancersPerGene(predictors, biosamples, "ENCODE_E2G")



#summary tab
#

stats_tab = read.table("Encode_E2G_stats_updated.txt", header = T, sep = '\t')
gene_tab = read.table("genesPerEnhancer_ENCODE_E2G_03.2023.txt",header=T,sep = '\t')
enhancer_tab = read.table("updated_03_2023/enhancersPerGene_ENCODE_E2G.txt",header=T,sep = "\t")

#just do table and mean for everything


median_enhancers_per_genes = c()
median_genes_per_enhancers = c()


for (biosample in unique(gene_tab$Biosample)){
  bio_e_tab = enhancer_tab[enhancer_tab$Biosample == biosample,]
  bio_g_tab = gene_tab[gene_tab$Biosample == biosample,]
  median_enhancers_per_genes = c(median_enhancers_per_genes,median(bio_e_tab$nEnhancers))
  median_genes_per_enhancers = c(median_genes_per_enhancers, median(bio_g_tab$nGenes))
}

med_tab = data.frame(cbind(unique(gene_tab$Biosample),median_genes_per_enhancers,median_enhancers_per_genes))
colnames(med_tab$V1) = "biosample"

#gene_tab =read.table("genes_expressed_by_biosample.txt",header=T,sep = '\t') #had all genes, even with 0 enhancer
gene_tab_all = read.table("genes_expressed_by_biosample.txt",header=T,sep = '\t')
#gene_per_bio = data.frame(table(gene_tab$Biosample))
#colnames(gene_per_bio) = c("biosample", "total_genes_expressed")
#merge1 = merge(gene_per_bio,stats_tab)
merge2 = merge(stats_tab,med_tab)
merge3 = merge2[,c("biosample","n_unique_enhancers","median_genes_per_enhancers","median_enhancers_per_genes","nEG_links")]

write.table(merge3,"stab19.txt", quote = F,sep = '\t',row.names = F)


genes_expressed_count = read.table("genes_expressed_by_biosample_count.txt",header=T,sep = '\t') #includes with 0

mean_enhancers = c()
biosamples = c()
#for (gene in paste0(unique(genes_expressed_count$Var1))){
for (i in 1:length(unique(genes_expressed_count$Var1))){
  cat(paste0(i,","))
  gene = paste0(unique(genes_expressed_count$Var1))[i]
  expressed = as.numeric(paste0(genes_expressed_count[genes_expressed_count$Var1 == gene,"Freq"]))
  enhancers_per_gene_tab = enhancer_tab[enhancer_tab$TargetGene == paste0(gene),]
  if (nrow(enhancers_per_gene_tab) < expressed){
    enhancers = c(enhancers_per_gene_tab$nEnhancers,rep(0,expressed - nrow(enhancers_per_gene_tab)) )
  }
  else{
    enhancers = enhancers_per_gene_tab$nEnhancers
  }
  mean_enhancers = c(mean_enhancers,mean(enhancers))
  biosamples = c(biosamples,length(enhancers))
}

tab = data.frame(cbind(paste0(unique(genes_expressed_count$Var1)),mean_enhancers,biosamples))
for (i in 1:nrow(tab)){
  if (tab$mean_enhancers[i] == "NaN"){
    tab$mean_enhancers[i] = 0
  }
}

colnames(tab) = c("gene","mean_enhancers_per_biosample","nbiosamples_expressed")

write.table(tab,"stab20.txt",quote = F,row.names = F,sep = '\t')

