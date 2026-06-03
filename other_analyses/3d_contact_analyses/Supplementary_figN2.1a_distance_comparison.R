library(dplyr)
library(readr)
library(ggplot2)

My_Theme2 = theme_bw()+theme(
  plot.title = element_text(size=18, face="bold",hjust = 0.5),
  axis.title = element_text(size = 18),
  axis.text.x = element_text(size = 18,margin = margin(t = 0, r = 0, b = 10, l = 10)),
  axis.text.y = element_text(size = 18,margin = margin(t = 0, r = 0, b = 10, l = 10)),
  legend.text = element_text(size = 18),
  legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

setwd("/Users/wangxi/Desktop/CIA/")

# load CRISPR data
crispr_url <- "https://github.com/EngreitzLab/CRISPR_comparison/raw/refs/heads/main/resources/crispr_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
crispr <- read_tsv(crispr_url, show_col_types = FALSE)

# load gene (TSS) universe file
tss_url <- "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/refs/heads/main/reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed"
tss <- read_tsv(tss_url, show_col_types = FALSE)

# filter CRISR data for valid element-gene interactions and genes in our "gene universe"
crispr <- crispr %>% 
  filter(ValidConnection == "TRUE") %>% 
  filter(measuredGeneSymbol %in% tss$name)

# this will leave 471 positives and 9885 negatives
table(crispr$Regulated)

df_prom <- data.frame(read.table('RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.hg38.bed',header=F))

df_crispri.distance <- crispr %>% filter(Regulated==TRUE) %>%
  merge(df_prom, by.x=c("measuredGeneSymbol"), by.y=c("V4")) %>%
  transmute(variable='CRISPRi', distance=abs(chromStart -(V2+V3)/2))

df_gwas.distance <- data.frame(read.table('FeatureSet/UKBiobank.ABCGene.anyabc_2.tsv',sep='\t',header=T)) %>%
  filter(truth=='TRUE') %>%
  merge(df_prom, by.x=c("TargetGene"), by.y=c("V4")) %>%
  transmute(variable='GWAS', distance=PromoterDistanceToBestSNP)

df_eqtl.distance <- data.frame(read.table('FeatureSet/GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv',sep='\t',header=F)) %>%
  `names<-`(c('chr','start','end','id','celltype','TargetGene','s1','s2')) %>%    
  merge(df_prom, by.x=c("TargetGene"), by.y=c("V4")) %>%
  transmute(variable='eQTL', distance=abs(start -(V2+V3)/2))

df_credset.distance <- rbind(df_crispri.distance, df_gwas.distance, df_eqtl.distance) 

df_credset.n <- df_credset.distance %>% 
  group_by(variable) %>%
  summarize(n=n()) %>% 
  mutate(name=paste0(variable, ' (n=', n,')'))

p<-rbind(df_crispri.distance, df_gwas.distance, df_eqtl.distance) %>%
  ggplot(aes(x=distance/1000, color=variable)) + 
  #geom_density(cex=1.2,alpha=0.2,show_guide = FALSE) +
  stat_density(cex=1.4, alpha=0.7, geom="line", position="identity") +
  My_Theme2 + 
  xlab('distance to TSS (kb)') + ylab('density') + xlim(0,1000) +
  scale_color_discrete(labels=df_credset.n$name)

p

ggsave(file=paste0("/Users/wangxi/Desktop/CIA/", "fig3a.svg"), plot=p, width=7, height=4)
