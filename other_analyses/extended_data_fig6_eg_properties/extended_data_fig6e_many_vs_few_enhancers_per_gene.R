#cd /oak/stanford/groups/engreitz/Users/ejagoda/encode_e2g_revision
#cp /oak/stanford/groups/engreitz/Users/kmualim/distalreg_figures/EnhancerPerGene_all.encodee2g.tsv.gz ./
#conda activate /oak/stanford/groups/engreitz/Users/abaskara/.conda/envs/R_env

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


#want to find, of genes expressed in at least half of the biosamples >=729, if enrichment bad, try 3/4 didn't reduce that much
#and then greater than 5 less than 5 

tab = read.table("EnhancerPerGene_all.encodee2g.tsv.gz",header = T)
#first get the cell frequency of the genes

gene_freq_tab = data.frame(table(tab$TargetGene))
g_50 = gene_freq_tab[gene_freq_tab$Freq >= 729,]

tab_use = tab[tab$TargetGene %in% g_50$Var1,] #14460 out of 20225

mean_enhancers_per_gene <- tab_use %>% 
  group_by(TargetGene) %>% 
  summarize(nEnhancers=mean(X0)) %>%
  as.data.frame()

write.table(mean_enhancers_per_gene,"mean_enhancers_per_gene_g729_biosamples.txt",quote = F,sep = '\t',row.names = F)
write.table(mean_enhancers_per_gene[mean_enhancers_per_gene$nEnhancers >= 5,"TargetGene"], "g5_enhancers_per_gene_g729_biosamples.txt", quote = F,sep = '\t',row.names = F)
write.table(mean_enhancers_per_gene[mean_enhancers_per_gene$nEnhancers < 1.5 ,"TargetGene"], "l1.5_enhancers_per_gene_g729_biosamples.txt", quote = F,sep = '\t',row.names = F)

##now on local

library(ggplot2)
mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10))

david_many_function = read.table("david_many_chart_F3A305911AF41734109331412_UP_KW_BIOLOGICAL_PROCESS .txt",header=T,sep = '\t')
david_few_function = read.table("david_few_chart_F3A305911AF41734109331412_UP_KW_BIOLOGICAL_PROCESS .txt",header=T,sep = '\t')
david_many = read.table("g5_enhancers_per_gene_g729_biosamples.txt",header=T,sep = '\t',fill = T)
david_few = read.table("l1.5_enhancers_per_gene_g729_biosamples.txt",header=T,sep = '\t',fill = T)



p <- ggplot(david_many_function, aes(reorder(Term,-Benjamini) , -log10(Benjamini))) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10)) +
  geom_bar(stat = 'identity', aes(fill = -log10(Benjamini)))
print(p) 

p1 <- ggplot(david_many_function[david_many_function$Benjamini < 0.05,], aes(reorder(Term,Fold.Enrichment) , Fold.Enrichment)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10))+
  labs(title = "KW_BIOLOGICAL_PROCESS\nMany Genes\nP-value adjusted < 0.05")
#print(p1) 

p2 <- ggplot(david_few_function[david_few_function$Benjamini < 0.05,], aes(reorder(Term,Fold.Enrichment) , Fold.Enrichment)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10))+
  labs(title = "KW_BIOLOGICAL_PROCESS\nFew Genes\nP-value adjusted < 0.05")
#print(p2) 


p3 <- ggplot(david_many_function[1:5,], aes(reorder(Term,Fold.Enrichment) , Fold.Enrichment)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10))+
  labs(title = "MANY")+
  mytheme
#print(p3) 

p4 <- ggplot(david_few_function[1:5,], aes(reorder(Term,Fold.Enrichment) , Fold.Enrichment)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10))+
  labs(title = "few")+
  mytheme
#print(p4) 

ggsave("all_David_plots.pdf",width = 10,height = 5)
grid.arrange(p1,p2,p3,p4,nrow = 2)
dev.off()
#maybe do volcanos
ggsave("all_David_plots.eps")
grid.arrange(p1,p2,p3,p4,nrow = 2)
dev.off()

ggsave(file = "Many_David_functional_sig_only_plot.svg",plot = p1,width = 10,height = 5)
ggsave(file = "Few_David_functional_sig_only_plot.svg",plot = p2,width = 10,height = 5)
ggsave(file = "Both David functional sig onlys.svg",plot = grid.arrange(p1,p2,nrow = 2))

david_many_function$type = "many"
david_few_function$type = "few"

david_few_function$Fold.Enrichment_for_combined = david_few_function$Fold.Enrichment * -1
david_many_function$Fold.Enrichment_for_combined = david_many_function$Fold.Enrichment

sig_few = david_few_function[david_few_function$Benjamini < 0.05,]
sig_many = david_many_function[david_many_function$Benjamini < 0.05,]
combined = rbind(sig_many[1:5,],sig_few[1:5,])
###okay maybe remove transcription from the text since it's both
#combined = rbind(david_many_function[david_many_function$Benjamini < 0.05,],sig_few[1:5,])
p = ggplot(combined, aes(reorder(Term,Fold.Enrichment_for_combined) , Fold.Enrichment_for_combined,fill = type)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + 
  theme(axis.text.y=element_text(angle=0, hjust=1, size = 10)) +
  geom_vline(xintercept = -1,linetype = "dashed")+
  mytheme

ggsave("david_many_few_functional_for_fig.pdf",width = 7,height = 4)
print(p)
dev.off()

ggsave(file = "david_many_few_functional_for_fig.pdf",plot = p)

tab = read.table("mean_enhancers_per_gene_g729_biosamples.txt",header = T,sep = '\t')
head(tab[order(tab$nEnhancers,decreasing = T),])
#TargetGene nEnhancers
#10963       SHOX  118.70270
#2755      CSF2RA   93.81273
#4894      GPR143   84.10092
#3648   ELFN1-AS1   82.32494
#9562      PSAPL1   74.72422
#6625   LINC01819   66.70975
##for locus plots
#top 5 many = 

#TargetGene nEnhancers
#3263      DNAJC16   1.005487
#11460       SNX15   1.004121
#7548       MRPS12   1.004118
#8453       ORMDL2   1.003432
#12646        TMX2   1.001374
#12647 TMX2-CTNND1   1.001374
