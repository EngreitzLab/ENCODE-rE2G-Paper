library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(scales)
library(ggpubr)

out_dir = '/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/MYC_enhancers/'
input_file = '/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fulco2016_DHSPeaks_noDistFilter/H3K27ac.table.txt'
df = read.table(input_file, header=TRUE, sep='\t')
# col names: 
# PerturbationName	chr	start	end	enhancer.name	H3K27ac.RPM.ctrl	H3K27ac.RPM.ctrl.sem	H3K27ac.KO.RPM	H3K27ac.KO.RPM.sem	enh.krab.dist	H3K27ac.RPM.ctrl.at.PerturbationSite	H3K27ac.RPM.ctrl.sem.at.PerturbationSite	H3K27ac.RPM.KO.at.PerturbationSite	H3K27ac.RPM.KO.sem.at.PerturbationSite
# perturbation = target site (where CRISPRi was landed, enhancer = where K27ac was measured)
contact = read.table('/oak/stanford/groups/engreitz/Users/kmualim/for_maya/MYC_enhancer_coordinates.contact_vals.bed', sep='\t', header=TRUE)
# col names: chr	start	end	PerturbationName	chrom1	start1	end1	name	observedSCALENormContact
# filter input to just e1-e7+ns1 pairs (e1_cripsri, etc.)
df.filtered = dplyr::filter(df, PerturbationName!='NA', enhancer.name!='NA')
#write.table(df.filtered, paste0(out_dir, "Fulco2016_DHSPeaks_H3K27ac_MYCEnhancers.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
gex = read.table('/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/MYC_enhancers/MycSynergyScores_CRISPRiScreen.tsv', sep='\t', header=TRUE)
# col names: loc1, loc2, synergyScore

# list of enhancer coordinates
enh.coord = dplyr::select(df.filtered, PerturbationName, chr, start, end) %>% distinct()
enh.coord = enh.coord[, c(2,3,4,1)]
#write.table(enh.coord, paste0(out_dir, "MYC_enhancer_coordinates.bed"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

# heat map of delta k27ac at each enhancer/enhancer intersection (log2FC Ci/WT)
df.filtered$H3K27ac.log2FC = with(df.filtered, log2(H3K27ac.KO.RPM/H3K27ac.RPM.ctrl))
df.filtered$H3K27ac.delta.RPM = with(df.filtered, H3K27ac.KO.RPM - H3K27ac.RPM.ctrl)

# rename factors
df.filtered$PerturbationName = factor(df.filtered$PerturbationName)
levels(df.filtered$PerturbationName) <- gsub("_crispri", "", levels(df.filtered$PerturbationName))
df.filtered$enhancer.name = factor(df.filtered$enhancer.name)
levels(df.filtered$enhancer.name) <- gsub("_crispri", "", levels(df.filtered$enhancer.name))

df.filtered = dplyr::filter(df.filtered, PerturbationName!="ns1", enhancer.name!="ns1")

# log2fc heatmap
hm = ggplot(df.filtered, aes(PerturbationName, enhancer.name, fill= H3K27ac.log2FC)) + 
  geom_tile() + 
  geom_text(aes(label = round(H3K27ac.log2FC, 2)), size=3) +
  xlab('CRISPRi target') + ylab('Enhancer') +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), limits = c(-1.75, max(df.filtered$H3K27ac.log2FC)), oob = scales::squish) +
  coord_fixed(1) + 
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
  
#ggsave(plot=hm, filename=paste0(out_dir, "H3K27ac_log2FC_heatmap_v3.pdf"), height=3, width=5)

# delta k27ac heatmap
hm2 = ggplot(df.filtered, aes(PerturbationName, enhancer.name, fill= H3K27ac.delta.RPM)) + 
  geom_tile() +   
  xlab('CRISPRi target') + ylab('Enhancer') +
  scale_fill_viridis_c(limits = c(-20, 0), oob = scales::squish) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

#ggsave(plot=hm2, filename=paste0(out_dir, "H3K27ac_delta_heatmap_v2.pdf"), height=3, width=5)

# contact heatmap
# add enhancer name by joining with filtered heatmap
df.filtered$location = df.filtered$enhancer.name
locs = dplyr::select(df.filtered, chr, start, end, location)
# contact$PerturbationName = factor(contact$PerturbationName)
# levels(contact$PerturbationName) <- gsub("_crispri", "", levels(contact$PerturbationName))'
contact = dplyr::left_join(contact, locs) %>%
  dplyr::left_join(locs, by = c("chrom1"="chr", "start1"="start", "end1"="end")) %>%
  dplyr::filter(!is.na(location.x), !is.na(location.y))

getPlContact3 <- function(dist, x0=5000, scale=5.522049023040865, gamma = 1.0044322152129266) {
  #return(exp(-1*GAMMA * log(pmax(dist, x0))) / exp(-1*GAMMA * log(x0)))
  return(exp(scale + -1*gamma * log(pmax(dist, x0))) / exp(scale + -1*gamma * log(x0)))
}
contact$dist = with(contact, abs((start+end)/2-(start1+end1)/2))
contact$pl.contact = getPlContact3(contact$dist)
#contact_grouped = dplyr::filter(contact, observedSCALENormContact!=0) %>% # get rid of 0 contact values
contact_grouped = dplyr::select(contact, location.x, location.y, observedSCALENormContact, pl.contact) %>%
  unique() %>% 
  group_by(location.x, location.y) %>%
  summarise(observedSCALENormContact=mean(observedSCALENormContact), pl.contact=mean(pl.contact))
contact_grouped$super.norm.contact = contact_grouped$observedSCALENormContact/max(contact_grouped$observedSCALENormContact)

locs_dist = dplyr::select(contact, location.x, location.y, chr, start, end, start1, end1, dist) %>%
  group_by(location.x, location.y) %>%
  summarise(dist=mean(dist))
write.table(locs_dist, paste0(out_dir, "locs_and_dist.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

#  dplyr::group_by(PerturbationName, location_key) %>% 
#  summarise(observedSCALENormContact=mean(observedSCALENormContact))
#write.table(contact, paste0(out_dir, "contact_values_all.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
#write.table(contact_grouped, paste0(out_dir, "contact_values_nonzero.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')



#contact_grouped$observedSCALENormContact[contact_grouped$observedSCALENormContact==0] = 0.1
#print(contact_grouped)
hm3 = ggplot(contact_grouped, aes(location.x, location.y, fill= log10(observedSCALENormContact+0.1))) + 
  geom_tile() +   
  xlab('Location 1') + ylab('Location 2') +
  geom_text(aes(label = round(observedSCALENormContact, 0)), size=2.5) +
  #scale_fill_viridis_c(limits = c(0, log10(500)), oob = scales::squish) +
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "white", high = muted("blue"), limits = c(-0.2, log10(500)), oob = scales::squish) +
  coord_fixed(1) + 
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
#ggsave(plot=hm3, filename=paste0(out_dir, "enhancer_contact_heatmap_v2.pdf"), height=3, width=6)

## correlation between gex, hic, and k27ac
# synergyScore, observedSCALENormContact, H3K27ac.log2FC
combined = dplyr::left_join(df.filtered, contact_grouped, by=c('PerturbationName'='location.x', 'enhancer.name'='location.y')) %>%
  dplyr::left_join(gex, by=c('PerturbationName'='loc1', 'enhancer.name'='loc2')) %>%
  dplyr::filter(!is.na(synergyScore))
write.table(combined, paste0(out_dir, "combined.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

# HiC vs gex
hic.vs.gex = ggplot(combined, aes(x = log10(observedSCALENormContact + 0.01), y = synergyScore)) +
  geom_point(color="black", size=0.7, alpha=0.7) +
  stat_cor(cor.coef.name='R', aes(label=..rr.label..)) +
  labs(x ='log10(Intact Hi-C contact)', y = 'Effect on gene expression\n(log2(Observed/expected\nunder additive model)') +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) + theme_classic()
ggsave(filename=paste0(out_dir, "contact_vs_gex.pdf"), plot=hic.vs.gex, width = 3, height = 3)


# HiC vs K27ac
hic.vs.k27ac = ggplot(combined, aes(x = log10(observedSCALENormContact + 0.01), y = H3K27ac.log2FC)) +
  geom_point(color="black", size=0.7, alpha=0.7) +
  stat_cor(cor.coef.name='R', aes(label=..rr.label..)) +
  labs(x ='log10(Intact Hi-C contact)', y = 'Effect on H3K27ac\n(log2(CRISPRi/control)') +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) + theme_classic()
ggsave(filename=paste0(out_dir, "contact_vs_k27ac.pdf"), plot=hic.vs.k27ac, width = 3, height = 3)

# K27ac vs gex
k27ac.vs.gex = ggplot(combined, aes(x = H3K27ac.log2FC, y = synergyScore)) +
  geom_point(color="black", size=0.7, alpha=0.7) +
  stat_cor(cor.coef.name='R', aes(label=..rr.label..)) +
  labs(x ='Effect on H3K27ac\n(log2(CRISPRi/control)', y = 'Effect on gene expression\n(log2(Observed/expected\nunder additive model)') +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) + theme_classic()
ggsave(filename=paste0(out_dir, "k27ac_vs_gex.pdf"), plot=k27ac.vs.gex, width = 3, height = 3)


  
  
  
  


