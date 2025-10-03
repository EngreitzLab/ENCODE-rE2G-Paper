library(data.table)
library(ggplot2)
library(dplyr)

base.dir <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Delaneau2019-chromQTL/data/"
out.dir <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Delaneau2019-chromQTL/analysis/"

snp <- fread(file.path(base.dir, "SNPs.txt.gz"))
colnames(snp) <- c("variant.chr","variant.position","variant.id","ref","alt")

qtl <- fread(file.path(base.dir, "LCL.cQTLs.txt.gz"))
colnames(qtl) <- c("qtl.id","variant.id","pval","beta","pval.adj")
qtl$histone.mark <- sapply(qtl$qtl.id, function(s) {strsplit(s, "_")[[1]][[1]]}, USE.NAMES = F)

regions <- fread(file.path(base.dir, "LCL.CRD.bed.gz"))
colnames(regions) <- c("region.chr","region.start","region.end","qtl.id","tree.id","parent1.id","parent2.id","isInternalNode")

all.data <- merge(qtl, snp, all.x = T, all.y = F) # by = 'variant.id' -> merge QTLs with their respective SNPs if posssible
all.data <- merge(all.data, regions, all.x = T, all.y = T, by = "qtl.id") # merge variants wih regions if have same qtl.id (variant affects this peak)

all.data$distance <- with(all.data, abs((region.end + region.start)/2 - variant.position)) # distance btw variant and center of region
all.data$variant.overlaps.region <- with(all.data, variant.chr == region.chr & 
                                          variant.position > region.start & variant.position < region.end) # does variant fall in region?

###############
# Write annotation files -- how do we these are all K27ac regions?
write.table(regions[, c("region.chr","region.start","region.end","qtl.id")], 
            file.path(base.dir, "processed/K27ac.regions.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(regions[regions$region.end - regions$region.start < 5000, c("region.chr","region.start","region.end","qtl.id")], 
            file.path(base.dir, "processed/K27ac.regions5kb.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(regions[regions$region.end - regions$region.start < 1000, c("region.chr","region.start","region.end","qtl.id")], 
            file.path(base.dir, "processed/K27ac.regions1kb.bed"), quote = F, sep = "\t", col.names = F, row.names = F)

write.table(subset(all.data, histone.mark == "H3K27ac")[, c("variant.chr","variant.position","variant.position","region.chr","region.start","region.end")],
            file.path(base.dir, "processed/K27ac.qtl.bedpe"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(subset(all.data, histone.mark == "H3K27ac" & !is.na(pval))[, c("variant.chr","variant.position","variant.position")],
            file.path(base.dir, "processed/K27ac.qtl.variants.bed"), quote = F, sep = "\t", col.names = F, row.names = F)

###############
# Diagnostic Plots
# Volcano plot
pdf(file.path(out.dir, "diagnostics/volcano_pval_vs_beta.pdf"))
with(subset(qtl, histone.mark == "H3K27ac"), plot(beta, -1*log10(pval), 
                                                  cex.axis = 1.5, cex.lab = 1.5, main = 'Volcano'))
dev.off()

#Region widths
pdf(file.path(out.dir, "diagnostics/cdf_region_widths.pdf"))
plot(ecdf(with(subset(all.data, histone.mark == "H3K27ac"), region.end - region.start)), 
     xlab = "Region width", ylab = "CDF", main = "Region Widths", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#CDF of distances
pdf(file.path(out.dir, "diagnostics/cdf_distances.pdf"))
plot(ecdf(log10(subset(all.data, histone.mark == "H3K27ac")$distance + 1)), 
     xlab = "distance (log10)", ylab = "CDF", main = "QTL distances")
dev.off()

pdf(file.path(out.dir, "diagnostics/beta_vs_distance.pdf"))
with(subset(all.data, histone.mark == "H3K27ac"), plot(distance, beta,
     cex.axis = 1.5, cex.lab = 1.5))
dev.off()

pdf(file.path(out.dir, "diagnostics/region_size_density.pdf"))
ggplot(subset(all.data, histone.mark == "H3K27ac"),
       aes(x = region.end - region.start,
           color = isInternalNode)) + 
  geom_line(stat = 'density')
dev.off()

###############
# make dataset of variants that are self-qtl and qtl for other peak
self.qtl <- subset(all.data, histone.mark == "H3K27ac" & !is.na(pval) & variant.overlaps.region & (region.end - region.start < 5e3))
freq.table <- as.data.table(table(subset(all.data, histone.mark == "H3K27ac" & !is.na(pval))$variant.id))
multiple.qtl <- subset(freq.table, N > 1)$V1
multiple.qtl2 <- subset(all.data, histone.mark == "H3K27ac" & !is.na(pval) & variant.id %in% self.qtl$variant.id & variant.id %in% multiple.qtl)
write.table(multiple.qtl2[, c("variant.chr","variant.position","variant.position","region.chr","region.start","region.end")],
            file.path(base.dir, "processed/K27ac.qtl.selfandmulti.bedpe"), quote = F, sep = "\t", col.names = F, row.names = F)

nrow(subset(all.data, histone.mark == "H3K27ac" & !is.na(pval))) #Number of H3K27ac qtl
nrow(subset(all.data, histone.mark == "H3K27ac" & !is.na(pval) & (region.end - region.start < 5e3)))
nrow(self.qtl) #Number of qtl in self-peak
sum(self.qtl$variant.id %in% multiple.qtl2$variant.id) #Number of self-qtl that are qtl for another peak

#
self.qtl2 <- self.qtl
colnames(self.qtl2) <- paste0(colnames(self.qtl2), ".of.self.region")

multiple.qtl3 <- merge(multiple.qtl2, self.qtl2, by.x = "variant.id", by.y = "variant.id.of.self.region", all.x =TRUE, all.y = FALSE)

pdf(file.path(out.dir, "self-qtls_distance_vs_beta.pdf"))
ggplot(multiple.qtl3,
       aes(x = log10(distance - distance.of.self.region + 1),
           y = beta, 
           group = variant.id,
           color = variant.id == "rs117157251")) + 
  geom_point() + 
  geom_line() + 
  labs(x = "Distance from variant to peak (log10)",
       y = "Effect Size",
       title = "Variants which are cQTL for >1 peak") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20)) + 
  scale_color_manual(values=c("#000000","#FF3333")) + 
  theme(legend.position = "none")
dev.off()

# ggplot(subset(multiple.qtl3, variant.id == "rs117157251"),
#        aes(x = log10(distance - distance.of.self.region + 1),
#            y = beta, #/beta.of.self.region * sign(beta)
#            group = variant.id,
#            color = variant.id == "rs117157251")) + 
#   geom_point() + 
#   #geom_line(data = multiple.qtl3, alpha = .75) + 
#   labs(x = "Distance from variant to peak (log10)",
#        y = "Effect Size",
#        title = "Variants which are cQTL for >1 peak") + 
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title = element_text(size = 20)) + 
#   scale_color_manual(values=c("#000000","#FF3333")) + 
#   theme(legend.position = "none") + 
#   geom_line(data = subset(multiple.qtl3, variant.id == "rs117157251"), size = 1.2, col = 'red')+
#   geom_point(data = subset(multiple.qtl3, variant.id == "rs117157251"), size = 3, col = 'red') + 
#   coord_cartesian(xlim = c(0, 5), ylim = c(-1.5, 2.2))

pdf(file.path(out.dir, "beta_vs_distance_slope_distribution.pdf"))
ggplot(multiple.qtl2[, summary(lm(abs(beta) ~ log10(distance)))$coefficients[2,1], by = 'variant.id'],
       aes(x = V1)) + 
  geom_histogram() + 
  labs(x = "Slope", y = 'Count', title = "Distribution of Slopes") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
dev.off()

pdf(file.path(out.dir, "self-qtls-not-mult_beta_vs_distance.pdf"))
ggplot(subset(all.data, histone.mark == "H3K27ac" & !is.na(pval) & variant.id %in% self.qtl$variant.id & !(variant.id %in% multiple.qtl)),
       aes(x = log10(distance),
           y = beta, 
           group = variant.id)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "Distance from QTL to Peak (log10)",
       y = "Beta",
       title = "Effect sizes of variant over its peaks") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
dev.off()

## Maya's plot(s)
## use multiple.qtl3 dataset (in peak it regulates and linked to other peaks)

## iterations on beta
multiple.qtl3$abs.beta = abs(multiple.qtl3$beta)
multiple.qtl3$beta.norm = with(multiple.qtl3, beta/beta.of.self.region)

## iterations on distance
#multiple.qtl3$log.dist.to.peak = with(multiple.qtl3, (log10(distance - distance.of.self.region + 1)))
multiple.qtl3$dist.to.peak = with(multiple.qtl3, (distance-distance.of.self.region))
#multiple.qtl3 = mutate(multiple.qtl3, binned.dist = cut_width(log.dist.to.peak, width=0.5, boundary=0))
# filter to distance > 1kb
multiple.qtl3 = dplyr::filter(multiple.qtl3, dist.to.peak>1000)

multiple.qtl3$binned.dist = ""
multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak<=10000] = "1-10kb"
n1 = filter(multiple.qtl3, binned.dist=="1-10kb") %>% nrow()
multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak<=10000] = paste0("1-10kb\nN=", n1)

multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak>10000 & multiple.qtl3$dist.to.peak<=100000] = "10-100kb"
n2 = filter(multiple.qtl3, binned.dist=="10-100kb") %>% nrow()
multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak>10000 & multiple.qtl3$dist.to.peak<=100000] = paste0("10-100kb\nN=", n2)


multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak>100000] = ">100kb"
n3 = filter(multiple.qtl3, binned.dist==">100kb") %>% nrow()
multiple.qtl3$binned.dist[multiple.qtl3$dist.to.peak>100000] = paste0(">100kb\nN=", n3)

multiple.qtl3$binned.dist = factor(multiple.qtl3$binned.dist)
multiple.qtl3$binned.dist = ordered(multiple.qtl3$binned.dist, levels = c(paste0("1-10kb\nN=", n1), 
                                                                          paste0("10-100kb\nN=", n2), 
                                                                          paste0(">100kb\nN=", n3)))

#pdf(file.path(out.dir, "test_fig_v1.pdf"), width=3.25, height=7)
fig = ggplot(multiple.qtl3, aes(x = binned.dist,
                          y = -beta.norm)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1), size=0.5, alpha=0.5) +
  geom_hline(yintercept=0, lty='solid') +
  geom_hline(yintercept=-1, lty='dashed', size=0.5, color='red') +
  labs(x="Distance from\nvariant to peak", y="Effect on nearby peak\n(relative to effect on primary peak)") + 
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave(file.path(out.dir, "main_fig4d.eps"), fig, height=3, width=2)
ggsave(file.path(out.dir, "main_fig4d.pdf"), fig, height=3, width=2)


## Wilcoxon rank sum test
g1 = filter(multiple.qtl3, binned.dist== paste0("1-10kb\nN=", n1))
g2 = filter(multiple.qtl3, binned.dist==paste0("10-100kb\nN=", n2))
g3 = filter(multiple.qtl3, binned.dist==paste0(">100kb\nN=", n3))

g1.v.g3 = wilcox.test(g1$beta.norm, g3$beta.norm)
g1.v.g2 = wilcox.test(g1$beta.norm, g2$beta.norm)
g2.v.g3 = wilcox.test(g2$beta.norm, g3$beta.norm)

print("Wilcoxon\n1 vs 3")
print(g1.v.g3)

print("1 vs 2")
print(g1.v.g2)

print("2 vs 3")
print(g2.v.g3)

## t-test
t1 = t.test(g1$beta.norm, mu=0)
t2 = t.test(g2$beta.norm, mu=0)
t3 = t.test(g3$beta.norm, mu=0)

print("Bin 1")
print(t1)
print("Bin 2")
print(t2)
print("Bin 3")
print(t3)
  
  
  



dev.off()