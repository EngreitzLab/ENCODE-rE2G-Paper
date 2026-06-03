## Plot properties distributions for Extended Data Figure 6

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(cowplot)
})

# load ENCODE-rE2G properties
combined_properties <- readRDS(snakemake@input[[1]])

# Distance to TSS distribution ---------------------------------------------------------------------

# calculate median and mean distance to TSS value
median_dist_to_tss <- median(combined_properties$distance_to_tss$distance_to_tss)
mean_dist_to_tss <- mean(combined_properties$distance_to_tss$distance_to_tss)

# re-code distance to bin data beyond 2Mb
dist_to_tss <- combined_properties$distance_to_tss %>% 
  mutate(dist_kb = distance_to_tss/1000) %>% 
  mutate(dist_kb_recoded = ifelse(dist_kb > 1000, 1025, dist_kb))

# plot distance to TSS distribution plot
p1 <- ggplot(dist_to_tss, aes(x = dist_kb_recoded)) +
  geom_histogram(breaks = c(seq(0, 1000, 50), 1050), fill = "#9b9b9b", color = "black") +
  geom_vline(xintercept = median_dist_to_tss / 1000, linetype = "dashed", color = "red") +
  labs(x = "Distance of E-G pairs to TSS (kb)", y = "# Enhancer-Gene Pairs") +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1050),
                     labels = c(0, 250, 500, 750, 1000, "1000-5000")) +
  theme_classic()

# calculate number of enhancer-gene interactions below 10 and 100kb (for main text)
below_10kb <- mean(combined_properties$distance_to_tss$distance_to_tss < 10000)
below_100kb <- mean(combined_properties$distance_to_tss$distance_to_tss < 100000)

# Number of genes per enhancer ---------------------------------------------------------------------

# calculate median and mean genes per enhancer
median_genes_per_enh <- median(combined_properties$genes_per_enh$target_genes)
mean_genes_per_enh <- mean(combined_properties$genes_per_enh$target_genes)

# re-code genes per enhancers to bin data beyond 10
genes_per_enh <- combined_properties$genes_per_enh %>% 
  mutate(target_genes_recoded = ifelse(target_genes > 10, 10.5, target_genes))

# make histogram plot
p2 <- ggplot(genes_per_enh, aes(x = target_genes_recoded)) +
  geom_histogram(breaks = c(0:10, 11), fill = "#9b9b9b", color = "black") +
  geom_vline(xintercept = median_genes_per_enh, linetype = "dashed", color = "red") +
  labs(x = "# Genes per Enhancer", y = "# Enhancers") +
  scale_x_continuous(breaks = 0:11, labels = c(0:10, "11-33")) +
  theme_classic()

# Number of enhancers per gene ---------------------------------------------------------------------

# calculate median and mean enhancers per gene
median_enh_per_gene <- median(combined_properties$enhancers_per_gene$enhancers)
mean_enh_per_gene <- mean(combined_properties$enhancers_per_gene$enhancers)

# re-code enhancer per gene to bin data beyond 10
enh_per_gene <- combined_properties$enhancers_per_gene %>% 
  mutate(enhancers_recoded = ifelse(enhancers > 40, 40.5, enhancers))

# Number of genes per enhancer
p3 <- ggplot(enh_per_gene, aes(x = enhancers_recoded)) +
  geom_histogram(breaks = c(0:40, 41), fill = "#9b9b9b", color = "black") +
  geom_vline(xintercept = median_enh_per_gene, linetype = "dashed", color = "red") +
  labs(x = "# Enhancers per Gene", y = "# Genes") +
  scale_x_continuous(breaks = c(seq(0, 40, 10), 41), labels = c(seq(0, 40, 10), "41-679")) +
  theme_classic()

# Create plot --------------------------------------------------------------------------------------

# save plots to files
fig <- plot_grid(p1, p2, p3, ncol = 1, labels = c("a", "b", "c"))
ggsave(fig, file = snakemake@output[[1]], height = 9, width = 4)
