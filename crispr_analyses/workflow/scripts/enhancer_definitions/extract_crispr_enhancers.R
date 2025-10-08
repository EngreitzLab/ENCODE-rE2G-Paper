## Extract unique enhancer coordinates from CRISPR data (positives only) and write to .bed file

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
})

# load crispr enhancers
crispr <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# extract CRISPR enhancer-gene pairs (Regulated == TRUE) and get unique enhancers
crispr_enh <- crispr %>% 
  filter(Regulated == TRUE) %>% 
  select(chrom, chromStart, chromEnd, class = Dataset) %>% 
  distinct()

# sort according to chromosomes and coordinates
chrs <- paste0("chr", c(seq(1:23), "X", "Y"))
crispr_enh <- crispr_enh %>% 
  mutate(chrom = factor(chrom, levels = chrs)) %>% 
  arrange(chrom, chromStart)

# write to output file
write_tsv(crispr_enh, file = snakemake@output[[1]], col_names = FALSE)
