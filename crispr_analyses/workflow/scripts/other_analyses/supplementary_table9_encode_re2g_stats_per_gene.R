## Compute ENCODE-rE2G summary statistics per gene

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(BiocParallel)
})

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  message("Registering parallel backend with ", snakemake@threads, " cores.")
  register(MulticoreParam(workers = snakemake@threads))
} else {
  message("Registering serial backend.")
  register(SerialParam())
}

# Define functions ---------------------------------------------------------------------------------

# function to compute summary statistics for each gene in a prediction file
compute_stats_file <- function(file, genes) {
  
  # load predictions
  pred <- fread(file, showProgress = FALSE)
  
  # check which genes are considered expressed in the given set of predictions
  genes <- genes %>% 
    mutate(expressed = gene %in% pred$TargetGene)
  
  # compute the number of predicted enhancers per gene
  enh_per_gene <- pred %>% 
    filter(isSelfPromoter == FALSE) %>% 
    select(enhancer = name, gene = TargetGene) %>% 
    count(gene, name = "enhancers")
  
  # create output table
  output <- genes %>% 
    left_join(enh_per_gene, by = "gene") %>% 
    replace_na(list(enhancers = 0))
    
  return(output)
}

# Compute summary statistics -----------------------------------------------------------------------

# get list of thresholded predictions for each sample
pred_files <- unlist(snakemake@config$encode_re2g_predictions$thresholded)

# load gene universe and extract simple list of unique gene symbols
genes <- fread(snakemake@input$genes)
genes <- distinct(select(genes, gene = name))

# calculate number of predicted enhancers per gene across all DNase-seq experiments
gene_stats <- bplapply(pred_files, FUN = compute_stats_file, genes = genes)
gene_stats <- bind_rows(gene_stats, .id = "sample")

# calculate summary statistics per gene
summary_stats <- gene_stats %>% 
  group_by(gene) %>% 
  summarize(expressed = sum(expressed),
            mean_enh = mean(enhancers))

# reformat for supplementary table
summary_stats <- summary_stats %>% 
  select(Gene = gene, `Expressed in n DNase-seq experiments` = expressed,
         `Mean enhancers per DNase-seq experiment` = mean_enh)

# save to output file
write_csv(summary_stats, file = snakemake@output[[1]], na = "")
