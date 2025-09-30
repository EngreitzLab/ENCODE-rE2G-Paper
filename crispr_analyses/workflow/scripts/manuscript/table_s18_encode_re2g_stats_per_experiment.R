## Compute ENCODE-rE2G summary statistics per DNase-seq experiment

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

# function to compute summary statistics for one prediction file
compute_stats_file <- function(file) {
  
  # load predictions
  pred <- fread(file, showProgress = FALSE)
  
  # get eg links without self promoters
  eg_links <- pred %>% 
    filter(isSelfPromoter == FALSE) %>% 
    select(enhancer = name, gene = TargetGene)
  
  # calculate the number of unique enhancers and total E-G links
  n_enh <- n_distinct(eg_links$enhancer)
  n_links <- nrow(eg_links)
  
  # calculate median number of predicted enhancers per gene and target genes per enhancers
  genes_per_enh <- count(eg_links, enhancer, name = "target_genes")
  enh_per_gene <- count(eg_links, gene, name = "enhancers")
  median_genes_per_enh <- median(genes_per_enh$target_genes)
  median_enh_per_gene <- median(enh_per_gene$enhancers)
  
  # create output table
  output <- tibble(
    `Unique enhancers` = n_enh,
    `Predicted E-G interactions` = n_links,
    `Median genes per enhancer` = median_genes_per_enh,
    `Median enhancers per gene` = median_enh_per_gene
  )
  
  return(output)
  
}

# Compute summary statistics -----------------------------------------------------------------------

# get list of thresholded predictions for each sample
pred_files <- unlist(snakemake@config$encode_re2g_predictions$thresholded)

# compute summary statistics for each prediction file
summary_stats <- bplapply(pred_files, FUN = compute_stats_file)

# combine into one table
summary_stats <- bind_rows(summary_stats, .id = "Sample")

# save to output file
write_csv(summary_stats, file = snakemake@output[[1]], na = "")
