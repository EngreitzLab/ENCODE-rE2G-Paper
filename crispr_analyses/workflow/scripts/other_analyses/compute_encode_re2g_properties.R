# Calculate properties for Extended Data Figure 6

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# input files
pred_files <- unlist(snakemake@input)
names(pred_files) <- sub("encode_e2g_predictions_(.+)_thresholded_predictions.tsv.gz", "\\1",
                         basename(pred_files))

# function to compute or extract properties from ENCODE-rE2G files
compute_properties <- function(pred_file) {
  
  message("Processing: ", pred_file)
  pred <- fread(pred_file)
  pred <- filter(pred, isSelfPromoter == FALSE)
  
  # compute number of genes per enhancer
  genes_per_enh <- pred %>% 
    count(name, name = "target_genes") %>% 
    rename(enhancer = name)
  
  # compute number of enhancers per gene
  enhancers_per_gene <- count(pred, TargetGene, name = "enhancers")
  
  # extract distance to TSS for each E-G pair
  distance_to_tss <- pred %>% 
    unite(col = "eg_pair", TargetGene, name, sep = "|") %>% 
    select(eg_pair, distance_to_tss = distanceToTSS.Feature) %>% 
    as.data.table()
  
  # combine into output list
  output <- list(
    genes_per_enh = genes_per_enh,
    enhancers_per_gene = enhancers_per_gene,
    distance_to_tss = distance_to_tss
  )

  return(output)
  
}

# compute properties for all predictions
properties <- lapply(pred_files, FUN = compute_properties)

# combine properties tables across all predictions
combined_properties <- properties %>%
  transpose() %>%
  imap(~bind_rows(.x, .id = "sample"))

# save to output file
saveRDS(combined_properties, file = snakemake@output[[1]])
