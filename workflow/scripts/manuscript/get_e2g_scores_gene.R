## Extract ENCODE-rE2G scores for a given gene in all samples

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(BiocParallel)
  library(GenomicRanges)
})

# Define functions ---------------------------------------------------------------------------------

# function to load predictions for one sample and extract scores for one gene
extract_scores_gene <- function(file, gene, score_col) {
  
  # load predictions and extract scores for desired gene
  pred <- fread(file)
  pred <- filter(pred, TargetGene == gene)
  
  # create GRanges object with scores for all enhancers for desired gene
  scores <- pred %>% 
    filter(class != "promoter") %>% 
    select(chr = `#chr`, start, end, score = all_of(score_col))
  
  # create GRanges track from scores data frame
  if (nrow(scores) > 0) {
    scores <- makeGRangesFromDataFrame(scores, keep.extra.columns = TRUE)
  } else {
    scores <- GRanges()
  }

  return(scores)
  
}

# Extract ENCODE-rE2G scores -----------------------------------------------------------------------

# load table with all ENCODE-rE2G predictions
e2g_metadata <- fread(snakemake@input[[1]])

# get prediction files for all ENCODE-rE2G samples
e2g_files <- structure(e2g_metadata$thresholded_predictions_path, names = e2g_metadata$Sample)

## TODO: FIX METADATA TABLE
e2g_files <- e2g_files[!is.na(e2g_files)]

# get scores for the specified gene for all files
gene <- toupper(snakemake@wildcards$gene)
register(MulticoreParam(workers = 4))
pred_scores <- bplapply(e2g_files, FUN = extract_scores_gene, gene = gene, score_col = "Score")

# save pred scores to .rds file
saveRDS(GRangesList(pred_scores), file = snakemake@output[[1]])
