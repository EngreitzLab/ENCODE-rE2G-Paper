## Plot DNase-seq signal of predicted ENCODE-rE2G E-G interactions as function of distance to TSS

# save.image("encode_re2g_activity_vs_dist.rda")
# stop()

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggrastr)
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

# extract activity, 3D contact, distance to TSS and E2G score only from thresholded predictions
extract_cols <- function(pred_file, enh_file) {
  
  # load prediction and enhancer files
  pred <- fread(pred_file, showProgress = FALSE)
  enh <- fread(enh_file, showProgress = FALSE)
  
  # remove any self promoters from predictions
  pred <- filter(pred, isSelfPromoter == FALSE)
  
  # extract relevant columns only from predictions file
  pred <- pred %>% 
    select(enhancer = name, TargetGene, score = Score, distToTSS = `distanceToTSS.Feature`,
           contact = `3DContact.Feature`)
  
  # extract enhancer activities from enhancer list file
  enh <- enh %>% 
    mutate(enhancer = paste0(chr, ":", start, "-", end)) %>% 
    select(enhancer, activity_base)
  
  # merge together into one output table
  output <- merge(pred, enh, by = "enhancer", all.x = TRUE)
  
  return(output)
  
}

# Compute summary statistics -----------------------------------------------------------------------

# get list of thresholded predictions for each sample
pred_files <- unlist(snakemake@config$encode_re2g_predictions$thresholded)

# get enhancer list files containing DNase-seq activity for all samples (and ensure correct order)
enh_files <- unlist(snakemake@config$encode_re2g_enhancer_lists)
enh_files <- enh_files[names(pred_files)]

# get activity, contact and distance for all samples
activities <- bpmapply(FUN = extract_cols, pred_files, enh_files, SIMPLIFY = FALSE)

# combine into one table
activities <- rbindlist(activities, idcol = "sample")

# bin pairs by distance to TSS
dist_bins <- c(0, 10, 50, 100, 5000)
activities <- activities %>%
  mutate(distToTSS = distToTSS / 1000) %>% 
  mutate(dist_bin = cut(distToTSS, breaks = dist_bins)) %>% 
  filter(!is.na(dist_bin))

# make violin + boxplots
p <- ggplot(activities, aes(x = dist_bin, y = activity_base)) +
  geom_violin(fill = "gray") +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.1) +
  labs(x = "Distance to TSS (kb)", y = "Enhancer activity\n(DNase-seq signal)",
       title = "Enhancer activity vs. distance to TSS\nENCODE-rE2G predicted E-G interactions") +
  theme_classic()

# save plot to output file
ggsave(p, filename = snakemake@output[[1]], width = 4, height = 3.75)
