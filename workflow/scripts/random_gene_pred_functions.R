## Functions to compute random gene within X distance baseline predictor

# function to compute performance of random gene predictor and create row to add to PR table
computeRandomGenePredPR <- function(merged, tss_annot, max_distance = 1e6) {
  
  # compute precision and recall for random predictor
  random_pr <- computeRandomGenePred(merged, tss_annot = tss_annot, max_distance = max_distance)
  
  # compute average precision and recall for PRC curve
  output <- random_pr %>% 
    summarize(precision = mean(precision, na.rm = TRUE), recall = mean(recall, na.rm = TRUE)) %>% 
    mutate(pred_uid = "baseline.randomExprGene", alpha = 1, F1 = NaN) %>% 
    relocate(pred_uid, alpha, precision, recall, F1)
  
  return(output)
  
}

# function to add random gene baseline predictor entry to pred_config file
addRandomGenePredConfig <- function(pred_config, distance = "1Mb", color = "black") {
  
  # pred_config entry for random gene baseline predictor
  random_gene_config <- data.table(
    pred_id = "baseline",
    pred_col = "randomExprGene",
    boolean = TRUE,
    alpha = NA,
    aggregate_function = "max",
    fill_value = 0,
    inverse_predictor = FALSE,
    pred_name_long = paste0("Random expr. gene (", distance, ")"),
    color = color,
    plot_crispr = TRUE,
    pred_uid = "baseline.randomExprGene"
  )
  
  # add random gene predictor to pred config
  pred_config <- bind_rows(pred_config, random_gene_config)
  
  return(pred_config)
  
}

# function to compute precision and recall for a predictor assigning each enhancer to exactly one
# gene (TSS) within a certain distance from the enhancer
computeRandomGenePred <- function(merged, tss_annot, max_distance = 1e6) {
  
  # add unique enhancer id to merged data
  merged <- mutate(merged, enh_id = factor(paste0(chrom, ":", chromStart, "-", chromEnd)))
  
  # create 1Mb windows for each TSS
  tss_windows <- resize(tss_annot, width = max_distance, fix = "center")
  
  # extract CRISPR elements and create GenomicRanges object
  expt_enh <- merged %>% 
    select(chrom, chromStart, chromEnd, enh_id) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # get expressed genes within specified distance of each element
  ovl <- findOverlaps(expt_enh, tss_windows)
  enh_tss_dist <- tibble(enh = expt_enh$enh_id[queryHits(ovl)],
                         gene = tss_windows$name[subjectHits(ovl)])
  
  # get any experimental positive genes for each element
  positives <- merged %>% 
    select(chrom, chromStart, chromEnd, name, measuredGeneSymbol, Regulated, enh_id) %>% 
    distinct() %>% 
    filter(Regulated == "TRUE")
  
  # convert tss within distance and positives to list of vectors containing genes per element
  enh_tss_dist <- split(enh_tss_dist$gene, f = enh_tss_dist$enh)
  positives <- split(positives$measuredGeneSymbol, f = positives$enh_id)
  positives <- positives[names(enh_tss_dist)]  # anxiety driven sorting
  
  # compute theoretical precision and recall of a random gene predictor for each enhancer
  random_pr <- mapply(random_gene_pr, enh_tss_dist, positives, SIMPLIFY = FALSE)
  random_pr <- bind_rows(random_pr, .id = "element")
  
  return(random_pr)
  
}

# function to compute precision and recall for a random predictor
random_gene_pr <- function(genes, pos) {
  
  # compute the proportion of positives in expressed genes, which corresponds to random precision
  if (length(pos) > 0) {
    if (length(genes) > 0) {
      precision <- mean(genes %in% pos)
    } else {
      precision <- 0
    }
  } else {
    precision <- NA
  }
  
  # compute recall, which corresponds to random precision / number of positives (NA if no positives
  # since computing a recall value for these enhancers doesn't make sense)
  if (length(pos) > 0) {
    recall <- precision / length(pos)
  } else {
    recall <- NA
  }
  
  # create output
  output <- tibble(precision, recall)
  
  return(output)
  
}

# Empirically simulate random gene predictor -------------------------------------------------------

# simulate random gene predictor n times and calculate mean precision and recall across iterations
simulateRandomGenePred <- function(merged, tss_annot, max_distance = 1e6, n = 100) {
  
  # add unique enhancer id to merged data
  merged <- mutate(merged, enh_id = factor(paste0(chrom, ":", chromStart, "-", chromEnd)))
  
  # create 1Mb windows for each TSS
  tss_windows <- resize(tss_annot, width = max_distance, fix = "center")
  
  # extract CRISPR elements and create GenomicRanges object
  expt_enh <- merged %>% 
    select(chrom, chromStart, chromEnd, enh_id) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # get expressed genes within specified distance of each element
  ovl <- findOverlaps(expt_enh, tss_windows)
  enh_tss_dist <- tibble(enh = expt_enh$enh_id[queryHits(ovl)],
                         gene = tss_windows$name[subjectHits(ovl)])
  
  # extract CRISPR E-G pairs including experimental outcome
  expt <- merged %>% 
    select(enh_id, measuredGeneSymbol, Regulated) %>% 
    distinct()
  
  # simulate random gene predictor n times
  sim_pred <- replicate(n, expr = simulate_random_gene_enh(expt, enh_tss_dist), simplify = FALSE)
  
  # combine into one data frame and calculate mean precision and recall across iterations
  output <- bind_rows(sim_pred, .id = "iteration") %>% 
    summarize(pred_uid = "baseline.randomExprGene", alpha = 1, precision = mean(precision),
              recall = mean(recall), F1 = NaN)

  return(output)
  
}

# function to simulate random gene predictor and extract precision and recall
simulate_random_gene_enh <- function(expt, enh_tss_dist) {
  
  # randomly pick one gene within distance per enhancer
  random_gene <- enh_tss_dist %>% 
    group_by(enh) %>% 
    sample_n(size = 1) %>% 
    mutate(pred_value = TRUE)
  
  # merge with CRISPR data
  expt <- expt %>% 
    left_join(random_gene, by = c("enh_id" = "enh", "measuredGeneSymbol" = "gene")) %>% 
    replace_na(replace = list(pred_value = FALSE))
  
  # make contingency table
  tab <- as.matrix(table(expt$pred_value, expt$Regulated))
  
  # extract true positives, false positives and false negatives
  tp <- tab["TRUE", "TRUE"]
  fp <- tab["TRUE", "FALSE"]
  fn <- tab["FALSE", "TRUE"]
  
  # calculate precision and recall
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  
  # create output data frame
  output <- tibble(precision, recall)
  
  return(output)
  
}
