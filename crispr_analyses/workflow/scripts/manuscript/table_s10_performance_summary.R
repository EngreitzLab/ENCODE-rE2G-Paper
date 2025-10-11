## Create table of performance summary for (almost) all predictors

suppressPackageStartupMessages({
  
  # attach required packages
  library(tidyverse)
  library(ROCR)
  library(caTools)
  
  # load required functions
  crispr_benchmark_dir <- snakemake@config$crispr_benchmark_dir
  source(file.path(crispr_benchmark_dir, "workflow/scripts/crisprComparisonLoadInputData.R"))
  source(file.path(crispr_benchmark_dir, "workflow/scripts/crisprComparisonPlotFunctions.R"))
  source(file.path(crispr_benchmark_dir, "workflow/scripts/crisprComparisonBootstrapFunctions.R"))
  
})

# set seed for RNG processes
set.seed(snakemake@params$seed)

# Load performance summaries from CRISPR benchmarking pipeline -------------------------------------

# samples to load
samples <- c("MainPredictors", "PublishedPredictors", "BaselinePredsDHS", "BaselinePredsABC")

# performance summary input files
infiles <- c(snakemake@input$main_perf_summary, snakemake@input$published_perf_summary,
             snakemake@input$baseline_dhs_perf_summary, snakemake@input$baseline_abc_perf_summary)
names(infiles) <- samples

# load all performance summaries
perf <- lapply(infiles, FUN = read_tsv, show_col_types = FALSE)

# Compute performance for enhancer activity predictors ---------------------------------------------

if (snakemake@params$compute_enhAct_performance == TRUE) {
  
  # load merged data and pred config file from enhancer activity benchmarking
  merged_enhAct <- fread(snakemake@input$merged_enhAct)
  pred_config_enhAct <- importPredConfig(snakemake@input$pred_config_enhAct)
  
  # only retain enhancer activity predictors
  pred_config_enhAct <- filter(pred_config_enhAct, pred_id == "EnhActAssayOnly")
  
  # process merged data for benchmarking analyses, including filtering for ValidConnection == TRUE
  merged_enhAct <- processMergedData(merged_enhAct, pred_config = pred_config_enhAct,
                                     filter_valid_connections = TRUE,
                                     include_missing_predictions = TRUE)
  
  # compute performance for all enhancer assay predictors
  perf_enhAct <- makePRSummaryTableBS(merged_enhAct, pred_config = pred_config_enhAct,
                                      pos_col = "Regulated", min_sensitivity = 0.7, R = 1000,
                                      conf = 0.95, ncpus = snakemake@threads)
  
  # add cell type and long predictor name to performance summary
  perf_enhAct <- perf_enhAct %>% 
    as_tibble() %>% 
    left_join(select(pred_config_enhAct, pred_uid, pred_name_long), by = "pred_uid") %>% 
    mutate(cell_type = "K562", .before = 1) %>% 
    relocate(pred_name_long, .after = "pred_uid")
  
  # add to performance summary list
  perf <- c(perf, list(EnhancerAssays = perf_enhAct))
  
}


# Create table with performance for all distal regulation predictors -------------------------------

# function to set performance at threshold columns to precision at minimum sensitivity if no threshold
# was provided
set_prec_threshold <- function(perf) {
  if (!"Threshold" %in% colnames(perf)) {
    perf <- perf %>% 
      mutate(Threshold = round(ThresholdMinSens, digits = 4), PrecThresh = PrecMinSens,
             PrecThresh_lowerCi = PrecMinSens_lowerCi,
             PrecThresh_upperCi = PrecMinSens_upperCi,
             RecallThresh = RecallMinSens, RecallThresh_lowerCi = RecallMinSens_lowerCi,
             RecallThresh_upperCi = RecallMinSens_upperCi)
  }
  return(perf)
}

# create perform at threshold columns for each benchmark if needed
perf <- lapply(perf, set_prec_threshold)

# combine performance summaries into one table
perf <- bind_rows(perf, .id = "sample")

# only pick one entry per predictor
perf <- perf %>% 
  group_by(pred_name_long) %>% 
  slice_head()

# add predictor category for easier to understand output
perf <- perf %>% 
  mutate(category = case_when(
    grepl("Baseline", sample) | grepl("baseline", pred_uid) ~ "Baseline predictor",
    sample == "EnhancerAssays" ~ "Enhancer assay ABC",
    TRUE ~ "Predictive model"
  )) %>% 
  mutate(category = factor(category, levels = c("Predictive model", "Baseline predictor",
                                                "Enhancer assay ABC")))

# reformat for supplementary table
output <- perf %>% 
  select(Predictor = pred_name_long, Category = category, AUPRC, `AUPRC lower CI` = AUPRC_lowerCi,
         `AUPRC upper CI` = AUPRC_upperCi, Threshold, `Precision at threshold` = PrecThresh,
         `Precision at threshold lower CI` = PrecThresh_lowerCi, 
         `Precision at threshold upper CI` = PrecThresh_upperCi) %>% 
  mutate(Threshold = abs(Threshold)) %>%  # convert to absolute value for inverse predictors
  distinct() %>% 
  arrange(Category, desc(AUPRC))
  
# write to output file
write_csv(output, file = snakemake@output[[1]])
