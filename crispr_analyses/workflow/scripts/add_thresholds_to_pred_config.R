
library(tidyverse)
library(here)

# load pred config file
pred_conig_file <- "/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/benchmarking_pred_config.tsv"
pred_config <- read_tsv(pred_conig_file, show_col_types = FALSE)

# performance summary files to include
comparisons <- c("MainPredictors", "PublishedPredictors", "BaselinePredsDHS", "BaselinePredsABC")
infiles <- here("CRISPR_benchmarks/results", comparisons, "performance_summary.txt")
names(infiles) <- comparisons

# load all performance summaries
perf <- lapply(infiles, FUN = read_tsv, show_col_types = FALSE)