## Combine estimated indirect effect ratios into one table for supplementary material

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
})

# all indirect stats files
stat_files <- snakemake@input[names(snakemake@input) != ""]

# load all indidrect effects stats files combine into one table
stats <- stat_files %>% 
  lapply(FUN = read_csv, show_col_types = FALSE) %>% 
  bind_rows(.id = "dataset")

# save to output file
write_csv(stats, file = snakemake@output$stats_all_dist)

# create alternative version of the table with only pairs within 1Mb from TSS
stats_1mb <- stats %>% 
  filter(max_dist == "1mb") %>% 
  select(-max_dist)

# save to output file
write_csv(stats, file = snakemake@output$stats_1mb_dist)
