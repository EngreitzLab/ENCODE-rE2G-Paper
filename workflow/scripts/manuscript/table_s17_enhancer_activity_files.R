## Create supplementary table with all ENCODE chromatin files used for enhancer activity ABC models

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
})

# load metadata table containing accession ids for used bigWig files
meta <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# filter for K562 metadata and reformat table for use as supplementary table
output <- meta %>% 
  filter(`Biosample term name` == "K562") %>% 
  select(`Experiment accession`, `File accession`, Assay, `Experiment target`, `Output type`) %>% 
  arrange(Assay)

# write to output file
write_csv(output, file = snakemake@output[[1]], na = "")
