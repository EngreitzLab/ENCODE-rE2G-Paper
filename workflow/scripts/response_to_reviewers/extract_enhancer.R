## Extract unique enhancers from combined K562 CRISPR training dataset, including information
## whether each enhancer is functional (has >0 hits) or not (0 hits). These are then used by the
## ENCODE functional characterization group to overlap with MPRA, STARR-seq etc.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(here)
})

# load combined K562 CRISPR data
crispr_file <- "/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/CRISPR_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
crispr <- read_tsv(crispr_file, show_col_types = FALSE)

# get list of unique elements and whether they have at least on CRISPR linked target gene
elements <- crispr %>% 
  mutate(name = paste0(chrom, ":", chromStart, "-", chromEnd)) %>% 
  group_by(chrom, chromStart, chromEnd, name, dataset) %>% 
  summarize(Regulated = any(Regulated == TRUE))

# save to output file
outfile <- here("analyses/manuscript_figures/revisions/response_to_reviewers/crispr_enhancers.tsv")
write_tsv(elements, file = outfile)
