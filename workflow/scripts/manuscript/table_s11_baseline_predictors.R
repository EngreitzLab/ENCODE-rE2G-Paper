## Create supplementary table listing input files and output synapse ids for baseline predictors

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
})

# load table containing input files for all baseline predictors
input_files <- read_tsv(snakemake@input$input_files, show_col_types = FALSE)

# load table containing ENCODE portal accessions for all ABC candidate elements
abc_elements_meta <- read_tsv(snakemake@input$abc_elements_meta, show_col_types = FALSE, skip = 1)

# load table containing ENCODE portal accessions for all ABC and ENCODE-rE2G full predictions
full_pred_meta <- read_tsv(snakemake@input$full_pred_meta, show_col_types = FALSE, skip = 1)

# load synapse manifest tables with synapse IDs for baseline predictor files
dhs_synapse <- read_tsv(snakemake@input$dhs_synapse_manifest, show_col_types = FALSE)
abc_synapse <- read_tsv(snakemake@input$abc_synapse_manifest, show_col_types = FALSE)

# Reformat input files -----------------------------------------------------------------------------

# extract accession of ABC elements for every baseline predictors sample
abc_elements <- abc_elements_meta %>% 
  select(sample = Aliases, ABC_peaks = Accession) %>% 
  mutate(sample = sub("^.+abc_predictions_(.+)_DNaseOnly_candidate_elements", "\\1", sample))

# extract accession of full ABC predictions for every baseline predictors sample
abc_predictions <- full_pred_meta %>% 
  filter(grepl("abc_predictions", Aliases)) %>% 
  select(sample = Aliases, ABC_preds = Accession) %>% 
  mutate(sample = sub("^.+abc_predictions_(.+)_DNaseOnly_full_predictions", "\\1", sample))

# replace ABC peaks and predictions columns containing local paths to files with ENCODE accessions
output <- input_files %>% 
  select(-c(ABC_peaks, ABC_preds)) %>% 
  left_join(abc_elements, by = "sample") %>% 
  left_join(abc_predictions, by = "sample")

# replace 'Null' for cases without bam files with NA
output <- output %>%
  mutate(DHS_bam = replace(DHS_bam, list = DHS_bam == "Null", values = NA_character_),
         H3K27ac_bam = replace(H3K27ac_bam, list = H3K27ac_bam == "Null", values = NA_character_))

# Add synapse IDs for baseline predictor files -----------------------------------------------------

# get ID for synapse directory containing all baseline predictor files for each sample
dhs_synapse <- dhs_synapse %>% 
  mutate(sample = basename(dirname(path))) %>% 
  select(sample, dhs_synapse_id = parent) %>% 
  distinct()

abc_synapse <- abc_synapse %>% 
  mutate(sample = basename(dirname(path))) %>%
  select(sample, abc_synapse_id = parent) %>% 
  distinct()

# add synapse ids for DHS and ABC based baseline predictors to output table
output <- output %>% 
  left_join(dhs_synapse, by = "sample") %>% 
  left_join(abc_synapse, by = "sample")

# create create synapse urls for each sample
output <- output %>% 
  mutate(
    dhs_synapse_url = if_else(!is.na(dhs_synapse_id),
                              true = paste0("https://www.synapse.org/#!Synapse:", dhs_synapse_id),
                              false = NA_character_),
    abc_synapse_url = if_else(!is.na(abc_synapse_id),
                              true = paste0("https://www.synapse.org/#!Synapse:", abc_synapse_id),
                              false = NA_character_))

# Reformat for output ------------------------------------------------------------------------------

# remove ABC peaks file from samples where no ABC baseline predictors were computed
output <- output %>% 
  mutate(ABC_peaks = replace(ABC_peaks, list = is.na(abc_synapse_id), values = NA_character_))

# reformat columns for output and sort
output <- output %>% 
  select(Sample = sample, `Cell type` = cell_type, `DHS peaks` = DHS_peaks, `ABC peaks` = ABC_peaks,
         `ABC predictions` = ABC_preds, `DHS bam file` = DHS_bam, `H3K27ac bam file` = H3K27ac_bam,
         `DHS synapse id` = dhs_synapse_id, `DHS synapse url` = dhs_synapse_url,
         `ABC synapse id` = abc_synapse_id, `ABC synapse url` = abc_synapse_url) %>% 
  arrange(`H3K27ac bam file`, `DHS bam file`, Sample)

# save table to output file
write_csv(output, file = snakemake@output[[1]], na = "")
