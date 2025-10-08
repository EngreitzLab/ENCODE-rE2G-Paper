## Create yaml file containing paths to all ENCODE-rE2G prediction files

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
})

# load ENCODE-rE2G metadata table
e2g_meta_file <- "/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv"
e2g_meta <- read_tsv(e2g_meta_file, show_col_types = FALSE)

## Make list of thresholded prediction files -------------------------------------------------------

# get local paths to all thresholded prediction files
thresholded_preds <- e2g_meta %>% 
  select(Sample, thresholded_predictions_path) %>% 
  deframe() %>% 
  as.list()

# create list with all thresholded ENCODE-rE2G predictions
encode_re2g_predictions <- list(encode_re2g_predictions = list(thresholded = thresholded_preds))

# save list to .yaml file to use as config file in snakemake workflow
write_yaml(encode_re2g_predictions, here("config/encode_re2g_predictions.yml"))

## Make list of corresponding enhancer list files --------------------------------------------------

# directory containing ABC enhancer list files
enh_dir <- "/oak/stanford/groups/engreitz/Users/atan5133/encode_dataset_processing/results/dhs_only"

# infer paths of enhancer list files
enh_files <- e2g_meta %>% 
  mutate(sample_in_path = sub("encode_e2g_predictions_(.+)_thresholded.+", "\\1",
                              basename(thresholded_predictions_path))) %>% 
  mutate(enhancer_list_path = file.path(enh_dir, sample_in_path,
                                        "Neighborhoods/EnhancerList.txt")) %>% 
  select(Sample, enhancer_list_path) %>% 
  deframe() %>% 
  as.list()

# create list with all enhancer list files
encode_re2g_enhancer_lists <- list(encode_re2g_enhancer_lists = enh_files)

# save list to .yaml file to use as config file in snakemake workflow
write_yaml(encode_re2g_enhancer_lists, here("config/encode_re2g_enhancer_lists.yml"))
