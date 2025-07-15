## Define default samples based on the heuristic to define the 352 biosample-treatment combinations
## in the preprint

library(tidyverse)
library(here)

# load ENCODE-rE2G metadata table
meta_file <- "/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv"
meta <- read_tsv(meta_file, show_col_types = FALSE)

# extract center information from lab and create ordered factor representing priorities
meta <- meta %>% 
  mutate(`DNase Center` = case_when(
    grepl("UW", `DNase Lab`) ~ "UW",
    grepl("Duke", `DNase Lab`) ~ "Duke",
    grepl("UMass", `DNase Lab`) ~ "UMass",
    TRUE ~ NA_character_
  )) %>% 
  mutate(`DNase Center` = factor(`DNase Center`, levels = c("UW", "Duke", "UMass"), ordered = TRUE))

# Pick one default experiment for each unique Biosample based on DNase center priority, no treatment
# if possible, the most recent experiment. If more than 1 left, pick random
default_experiments <- meta %>% 
  mutate(treated = !is.na(`Biosample treatments`)) %>% 
  group_by(`Biosample term name`) %>% 
  mutate(random_pick = sample(row_number())) %>% 
  arrange(`Biosample term name`, `DNase Center`, treated, desc(`DNase Experiment date released`),
          random_pick) %>% 
  slice_head(n = 1) %>% 
  pull(Sample)

# add this information to the metadata table and save for output
