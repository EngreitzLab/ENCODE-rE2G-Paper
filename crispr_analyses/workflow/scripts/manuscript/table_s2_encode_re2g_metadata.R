## Create supplementary table with metadata including ENCODE portal accessions for all ENCODE-rE2G
## and ENCODE-rE2G_Extended predictions

suppressPackageStartupMessages({
  library(tidyverse)
})

# Create metadata table for ENCODE-rE2G predictions ------------------------------------------------

# load ENCODE-rE2G metadata table
e2g_meta <- read_tsv(snakemake@input$e2g_meta, show_col_types = FALSE)

# load ENCODE-rE2G metadata from ENCODE portal
e2g_portal <- read_tsv(snakemake@input$e2g_portal_meta, show_col_types = FALSE, skip = 1)

# get annotation object accession number and input DNA-seq experiment for ENCODE-rE2G predictions
e2g_accessions <- e2g_portal %>% 
  filter(!grepl("extended", Aliases)) %>% 
  select(`ENCODE-rE2G accession` = Accession,
         `DNase Experiment accession` = `Experimental input`) %>% 
  mutate(`DNase Experiment accession` = basename(`DNase Experiment accession`))

# add accessions to ENCODE-rE2G metadata table
e2g_meta <- e2g_meta %>% 
  left_join(e2g_accessions, by = "DNase Experiment accession") %>% 
  relocate(`ENCODE-rE2G accession`, .before = 1)

# select columns to include in supplementary table
e2g_meta <- e2g_meta %>% 
  select(-c(`Default Sample`, ends_with("_path"), ends_with("_mitra"), ends_with("_synapse"),
            dnase_bigwig_url))

# Create metadata table for ENCODE-rE2G_Extended predictions ---------------------------------------

# load list of additional ENCODE-rE2G_Extended assays
extended_assays <- read_tsv(snakemake@input$extended_assays, show_col_types = FALSE)

# extract accessions extended model
e2g_ext_accessions <- e2g_portal %>% 
  mutate(Aliases = strsplit(Aliases, ",")) %>% 
  unnest(Aliases) %>% 
  filter(grepl("extended", Aliases)) %>% 
  mutate(`Cell type` = sub(".+predictions_(.+)", "\\1", Aliases)) %>% 
  select(`ENCODE-rE2G accession` = Accession, `Cell type`)

# add input assays
e2g_ext_assays <- e2g_ext_accessions %>% 
  left_join(extended_assays, by = "Cell type") %>% 
  select(-`Cell type`)

# combine with other DNase-seq metadata in common with the basic ENCODE-rE2G model
e2g_ext_meta <- e2g_meta %>% 
  select(-c(`ENCODE-rE2G accession`, Sample)) %>% 
  left_join(e2g_ext_assays, ., by = "DNase Experiment accession")

# correctly order columns for supplementary tables
table_cols <- c(colnames(e2g_meta), setdiff(colnames(e2g_ext_meta), colnames(e2g_meta)))
e2g_ext_meta <- select(e2g_ext_meta, all_of(table_cols))

# Save to output files -----------------------------------------------------------------------------

# create combined metadata table for both ENCODE-rE2G and ENCODE-rE2G_extended
e2g_meta_combined <- bind_rows(Extended = e2g_ext_meta, `DNase-only` = e2g_meta, .id = "Model")

# save metadata tables to files
write_csv(e2g_meta, file = snakemake@output$dnase_table, na = "")
write_csv(e2g_ext_meta, file = snakemake@output$extended_table, na = "")
write_csv(e2g_meta_combined, file = snakemake@output$combined_table, na = "")
