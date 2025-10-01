## Process ENCODE DNase-seq metadata to select bam files for each released human DNase-seq
## experiment as input for ENCODE-rE2G

library(tidyverse)

# load DNase-seq metadata obtained from the ENCODE portal
bam   <- read_tsv("dnase_metadata_bam_230927.tsv.gz",   show_col_types = FALSE)
fastq <- read_tsv("dnase_metadata_fastq_230927.tsv.gz", show_col_types = FALSE)

# only retain bam files with reads mapped to GRCh38 genome build
bam <- filter(bam, `File assembly` == "GRCh38")

# get run type (paired- or single-ended) and lab from fastq metadata
run_types <- distinct(select(fastq, `File accession`, `Run type`, Lab))

# filter for bam files produced by the ENCODE4 DNase-seq pipeline only and remove archived files
bam <- bam %>% 
  filter(`File Status` == "released") %>% 
  filter(grepl("^ENCODE4.+", `File analysis title`))

# extract relevant DNase-seq columns from bam metadata 
meta <- bam %>% 
  select(`Donor(s)`, starts_with("Biosample"), `Technical replicate(s)`, `Experiment accession`,
         `Experiment date released`, `Output type`, `File accession`, `Derived from`)

# split "Derived from" files into list and convert to long format to be able to merge with run types
meta <- meta %>% 
  mutate(`Derived from` = strsplit(`Derived from`, split = ", ")) %>% 
  unnest(cols = `Derived from`) %>% 
  mutate(`Derived from` = basename(`Derived from`))

# merge with run types, which will add information for all unfiltered alignments derived from fastq
meta_unfilt <- inner_join(meta, run_types, by = c("Derived from" = "File accession"))

# many bam files are created from combining multiple fastq files, sometimes with different run types
# and from different labs, which are combined into single strings
meta_unfilt <- meta_unfilt %>% 
  group_by(across(-c(`Derived from`, `Run type`, Lab))) %>% 
  summarize(`Run type` = paste(sort(unique(`Run type`)), collapse = " & "),
            Lab = paste(sort(unique(Lab)), collapse = " & "),
            .groups = "drop")

# extract run types of all unfiltered alignments from just created table
unfilt_run_types <- select(meta_unfilt, `File accession`, `Run type`, Lab)

# add run type for all filtered alignments based on run types of unfiltered aligments they are
meta_filt <- meta %>% 
  inner_join(unfilt_run_types, by = c("Derived from" = "File accession")) %>% 
  group_by(across(-c(`Derived from`, `Run type`, Lab))) %>% 
  summarize(`Run type` = paste(sort(unique(`Run type`)), collapse = "/"),
            Lab = paste(sort(unique(Lab)), collapse = "; "),
            .groups = "drop")

# combine metadata for unfiltered and filtered alignments into one table
meta <- arrange(bind_rows(meta_unfilt, meta_filt), `Experiment accession`)

# for single-ended reads only retain unfiltered alignments and for paired-ended retain only filtered
# alignments. cases with mixed run types are considered single-ended
meta <- meta %>% 
  filter((`Run type` %in% c("single-ended", "paired-ended & single-ended") &
            `Output type` == "unfiltered alignments") |
           (`Run type` == "paired-ended" & `Output type` == "alignments"))

# multiple bam files per experiment should now represent replicates for that experiment and are
# now combined into one comma separated string per experiment. Donor(s), Output type and Run type
# also get combined into comma separated strings
meta <- meta %>% 
  group_by(across(-c(`Donor(s)`, `Technical replicate(s)`, `Output type`, `File accession`,
                     `Run type`))) %>% 
  summarize(`Donor(s)` = paste(unique(`Donor(s)`), collapse = ", "),
            `File accession` = paste(`File accession`, collapse = ", "),
            `Output type` = paste(`Output type`, collapse = ", "),
            `Run type` = paste(`Run type`, collapse = ", "),
            .groups = "drop")

# write processed DNase-seq metadata to file
write_tsv(meta, file = "encode_re2g_dnase_input_metadata.tsv")
