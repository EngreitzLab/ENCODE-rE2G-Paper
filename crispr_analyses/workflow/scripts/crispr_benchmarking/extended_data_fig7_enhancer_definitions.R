## Plot overlaps of ENCODE-rE2G predicted and CRISPR enhancers with cCRE ELS elements and chromHMM
## enhancer states

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# Load all data ------------------------------------------------------------------------------------

# load overlaps
enh_ccre_ovl <- read_tsv(snakemake@input$e2g_ccres, col_names = FALSE, show_col_types = FALSE)
crispr_ccre_ovl <- read_tsv(snakemake@input$crispr_ccres, col_names = FALSE, show_col_types = FALSE)
enh_chmm_ovl <- read_tsv(snakemake@input$e2g_chromhmm, col_names = FALSE, show_col_types = FALSE)
crispr_chmm_ovl <- read_tsv(snakemake@input$crispr_chromhmm, col_names = FALSE,
                            show_col_types = FALSE)

# compute overlaps ---------------------------------------------------------------------------------

# function to compute
get_overlaps <- function(ovl, feat_class_col, filt_enh_classes = NULL, filt_feat_class = NULL) {
  
  # reformat to simple table
  ovl_simple <- ovl %>% 
    select(chr = 1, start = 2, end = 3, class = 4, feat_class = all_of(feat_class_col)) %>% 
    distinct()
  
  # filter based on class if specified
  if (!is.null(filt_enh_classes)) {
    ovl_simple <- filter(ovl_simple, class %in% filt_enh_classes)
  }
  
  # filter for feature classes if specified
  if (!is.null(filt_feat_class)) {
    ovl_simple <- ovl_simple %>% 
      mutate(feat_class = if_else(feat_class %in% filt_feat_class, true = feat_class, false = "."))
  }
  
  return(ovl_simple)
  
}

# function to compute percentage overlap
compute_perc_overlap <- function(ovl, feat_class_hierarchy) {
  
  # if there are multiple overlaps, pick only one
  ovl <- ovl %>% 
    mutate(feat_class = factor(feat_class, levels = feat_class_hierarchy)) %>% 
    group_by(chr, start, end) %>% 
    slice_max(order_by = feat_class, n = 1)
  
  # count the percentage of enhancers overlapping each cCRE class
  perc_ovl_classes <- ovl %>% 
    group_by(feat_class) %>% 
    summarize(sum_enhancers = n()) %>% 
    mutate(perc_enhancers = sum_enhancers / sum(sum_enhancers))
  
  return(perc_ovl_classes)
  
}

# process overlaps for cCREs
enh_ccre_ovl_proc <- get_overlaps(enh_ccre_ovl, feat_class_col = 14,
                                  filt_enh_classes = c("genic", "intergenic"),
                                  filt_feat_class = c("pELS", "dELS"))

crispr_ccre_ovl_proc <- get_overlaps(crispr_ccre_ovl, feat_class_col = 14,
                                     filt_feat_class = c("pELS", "dELS"))

# chromHMM enhancer states
chromHMM_enh <- c("EnhWk", "EnhG2", "EnhG1", "EnhBiv", "EnhA2", "EnhA1")

# process overlaps for chromHMM
enh_chmm_ovl_proc <- get_overlaps(enh_chmm_ovl, feat_class_col = 8,
                                  filt_enh_classes = c("genic", "intergenic"),
                                  filt_feat_class = chromHMM_enh)

crispr_chmm_ovl_proc <- get_overlaps(crispr_chmm_ovl, feat_class_col = 8,
                                     filt_feat_class = chromHMM_enh)

# compute percentage of overlap for cCREs
enh_ccre_perc <- compute_perc_overlap(enh_ccre_ovl_proc, feat_class_hierarchy = c(".", "pELS", "dELS"))
crispr_ccre_perc <- compute_perc_overlap(crispr_ccre_ovl_proc, feat_class_hierarchy = c(".", "pELS", "dELS"))

# compute percentage of overlap for chromHMM states
enh_chmm_perc <- compute_perc_overlap(enh_chmm_ovl_proc, feat_class_hierarchy = c(".", chromHMM_enh))
crispr_chmm_perc <- compute_perc_overlap(crispr_chmm_ovl_proc, feat_class_hierarchy = c(".", chromHMM_enh))

# combine into tables
perc_ovl_ccre <- bind_rows(`ENCODE-rE2G` = enh_ccre_perc, CRISPR = crispr_ccre_perc, .id = "dataset")
perc_ovl_chmm <- bind_rows(`ENCODE-rE2G` = enh_chmm_perc, CRISPR = crispr_chmm_perc, .id = "dataset")

# remove non-overlaps for plot
perc_ovl_ccre <- filter(perc_ovl_ccre, feat_class != ".")
perc_ovl_chmm <- filter(perc_ovl_chmm, feat_class != ".")

# set all factor levels
perc_ovl_ccre <- mutate(perc_ovl_ccre, feat_class = factor(feat_class, levels = c("pELS", "dELS")))
perc_ovl_chmm <- mutate(perc_ovl_chmm, feat_class = factor(feat_class, levels = chromHMM_enh))

# calculate the summed up percentages
perc_ovl_ccre_total <- perc_ovl_ccre %>%
  group_by(dataset) %>%
  summarize(total_perc = sum(perc_enhancers)) %>% 
  mutate(label = paste0(round(total_perc * 100, digits = 2), "%"))

perc_ovl_chmm_total <- perc_ovl_chmm %>%
  group_by(dataset) %>%
  summarize(total_perc = sum(perc_enhancers)) %>% 
  mutate(label = paste0(round(total_perc * 100, digits = 2), "%"))

# make plot for cCREs
p1 <- ggplot(perc_ovl_ccre, aes(x = dataset, y = perc_enhancers)) +
  geom_bar(stat = "identity", aes(fill = feat_class)) +
  geom_text(data = perc_ovl_ccre_total, aes(x = dataset, y = total_perc, label = label),
            vjust = -0.5) +
  labs(y = "Percentage of\ninferred enhancers", title = "ENCODE4 cCRES (ELS)",
       fill = "cCRE ELS\nclass") +
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank())

p2 <- ggplot(perc_ovl_chmm, aes(x = dataset, y = perc_enhancers)) +
  geom_bar(stat = "identity", aes(fill = feat_class)) +
  geom_text(data = perc_ovl_chmm_total, aes(x = dataset, y = total_perc, label = label),
            vjust = -0.5) +
  labs(y = "Percentage of\ninferred enhancers", title = "ChromHMM enhancer states",
       fill = "ChromHMM\nstate") +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank())

# CRISPR data by dataset ---------------------------------------------------------------------------

# split CRSIPR overlaps by dataset
crispr_ccre_ovl <- split(crispr_ccre_ovl, f = crispr_ccre_ovl$X4)
crispr_chmm_ovl <- split(crispr_chmm_ovl, f = crispr_chmm_ovl$X4)

# process overlaps
crispr_ccre_ovl <- lapply(crispr_ccre_ovl, FUN = get_overlaps, feat_class_col = 14,
                          filt_feat_class = c("pELS", "dELS"))
crispr_chmm_ovl <- lapply(crispr_chmm_ovl, FUN = get_overlaps, feat_class_col = 8,
                          filt_feat_class = chromHMM_enh)

# compute percentage of overlap
crispr_ccre_perc <- lapply(crispr_ccre_ovl, FUN = compute_perc_overlap,
                           feat_class_hierarchy = c(".", "pELS", "dELS"))
crispr_chmm_perc <- lapply(crispr_chmm_ovl, FUN = compute_perc_overlap,
                           feat_class_hierarchy = c(".", chromHMM_enh))

# combine into tables
crispr_ccre_perc <- bind_rows(crispr_ccre_perc, .id = "dataset")
crispr_chmm_perc <- bind_rows(crispr_chmm_perc, .id = "dataset")

# remove non-overlaps for plot
crispr_ccre_perc <- filter(crispr_ccre_perc, feat_class != ".")
crispr_chmm_perc <- filter(crispr_chmm_perc, feat_class != ".")

# set all factor levels
crispr_ccre_perc <- mutate(crispr_ccre_perc, feat_class = factor(feat_class, levels = c("pELS", "dELS")))
crispr_chmm_perc <- mutate(crispr_chmm_perc, feat_class = factor(feat_class, levels = chromHMM_enh))

# set better CRISPR dataset labels
crispr_ccre_perc <- crispr_ccre_perc %>% 
  mutate(dataset = case_when(
    dataset == "Nasser2021" ~ "Nasser et al., 2021",
    dataset == "Gasperini2019" ~ "Gasperini et al., 2019",
    dataset == "Schraivogel2020" ~ "Schraivogel et al., 2020"
  )) %>% 
  mutate(dataset = fct_relevel(dataset, "Nasser et al., 2021", "Gasperini et al., 2019"))

crispr_chmm_perc <- crispr_chmm_perc %>% 
  mutate(dataset = case_when(
    dataset == "Nasser2021" ~ "Nasser et al., 2021",
    dataset == "Gasperini2019" ~ "Gasperini et al., 2019",
    dataset == "Schraivogel2020" ~ "Schraivogel et al., 2020"
  )) %>% 
  mutate(dataset = fct_relevel(dataset, "Nasser et al., 2021", "Gasperini et al., 2019"))

# calculate the summed up percentages
crispr_ccre_perc_total <- crispr_ccre_perc %>%
  group_by(dataset) %>%
  summarize(total_perc = sum(perc_enhancers)) %>% 
  mutate(label = paste0(round(total_perc * 100, digits = 2), "%"))

crispr_chmm_perc_total <- crispr_chmm_perc %>%
  group_by(dataset) %>%
  summarize(total_perc = sum(perc_enhancers)) %>% 
  mutate(label = paste0(round(total_perc * 100, digits = 2), "%"))

# make datset stratified plots
p3 <- ggplot(crispr_ccre_perc, aes(x = dataset, y = perc_enhancers)) +
  geom_bar(stat = "identity", aes(fill = feat_class)) +
  geom_text(data = crispr_ccre_perc_total, aes(x = dataset, y = total_perc, label = label),
            vjust = -0.5) +
  labs(y = "Percentage of\nCRISPR enhancers", title = "ENCODE4 cCRES (ELS)",
       fill = "cCRE ELS\nclass") +
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank())

p4 <- ggplot(crispr_chmm_perc, aes(x = dataset, y = perc_enhancers)) +
  geom_bar(stat = "identity", aes(fill = feat_class)) +
  geom_text(data = crispr_chmm_perc_total, aes(x = dataset, y = total_perc, label = label),
            vjust = -0.5) +
  labs(y = "Percentage of\nCRISPR enhancers", title = "ChromHMM enhancer states",
       fill = "ChromHMM\nstate") +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank())

# arrange plots into figure
fig <- plot_grid(p1, p2, p3, p4, labels = c("a", "b", "c", "d"), nrow = 2)

# save figure to .pdf file
ggsave(fig, filename = snakemake@output[[1]], height = 7, width = 8)
