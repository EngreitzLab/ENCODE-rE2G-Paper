### --- DEFINE CHROMATIN CATEGORIES ---
# get thresholds based on genome-wide elements
get_category_thresholds <- function(enh, quantiles) {
    enh_ctcf <- enh %>% filter(CTCF_peak_overlap == 1) %>% 
        select(cell_type, CTCF.H3K27ac.ratio.CTCF_peak = CTCF.H3K27ac.ratio) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))
    
    enh_h3k27ac <- enh %>% filter(H3K27ac_peak_overlap == 1) %>% 
        select(cell_type, H3K27ac.RPM.H3K27ac_peak = H3K27ac.RPM) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_h3k27me3 <- enh %>% filter(H3K27me3_peak_overlap == 1) %>% 
        select(cell_type, H3K27me3.RPM.expandedRegion.H3K27me3_peak = H3K27me3.RPM.expandedRegion) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_other <- enh %>% 
        select(cell_type, CTCF.RPM, H3K27ac.RPM, H3K27ac.RPM.expandedRegion, DHS.RPM) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    res <- rbind(enh_ctcf, enh_h3k27ac, enh_h3k27me3, enh_other)
    
    return(res)
}

# get table of quantile values from genome-wide elements
get_threshold_key <- function(thresholds, feature_col, quantile_this) {
    filt <- thresholds %>% filter(feature == feature_col, quantile == quantile_this)
    key <- setNames(filt$value, filt$cell_type)
    return(key)
}

# categorize elements
categorize_elements <- function(enh, thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) {
    ### CATEGORIZATION LOGIC ###
    # if element overlaps H3K27ac peak:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion < 90% --> H3K27ac
    # else:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion > 50% --> H3K27ac
        # if element overlaps CTCF peak --> CTCF element
        # if element overlaps H3K27me3 peak --> H3K27me3 element
        # else: No H3K27ac
        

    key_high <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_high)
    key_low <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_low)

    enh <- enh %>%
         mutate(element_category = case_when(
                H3K27ac_peak_overlap == 1 & H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac_peak_overlap == 1 ~ "H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_low[cell_type] ~ "H3K27ac",
                CTCF_peak_overlap == 1 ~ "CTCF element",
                H3K27me3_peak_overlap == 1 ~ "H3K27me3 element",
                TRUE ~ "No H3K27ac"))
    return(enh)
}

# summarize properties of genome-wide element-gene pairs
annotate_genomewide_pairs <- function(enh, e2g_files, cell_types, remove_promoters, distance_threshold) {
    res_list <- vector("list", length(cell_types))
    res_list_genes <- vector("list", length(cell_types))

    for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        print(ct)

        enh_ct <- filter(enh, cell_type == ct)
        pairs_file <- e2g_files[ct]; print(pairs_file)

        pairs <- fread(pairs_file, sep = "\t") %>% 
            select(chr, start, end, class, TargetGene, ubiquitousExpressedGene, distanceToTSS, CellType) %>% 
            filter(distanceToTSS < distance_threshold) %>%
            mutate(distance_category = case_when(distanceToTSS < 10e3 ~  "0-10 kb",
                                             distanceToTSS < 100e3 ~ "10-100 kb",
                                             distanceToTSS < 250e3 ~ "100-250 kb",
                                             distanceToTSS < 1000e3 ~ "250 kb-1 Mb",
                                             distanceToTSS < 2000e3 ~ "1 Mb-2 Mb",
                                             TRUE ~ ">2 Mb"),
                    ubiq_category = ifelse(ubiquitousExpressedGene %in% c("True", TRUE), "Ubiq. expr. gene", "Other gene"))

        if (remove_promoters) {
            pairs <- filter(pairs, class != "promoter")
        }

        res_list[[i]] <- left_join(pairs, enh_ct, by = c("chr", "start", "end")) %>% 
            group_by(cell_type, element_category, distance_category, ubiq_category) %>% 
            summarize(n_pairs = n())

        res_list_genes[[i]] <- pairs %>% select(cell_type = CellType, TargetGene, ubiq_category) %>% distinct() %>% 
            group_by(cell_type, ubiq_category) %>%
            summarize(n_e2g_genes = n())
    }

    res <- rbindlist(res_list) %>% as.data.frame()
    res_genes <- rbindlist(res_list_genes) %>% as.data.frame()

    return(list(res, res_genes))
}
