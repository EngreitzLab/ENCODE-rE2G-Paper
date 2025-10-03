# Load required packages
suppressPackageStartupMessages({
	library(ggplot2)
	library(data.table)
	library(dplyr)
	library(tidyr)
})



make_gene_annotation_file <- function(abc_gene_file, promoter_class_file, corr_file, out_dir) {
	promoter_class <- fread(promoter_class_file) %>%
		dplyr::select(TargetGene, ubiquitousExpressedGene = is_ubiquitous_uniform, P2PromoterClass)
	gene_corr <- fread(corr_file) %>%
		dplyr::select(TargetGene, averageCorrWeighted) %>%
		left_join(promoter_class, by = "TargetGene") %>%
		mutate(CorrTop25 = (averageCorrWeighted >= quantile(averageCorrWeighted, 0.75, na.rm = TRUE)))
	abc_genes <- fread(abc_gene_file, col.names = c("chr", "start", "end", "TargetGene", "score", "strand", "Ensembl_ID", "gene_type")) %>%
		left_join(gene_corr, by = "TargetGene")

	fwrite(abc_genes, file.path(out_dir, "merged_gene_annotations.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	return(abc_genes)
}

# calculate odds of variants with constitutive eGene in the ABC partition versus specific eGene in Other partition
calculate_eqtl_odds_one_property <- function(property, eqtl) {
	# get eQTL counts the occurrences of each (partition, eGene property) combination
	counts <- eqtl %>%
		group_by(partition, !!sym(property)) %>%
		summarise(count = n(), .groups = "drop") %>%
		pivot_wider(names_from = !!sym(property), values_from = count, values_fill = 0)
	colnames(counts) <- c("partition", "specific", "constitutive")

	# run fisher's test
	A <- counts$constitutive[counts$partition == "ABC"]    # eQTL in ABC, linked to const gene
	B <- counts$specific[counts$partition == "ABC"]    # eQTL in ABC, linked to spec gene
	C <- counts$constitutive[counts$partition == "Other"]  # eQTL in Other, linked to const gene
	D <- counts$specific[counts$partition == "Other"]  # eQTL in Other, linked to spec gene
	fisher_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE) #, dimnames = list(c("ABC_partition", "Other_partition"), c("specific_eGene", "constitutive_eGene")))
	fisher_test_result <- fisher.test(fisher_table)
	#print(fisher_test_result)

	res_df <- data.frame(eGene_class = property,
		odds_ratio = fisher_test_result$estimate,
		p_value = fisher_test_result$p.value,
		CI95_low = fisher_test_result$conf.int[1],
  		CI95_high = fisher_test_result$conf.int[2],
		counts_ABC_constitutive = A,
		counts_ABC_specific = B,
		counts_Other_constitutive = C,
		counts_Other_specific = D)
}


# calculate the likelihood of eQTLs in the ABC region being linked to eGenes of one class versus another
run_eqtl_analysis <- function(eqtl_file, gene_annot_file, gene_classes, out_dir) {
	gene_annot <- fread(gene_annot_file, sep = "\t")
	eqtl <- fread(eqtl_file, sep = "\t", col.names = c("chr", "start", "end", "TargetGene", "partition")) %>%
		mutate(partition = ifelse(partition == "ABC", "ABC", "Other")) %>%
		left_join(gene_annot, by = "TargetGene")

	# calculate
	res <- lapply(names(gene_classes), calculate_eqtl_odds_one_property, eqtl) %>% 
		rbindlist() %>% as_tibble() %>%
		mutate(description = gene_classes[eGene_class],
			description = factor(description, levels = rev(gene_classes), ordered = TRUE))
	fwrite(res, file.path(out_dir, "eqtl_gene_class_by_ABC_overlap_odds.tsv"), sep = "\t", col.names = TRUE)

	# plot
	p <- ggplot(res, aes(x = odds_ratio, y = description)) +
		geom_vline(aes(xintercept = 1), linetype = "dashed", color = "#c5cad7") +
		geom_errorbarh(aes(xmin = CI95_low, xmax = CI95_high), height = 0.25) +
		geom_point(shape = 16, size = 4) +
		xlim(c(0.5, 1.05)) +
		labs(x = "Odds of variant in ABC enhancer\nhaving eGene in given class",
			y = "Class of eGene") +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"))
	ggsave(file.path(out_dir, "eqtl_gene_class_by_ABC_overlap.pdf"), p, width = 4, height = 3)
	
}

# calculate odds of variants with constitutive eGene in the ABC partition versus specific eGene in Other partition
calculate_crispr_odds_one_property <- function(property, crispr) {
	# get counts the occurrences of each (Regulated, gene property) combination
	counts <- crispr %>%
		mutate(!!sym(property) := replace_na(!!sym(property), FALSE)) %>%
		group_by(Regulated, !!sym(property)) %>%
		summarise(count = n(), .groups = "drop") %>%
		pivot_wider(names_from = !!sym(property), values_from = count, values_fill = 0)
	colnames(counts) <- c("regulated", "specific", "constitutive")

	# run fisher's test
	A <- counts$constitutive[counts$regulated]    # regulated EG pair, linked to const gene
	B <- counts$specific[counts$regulated]    # regulated EG pair, linked to spec gene
	C <- counts$constitutive[!counts$regulated]  # not-regulated EG pair, linked to const gene
	D <- counts$specific[!counts$regulated]  # not-regulated EG pair, linked to spec gene
	fisher_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE) 
	fisher_test_result <- fisher.test(fisher_table)
	print(fisher_test_result)

	res_df <- data.frame(target_gene_class = property,
		odds_ratio = fisher_test_result$estimate,
		p_value = fisher_test_result$p.value,
		CI95_low = fisher_test_result$conf.int[1],
  		CI95_high = fisher_test_result$conf.int[2],
		counts_regulated_constitutive = A,
		counts_regulated_specific = B,
		counts_not_regulated_constitutive = C,
		counts_not_regulated_specific = D)
}


run_crispr_analysis <- function(crispr_file, gene_annot_file, gene_classes, out_dir) {
	gene_annot <- fread(gene_annot_file, sep = "\t")
	crispr <- fread(crispr_file, sep = "\t") %>%
		select(name, measuredGeneSymbol, Regulated) %>% distinct() %>%
		mutate(TargetGene = measuredGeneSymbol) %>%
		left_join(gene_annot, by = "TargetGene")

	# calculate
	res <- lapply(names(gene_classes), calculate_crispr_odds_one_property, crispr) %>% 
		rbindlist() %>% as_tibble() %>%
		mutate(description = gene_classes[target_gene_class],
			description = factor(description, levels = rev(gene_classes), ordered = TRUE))
	fwrite(res, file.path(out_dir, "crispri_target_gene_class_by_regulated_odds.tsv"), sep = "\t", col.names = TRUE)

	# plot
	p <- ggplot(res, aes(x = odds_ratio, y = description)) +
		geom_vline(aes(xintercept = 1), linetype = "dashed", color = "#c5cad7") +
		geom_errorbarh(aes(xmin = CI95_low, xmax = CI95_high), height = 0.25) +
		geom_point(shape = 16, size = 4) +
		xlim(c(0, 1.05)) +
		labs(x = "Odds of CRISPRi-validated enhancer\nhaving target gene in given class",
			y = "Class of target gene") +
		theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"), axis.title = element_text(size = 8),
			axis.ticks = element_line(color = "#000000"))
	ggsave(file.path(out_dir, "crispri_target_gene_class_by_regulated.pdf"), p, width = 4, height = 3)
	
}

format_crispr_data_for_models <- function(crispr_file, gene_annot, out_dir) {
	gene_annot <- fread(gene_annot_file, sep = "\t") %>%
		select(TargetGene, ubiquitousExpressedGene, P2PromoterClass, averageCorrWeighted)

	crispr <- fread(crispr_file, sep = "\t") %>%
		filter(pred_uid == "ABCfull.ABC.Score") %>%
		select(chr = chrom, start = chromStart, end = chromEnd, name, TargetGene = measuredGeneSymbol,
			Regulated, ABC.Score = pred_value) %>% 
		distinct() %>%
		left_join(gene_annot, by = "TargetGene") %>%
		mutate(across(c(ubiquitousExpressedGene, P2PromoterClass), ~ ifelse(is.na(.), 0, as.integer(.))),
			averageCorrWeighted = replace_na(averageCorrWeighted, 0))

	fwrite(crispr, file.path(out_dir, "CRISPR_training_data_for_models.tsv.gz"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

#### MAIN

## file paths
promoter_class_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv"
corr_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/external_features/EnhancerWeightedCorr_acrossCellTypes_RefSeqCurated.tsv"
abc_gene_file <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/ENCODE_rE2G/reference/CollapsedGeneBounds.hg38.TSS500bp.bed"
eqtl_file_initial <- "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GTExVariants/GTEx_30tissues_hg38.SUSIE.cs.expressedTPM1.tsv.gz"
crispr_with_abc_file <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz"
#dataset	chrom	chromStart	chromEnd	name	EffectSize	chrTSS	startTSS	endTSS	measuredGeneSymbol	Significant	pValueAdjusted	PowerAtEffectSize25	ValidConnection	ExperimentCellType	Reference	Regulated	PowerAtEffectSize10	PowerAtEffectSize15	PowerAtEffectSize20	PowerAtEffectSize50	pair_uid	merged_uid	merged_start	merged_end	expressed	PredictionCellType	pred_uid	pred_id	pred_col	pred_value	Prediction

out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0129_promoter_class_analysis"
out_ref <- file.path(out_dir, "reference"); dir.create(out_ref, showWarnings = FALSE)
out_res <- file.path(out_dir, "eqtl_crispr_results"); dir.create(out_res, showWarnings = FALSE)

eqtl_file <- file.path(out_ref, "GTEx_variants.any_tissue.PIP0.5.txt.gz") # chr start end gene
gene_annot_file <- file.path(out_ref, "merged_gene_annotations.tsv")

gene_classes <- c(ubiquitousExpressedGene = "Ubiquitously-expressed\nversus not",
	P2PromoterClass = "P2 promoter\nversus not P2 promoter",
	CorrTop25 = "Top 25% enhancer correlation\nversus bottom 75%")

## RUN

#gene_annot <- make_gene_annotation_file(abc_gene_file, promoter_class_file, corr_file, out_ref)
#run_eqtl_analysis(eqtl_file, gene_annot_file, gene_classes, out_res)
#run_crispr_analysis(crispr_with_abc_file, gene_annot_file, gene_classes, out_res)
format_crispr_data_for_models(crispr_with_abc_file, gene_annot, out_ref)