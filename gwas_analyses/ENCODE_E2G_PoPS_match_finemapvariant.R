library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--trait", type="character", default = "data/traits.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

trait = opt$trait

#trait = "AFib"

################## Number of credible sets per trait ###################################################

eg_preds = data.frame(fread("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/ENCODE_E2G2024_ALL_EG_predictions.txt"))
gr1 = GRanges(seqnames = eg_preds$chr,
              ranges = IRanges(start=eg_preds$start, end = eg_preds$end),
              gene = eg_preds$TargetGene,
              score = eg_preds$Score)

pops_tabb = data.frame(fread(paste0("/data/deyk/kushal/PolyPerturb/codes/PoPS/supplemental_data/", "PoPS_FullResults.txt.gz")))

ll = list.files("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets")

binary_encode_pops = c()

pops_tabb_trait = pops_tabb[which(pops_tabb$trait == trait), ]
gr1_pops = GRanges(seqnames = paste0("chr", pops_tabb_trait$chromosome),
                ranges = IRanges(start=pmax(pops_tabb_trait$tss-500000, 0),
                                 end = pops_tabb_trait$tss + 500000),
                gene = pops_tabb_trait$gene,
                score = pops_tabb_trait$pops_score)


dff2 = data.frame(fread(paste0("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets",
                                 "/", trait, "/", "variant.list.txt")))
dff2 = dff2[which(dff2$pip > 0.1), ]
uu = unique(dff2$rsid)

bimlist = list()
for(numchr in 1:22){
  xx = data.frame(fread(paste0("/data/deyk/kushal/extras/BIMS_hg38/", "1000G.EUR.QC.", numchr, ".bim")))
  colnames(xx) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
  bimlist[[numchr]] = xx
  cat("Processing bimfile for chr:", numchr, "\n")
}
names(bimlist) = paste0("chr", 1:22)
bimfile_merged = do.call(rbind, bimlist)
bimfile_merged_trait = bimfile_merged[which(bimfile_merged$SNP %in% uu == T), ]


gr2 = GRanges(seqnames = paste0("chr", bimfile_merged_trait$CHR),
              ranges = IRanges(start=bimfile_merged_trait$BP-1,
                               end = bimfile_merged_trait$BP),
              snp_name = bimfile_merged_trait$SNP)

cc=findOverlaps(gr2, gr1, type = "any", select = "all")
xtabb = cbind.data.frame(gr2[queryHits(cc)], gr1[subjectHits(cc)])
colnames(xtabb) = c("var.seqnames", "var.start", "var.end", "var.width", "var.strand",
                    "snp_name", "seqnames", "start", "end", "width", "strand", "gene", "score")

cc2=findOverlaps(gr2, gr1_pops, type = "any", select = "all")
xtabb2 = cbind.data.frame(gr2[queryHits(cc2)], gr1_pops[subjectHits(cc2)])
colnames(xtabb2) = c("var.seqnames", "var.start", "var.end", "var.width", "var.strand",
                    "snp_name", "seqnames", "start", "end", "width", "strand", "gene", "score")


library(dplyr)
top1 <- xtabb %>%
  group_by(snp_name) %>%
  slice_max(order_by = score, n = 2, with_ties = FALSE) %>%
  summarise(top_genes1 = list(unique(gene)))

top2 <- xtabb2 %>%
  group_by(snp_name) %>%
  slice_max(order_by = score, n = 2, with_ties = FALSE) %>%
  summarise(top_genes2 = list(unique(gene)))

# 2. Join by snp_name
top_joined <- inner_join(top1, top2, by = "snp_name")

# 3. Filter for SNPs with at least one overlapping gene in top 2 sets
top_joined_with_overlap <- top_joined %>%
  rowwise() %>%
  filter(length(intersect(top_genes1, top_genes2)) > 0) %>%
  ungroup()

outdff = data.frame("Trait" = trait,
                    "Hits_Combo" = nrow(top_joined_with_overlap),
                    "Hits_ENCODE_E2G" = nrow(top1),
                    "Hits_PoPS" = nrow(top2),
                    "All" = length(uu))


write.table(outdff,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_PoPS_Finemap_Variant_part2/",
                          trait, "_variant_overlap_ENCODE_PoPS.txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)
