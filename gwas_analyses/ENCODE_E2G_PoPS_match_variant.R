library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)
library(dplyr)

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

pops_tabb_trait = pops_tabb[which(pops_tabb$trait == trait), ]
gr1_pops = GRanges(seqnames = paste0("chr", pops_tabb_trait$chromosome),
                ranges = IRanges(start=pmax(pops_tabb_trait$tss-500000, 0),
                                 end = pops_tabb_trait$tss + 500000),
                gene = pops_tabb_trait$gene,
                score = pops_tabb_trait$pops_score)


bimlist = list()
for(numchr in 1:22){
  xx = data.frame(fread(paste0("/data/deyk/kushal/extras/BIMS_hg38/", "1000G.EUR.QC.", numchr, ".bim")))
  colnames(xx) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
  bimlist[[numchr]] = xx
  cat("Processing bimfile for chr:", numchr, "\n")
}
names(bimlist) = paste0("chr", 1:22)
bimfile_merged = do.call(rbind, bimlist)

gr2 = GRanges(seqnames = paste0("chr", bimfile_merged$CHR),
                ranges = IRanges(start=bimfile_merged$BP-1,
                                 end = bimfile_merged$BP),
              snp_name = bimfile_merged$SNP)

cc=findOverlaps(gr2, gr1, type = "any", select = "all")
xtabb = cbind.data.frame(gr2[queryHits(cc)], gr1[subjectHits(cc)])
colnames(xtabb) = c("var.seqnames", "var.start", "var.end", "var.width", "var.strand",
                    "snp_name", "seqnames", "start", "end", "width", "strand", "gene", "score")

cc2=findOverlaps(gr2, gr1_pops, type = "any", select = "all")
xtabb2 = cbind.data.frame(gr2[queryHits(cc2)], gr1_pops[subjectHits(cc2)])
colnames(xtabb2) = c("var.seqnames", "var.start", "var.end", "var.width", "var.strand",
                    "snp_name", "seqnames", "start", "end", "width", "strand", "gene", "score")

top_list = list()
for(numchr in 1:22){
  xtabb_chr = xtabb[which(xtabb$var.seqnames == paste0("chr", numchr)), ]
  top1 <- xtabb_chr %>%
    group_by(snp_name) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    select(snp_name, gene_top1 = gene)

  xtabb2_chr = xtabb2[which(xtabb2$var.seqnames == paste0("chr", numchr)), ]

  top2 <- xtabb2_chr %>%
    group_by(snp_name) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    select(snp_name, gene_top2 = gene)

  top_joined <- inner_join(top1, top2, by = "snp_name") %>%
    filter(gene_top1 == gene_top2)

  top_list[[numchr]] = as.data.frame(top_joined)
}
top_dff = do.call(rbind, top_list)

write.table(top_dff,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_PoPS_Variant/",
                          trait, "_variant_overlap_ENCODE_PoPS.txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)


