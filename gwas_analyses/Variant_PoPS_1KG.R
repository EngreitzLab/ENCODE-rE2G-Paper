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

pops_traits_results = data.frame(fread("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/PoPS_FullResults.txt.gz"))
pops_traits_results_temp = pops_traits_results[which(pops_traits_results$trait == trait), ]

pops_scores = pops_traits_results_temp$pops_score
names(pops_scores) = pops_traits_results_temp$gene

pops_scores_order = rank(-pops_scores)
names(pops_scores_order) = pops_traits_results_temp$gene

pops_traits_results_temp$pops_rank = pops_scores_order

bimpooltabb1 = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0("/data/deyk/extras/1000G_BIMS_hg38/1000G.EUR.QC.", numchr, ".bim")))
  bimpooltabb1 = rbind(bimpooltabb1, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}
colnames(bimpooltabb1) = c("CHR", "SNP", "CM", "BP", "A1", "A2")


bimpooltabb2 = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0("/data/deyk/extras/1000G_BIMS/1000G.EUR.QC.", numchr, ".bim")))
  bimpooltabb2 = rbind(bimpooltabb2, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}
colnames(bimpooltabb2) = c("CHR", "SNP", "CM", "BP", "A1", "A2")

cc = intersect(bimpooltabb2$SNP, bimpooltabb1$SNP)

merged_bimpooltabb = cbind.data.frame(bimpooltabb1$CHR[match(cc, bimpooltabb1$SNP)],
                                      bimpooltabb1$SNP[match(cc, bimpooltabb1$SNP)],
                                      bimpooltabb1$BP[match(cc, bimpooltabb1$SNP)],
                                      bimpooltabb2$BP[match(cc, bimpooltabb2$SNP)],
                                      bimpooltabb1$A1[match(cc, bimpooltabb1$SNP)],
                                      bimpooltabb1$A2[match(cc, bimpooltabb1$SNP)])

colnames(merged_bimpooltabb) = c("CHR", "SNP", "BP_hg38", "BP_hg19", "A1", "A2")


gr2 = GRanges(seqnames = paste0("chr", merged_bimpooltabb$CHR),
              ranges = IRanges(start= merged_bimpooltabb$BP_hg19 - 1, end = merged_bimpooltabb$BP_hg19 + 1),
              snp_name = merged_bimpooltabb$SNP)

gr3 = GRanges(seqnames = paste0("chr", pops_traits_results_temp$chromosome),
              ranges = IRanges(start= pops_traits_results_temp$start - 500000, end = pops_traits_results_temp$end + 500000),
              pops_score = pops_traits_results_temp$pops_score,
              pops_rank = pops_traits_results_temp$pops_rank,
              gene = pops_traits_results_temp$gene,
              ensembl = pops_traits_results_temp$ensgid)

cc=findOverlaps(gr2, gr3, type = "any", select = "all")

if(length(cc) > 0){
  dff = cbind.data.frame(gr2[queryHits(cc)],
                         gr3[subjectHits(cc)])
}
dff$BP_hg19 = (dff$start + dff$end)/2

write.table(dff,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG/",
                          trait, "_Variant_gene.txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)








