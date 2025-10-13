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

dff2 = data.frame(fread(paste0("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets",
                               "/", trait, "/", "variant.list.txt")))
dff2 = dff2[which(dff2$pip > 0.1), ]
uu = unique(dff2$rsid)


bimpooltabb1 = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0("/data/deyk/extras/1000G_BIMS_hg38/1000G.EUR.QC.", numchr, ".bim")))
  bimpooltabb1 = rbind(bimpooltabb1, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}
colnames(bimpooltabb1) = c("CHR", "SNP", "CM", "BP", "A1", "A2")

finemap_tabb = bimpooltabb1
finemap_tabb2 = finemap_tabb[match(intersect(uu, finemap_tabb$SNP), finemap_tabb$SNP), ]


pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions"
outnames = list.files(pred_cell, pattern = "thresholded_predictions.tsv.gz")

top_biosamples = read.table(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Top10_Biosamples_Per_Trait/",
                                   trait, "_biosamples_ranked_enrichment.txt"), header=T)
outnames2 = top_biosamples[order(top_biosamples[,2], decreasing = T)[1:10], 1]

preds_tabb_top10 = list()
for(numl in 1:length(outnames2)){
  preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                                       "/", outnames2[numl])))
  preds_tabb$chr_num = as.integer(as.character(sapply(preds_tabb$X.chr, function(zz) return(strsplit(zz, "chr")[[1]][2]))))
  preds_tabb = preds_tabb[(preds_tabb$chr_num %in% 1:22), ]
  preds_tabb_top10[[numl]] = preds_tabb
}
preds_tabb_top10_dff = do.call(rbind, preds_tabb_top10)
eg_preds_finemap = preds_tabb_top10_dff[, c("X.chr", "start", "end", "TargetGene", "Score")]
colnames(eg_preds_finemap) = c("chr", "start", "end", "TargetGene", "Score")


eg_preds = data.frame(fread("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/ENCODE_E2G2024_ALL_EG_predictions.txt"))


gr1 = GRanges(seqnames = eg_preds$chr,
              ranges = IRanges(start= eg_preds$start - 1, end = eg_preds$end + 1),
              targetgene = eg_preds$TargetGene,
              score = eg_preds$Score)

gr2 = GRanges(seqnames = paste0("chr", finemap_tabb2$CHR),
              ranges = IRanges(start= finemap_tabb2$BP - 1, end = finemap_tabb2$BP + 1),
              snpname = finemap_tabb2$SNP)

cc=findOverlaps(gr2, gr1, type = "any", select = "all")

outdff1 = cbind.data.frame(seqnames(gr1)[subjectHits(cc)],
                           start(ranges(gr1))[subjectHits(cc)],
                           end(ranges(gr1))[subjectHits(cc)],
                           gr1$targetgene[subjectHits(cc)],
                           gr1$score[subjectHits(cc)],
                           gr2$snpname[queryHits(cc)])
colnames(outdff1) = c("CHR", "start", "end", "TargetGene", "Score", "SNP")



gr3 = GRanges(seqnames = eg_preds_finemap$chr,
              ranges = IRanges(start= eg_preds_finemap$start - 1, end = eg_preds_finemap$end + 1),
              targetgene = eg_preds_finemap$TargetGene,
              score = eg_preds_finemap$Score)

gr2 = GRanges(seqnames = paste0("chr", finemap_tabb2$CHR),
              ranges = IRanges(start= finemap_tabb2$BP - 1, end = finemap_tabb2$BP + 1),
              snpname = finemap_tabb2$SNP)

cc2=findOverlaps(gr2, gr3, type = "any", select = "all")

outdff2 = cbind.data.frame(seqnames(gr3)[subjectHits(cc2)],
                           start(ranges(gr3))[subjectHits(cc2)],
                           end(ranges(gr3))[subjectHits(cc2)],
                           gr3$targetgene[subjectHits(cc2)],
                           gr3$score[subjectHits(cc2)],
                           gr2$snpname[queryHits(cc2)])
colnames(outdff2) = c("CHR", "start", "end", "TargetGene", "Score", "SNP")


pops_tabb = data.frame(fread(paste0("/data/deyk/kushal/PolyPerturb/codes/PoPS/supplemental_data/", "PoPS_FullResults.txt.gz")))
pops_tabb_trait = pops_tabb[which(pops_tabb$trait == trait), ]

pops_score1 = mean(pops_tabb_trait$pops_score[match(intersect(unique(outdff1$TargetGene), pops_tabb_trait$gene),
                      pops_tabb_trait$gene)])
pops_score2 = mean(pops_tabb_trait$pops_score[match(intersect(unique(outdff2$TargetGene), pops_tabb_trait$gene),
                                      pops_tabb_trait$gene)])

merged_dff = data.frame("Gene_set" = c("ENCODE_E2G_ALL", "ENCODE_E2G_Top10"),
                        "PoPS_score" = c(pops_score1, pops_score2))

write.table(merged_dff,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/PoPS_scores_Top10_Biosamples_Per_Trait/",
                          trait, "_pops_scores_top10biosamples.txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)


