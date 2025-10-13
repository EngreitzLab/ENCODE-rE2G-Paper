library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--tissue", type="character", default = "data/tissue.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

tissue = opt$tissue

# tissues = list.files("/n/groups/price/kushal/ENCODE/data/GTEx_finemapping/MaxPIP_BY_TISSUE", pattern = ".rsid")
# outnames = as.character(unlist(sapply(tissues, function(xx) return(strsplit(xx, ".rsid")[[1]][1]))))
# write.table(outnames, file = "/n/groups/price/kushal/ENCODE/code/Figure4/ENCODE_E2G_biosamples/tissues_record.txt",
#                     row.names = F, col.names = F, sep = "\t", quote=F)
#

# tissue = "Liver.rsid"
maxcpp_tabb = data.frame(fread(paste0("/data/deyk/GTEx/SuSIE_finemapping/MaxPIP_BY_TISSUE", "/",
                                       tissue)))

bimpooltabb1 = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0("/data/deyk/extras/1000G_BIMS_hg38/1000G.EUR.QC.", numchr, ".bim")))
  bimpooltabb1 = rbind(bimpooltabb1, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}
colnames(bimpooltabb1) = c("CHR", "SNP", "CM", "BP", "A1", "A2")

common_snps = intersect(maxcpp_tabb$SNP, bimpooltabb1$SNP)
finemap_tabb = bimpooltabb1[match(common_snps, bimpooltabb1$SNP), ]
finemap_tabb$MaxCPP = maxcpp_tabb$MaxCPP[match(common_snps, maxcpp_tabb$SNP)]


pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions"
outnames = list.files(pred_cell, pattern = ".tsv.gz")

pooled_dff_list = list()

for(numl in 1:length(outnames)){
  preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                                       "/", outnames[numl])))

  gr1 = GRanges(seqnames = preds_tabb$X.chr,
                ranges = IRanges(start= preds_tabb$start, end = preds_tabb$end),
                targetgene = preds_tabb$TargetGene,
                celltype = preds_tabb$CellType,
                score = preds_tabb$Score)

  gr2 = GRanges(seqnames = paste0("chr", finemap_tabb$CHR),
                ranges = IRanges(start= finemap_tabb$BP - 1, end = finemap_tabb$BP + 1),
                snpname = finemap_tabb$SNP,
                maxpip = finemap_tabb$MaxCPP)

  cc=findOverlaps(gr2, gr1, type = "any", select = "all")

  if(length(cc) > 0){
    dff = cbind.data.frame(seqnames(gr2)[queryHits(cc)],
                           start(ranges(gr2))[queryHits(cc)],
                           end(ranges(gr2))[queryHits(cc)],
                           gr2$snpname[queryHits(cc)],
                           tissue,
                           gr2$maxpip[queryHits(cc)],
                           gr1$targetgene[subjectHits(cc)],
                           gr1$celltype[subjectHits(cc)],
                           gr1$score[subjectHits(cc)])
    colnames(dff) = c("CHR", "start", "end", "SNP", "Tissue",  "MaxPIP", "TargetGene", "Cell-type", "Score")
    pooled_dff_list[[numl]] = dff
  }
  cat("We are at biosample:", outnames[numl], "\n")
}

pooled_dff = do.call(rbind, pooled_dff_list)


fwrite(pooled_dff, file = paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GTEx_ENCODE_rE2G/", "/",
                                 tissue,  "_GTEX_ENCODE_E2G_all_biosamples.txt"),
       row.names=F, col.names=T, sep = "\t", quote=F)
