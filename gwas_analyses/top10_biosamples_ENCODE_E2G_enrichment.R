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

enrvec = c()
for(numl in 1:length(outnames)){
  preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                                       "/", outnames[numl])))
  preds_tabb$chr_num = as.integer(as.character(sapply(preds_tabb$X.chr, function(zz) return(strsplit(zz, "chr")[[1]][2]))))
  preds_tabb = preds_tabb[(preds_tabb$chr_num %in% 1:22), ]

  gr1 = GRanges(seqnames = preds_tabb$chr_num,
                ranges = IRanges(start= preds_tabb$start - 1, end = preds_tabb$end + 1),
                targetgene = preds_tabb$TargetGene,
                celltype = preds_tabb$CellType,
                score = preds_tabb$Score)

  gr2 = GRanges(seqnames = finemap_tabb$CHR,
                ranges = IRanges(start= finemap_tabb$BP - 1, end = finemap_tabb$BP + 1),
                snpname = finemap_tabb$SNP)

  cc=findOverlaps(gr2, gr1, type = "any", select = "all")
  frac1 = length(unique(queryHits(cc)))/nrow(finemap_tabb)

  gr3 = GRanges(seqnames = finemap_tabb2$CHR,
                ranges = IRanges(start= finemap_tabb2$BP - 1, end = finemap_tabb2$BP + 1),
                snpname = finemap_tabb2$SNP)

  cc2=findOverlaps(gr3, gr1, type = "any", select = "all")
  frac2 = length(unique(queryHits(cc2)))/nrow(finemap_tabb2)

  if(frac1 == 0 | frac2 == 0){
    enrichment = 1
  }else{
    enrichment = frac2/frac1
  }

  enrvec = c(enrvec, enrichment)
  cat("We are at biosample:", numl, "\n")
}

outdff = cbind.data.frame(outnames, enrvec)
colnames(outdff) = c("Biosample", "ENR")
write.table(outdff,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Top10_Biosamples_Per_Trait/",
                          trait, "_biosamples_ranked_enrichment.txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)



