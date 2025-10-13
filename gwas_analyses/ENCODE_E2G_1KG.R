library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--numchr", type="character", default = "data/traits.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

numchr = as.numeric(opt$numchr)

bimdf = data.frame(fread(paste0("/data/deyk/extras/1000G_BIMS_hg38/1000G.EUR.QC.", numchr, ".bim")))

gr2 = GRanges(seqnames = paste0("chr", bimdf[,1]),
              ranges = IRanges(start= bimdf[,4] - 1, end = bimdf[,4] + 1),
              snp_name = bimdf[,2])

pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions"
logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
outnames = as.character(unlist(sapply(logreg_dnase_biosamples, function(xx) return(strsplit(xx, "_DNase")[[1]][1]))))
logreg_dnase_biosamples_of_interest = logreg_dnase_biosamples

pooled_dff_list = c()
for(numl in 1:length(logreg_dnase_biosamples_of_interest)){
  preds_tabb = data.frame(fread(paste0(pred_cell, "/", logreg_dnase_biosamples_of_interest[numl])))
  preds_tabb = preds_tabb[which(preds_tabb$X.chr %in% paste0("chr", 1:22) == T), ]
  chr_number2 = as.numeric(unlist(sapply(preds_tabb$X.chr, function(xx) return(as.character(strsplit(xx, "chr")[[1]][2])))))

  gr3 = GRanges(seqnames = paste0("chr", chr_number2),
                ranges = IRanges(start= preds_tabb$start, end = preds_tabb$end),
                gene = preds_tabb$TargetGene,
                score = preds_tabb$Score)

  cc=findOverlaps(gr2, gr3, type = "any", select = "all")
  if(length(cc) > 0){
    dff = cbind.data.frame(gr2[queryHits(cc)],
                           gr3[subjectHits(cc)],
                           logreg_dnase_biosamples_of_interest[numl])
    pooled_dff_list[[numl]] = dff
  }
  cat("We are at cell-type:", numl, "\n")
}

pooled_dff = do.call(rbind, pooled_dff_list)
pooled_dff$BP_hg38 = (pooled_dff$start + pooled_dff$end)/2
pooled_dff2 = pooled_dff
colnames(pooled_dff2)[14] = "Biosample"

write.table(pooled_dff2,
            file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG/",
                          "ENCODE_E2G_gene_Allcelltypes.", numchr, ".txt"), row.names = F,
            col.names = T, sep = "\t", quote=F)






