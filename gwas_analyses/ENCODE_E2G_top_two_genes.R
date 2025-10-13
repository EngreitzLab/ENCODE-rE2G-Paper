library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)
library(dplyr)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--biosample", type="character", default = "data/traits.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

biosample = opt$biosample
# biosample = paste0(biosample, "_DNaseOnly_features.predictions.tsv.gz")
# biosample = "tibial_artery_ENCDO451RUA_ENCFF217PGC"


# biosamples = read.table("/data/deyk/ENCODE/ENCODE_E2G/bio_modes.txt", header=F)[,1]
# biosamples2 = paste0(as.character(sapply(biosamples, function(xx) return(strsplit(xx, "LogReg_DNase_")[[1]][2]))),
#                      "_DNaseOnly_features.predictions.tsv.gz")
# biosamples3 = paste0(as.character(sapply(biosamples, function(xx) return(strsplit(xx, "LogReg_DNase_")[[1]][2]))))
# write.table(biosamples3, file = "/data/deyk/ENCODE/ENCODE_E2G/bio_modes2.txt",
#             row.names=F, col.names=F, sep = "\t", quote=F)


library(tidyr)

ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG")
for(numchr in 1:22){
  tabb = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG",
                          "/",  "ENCODE_E2G_gene_Allcelltypes.", numchr, ".txt")))
  tabb = tabb[which(tabb$Biosample == paste0(biosample)), ]
  xx = tapply(tabb$score, paste0(tabb$snp_name, ":", tabb$gene), max)
  rsids = as.character(sapply(names(xx), function(y) return(strsplit(y, ":")[[1]][1])))
  genes = as.character(sapply(names(xx), function(y) return(strsplit(y, ":")[[1]][2])))
  outdff = cbind.data.frame(rsids, genes, xx)
  colnames(outdff) = c("SNP", "Gene", "value")
  rr = outdff %>% pivot_wider(names_from = Gene, values_from = value)
  rr = as.data.frame(rr)
  rr2 = rr[,-1]
  rownames(rr2) = rr[,1]


  genesXX = apply(rr2, 1, function(x){
    cc = which(!is.na(x))
    nn = names(x)[!is.na(x)]
    yy = x[!is.na(x)]
    zz = nn[order(yy, decreasing = T)[1:min(2, length(cc))]]
    if(length(zz) == 1){
      zz = c(zz, NA)
    }
    return(zz)
  })

  final_dff = cbind.data.frame(rownames(rr2), t(genesXX))
  colnames(final_dff) = c("SNP", "Gene1", "Gene2")
  if(!dir.exists(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG_records_top2/",
                        opt$biosample))){
    dir.create(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG_records_top2/",
                      opt$biosample))
  }
  write.table(final_dff, file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG_records_top2/",
                                       opt$biosample, "/",  "ENCODE_E2G_1KG_records.", numchr, ".txt"),
              row.names = F, col.names = T, sep = "\t", quote=F)
  cat("We are at chr:", numchr, "\n")
}
