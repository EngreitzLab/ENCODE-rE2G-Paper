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

library(tidyr)

ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG")

chunk <- function(x,n){
  numOfVectors <- floor(length(x)/n)
  elementsPerVector <- c(rep(n,numOfVectors-1),n+length(x) %% n)
  elemDistPerVector <- rep(1:numOfVectors,elementsPerVector)
  split(x,factor(elemDistPerVector))
}

tabb_pre = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG",
                                    "/", trait, "_Variant_gene.txt")))
mm = paste0(tabb_pre$snp_name, ":", tabb_pre$gene)
tabb = tabb_pre[match(unique(mm), mm), ]

outdff = cbind.data.frame(tabb$snp_name, tabb$gene, tabb$pops_score)
colnames(outdff) = c("SNP", "Gene", "value")

usnps = unique(outdff$SNP)

cc = chunk(usnps,floor(length(usnps)/100))

final_dff = list()
for(numc in 1:length(cc)){
  outdff2 = outdff[which(outdff$SNP %in% cc[[numc]]), ]
  rr = outdff2 %>% pivot_wider(names_from = Gene, values_from = value)
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

  final_dff_chunk = cbind.data.frame(rownames(rr2), t(genesXX))
  colnames(final_dff_chunk) = c("SNP", "Gene1", "Gene2")
  final_dff[[numc]] = final_dff_chunk
  cat("We are at chunk:", numc, "\n")
}

final_dff_frame = do.call(rbind, final_dff)

write.table(final_dff_frame, file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG_records_top2/", trait,
                                     "_PoPS_1KG_records.", ".txt"),
            row.names = F, col.names = T, sep = "\t", quote=F)




