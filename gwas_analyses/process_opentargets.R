library(data.table)
library(R.utils)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

numl <- as.numeric(toString(args[1]))

library(data.table)
merged_dff = list()

for(numchr in 1:22){
  dff = data.frame(fread(paste0("/data/deyk/opentargets/l2g_OpenTargets/OpenTargets.chr", numchr, ".data.gz")))
  dff2 = dff[which(dff[,8] > numl/10), ]
  merged_dff[[numchr]] = dff2
  cat("We are at chr:", numchr, "\n")
}
merged_dff_tabb = do.call(rbind, merged_dff)

fwrite(merged_dff_tabb, file = paste0("/data/deyk/opentargets/opentargets_filtered/opentargets_filtered_score", numl*10, ".txt"),
       row.names=F, col.names=T, sep = "\t", quote=F)



