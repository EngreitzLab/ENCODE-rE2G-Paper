ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2_trial2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 91, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2_trial2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/EE2calc_traits_x_celltypes_2024_trial3.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)



ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_only_ENCODE_E2G_NoPoPS_top2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 91, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_only_ENCODE_E2G_NoPoPS_top2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/EE2calc_noPoPS_traits_x_celltypes_2024_trial1.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)


ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2_newbase")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 91, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2_newbase", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/EE2calc_newbase_traits_x_celltypes_2024_trial1.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)



ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 91, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/EE2calc_traits_x_celltypes_2024_trial1.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)


ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS_top2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 91, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS_top2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recallcalc_traits_x_celltypes_2024_trial1.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)




ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 93, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/EE2calc_traits_x_celltypes_2024_trial2.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)


ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 93, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recallcalc_traits_x_celltypes_2024_trial2.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)


ll = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS2")
ll2 = as.character(sapply(ll, function(z) return(strsplit(strsplit(z, "_thresholded")[[1]][1], "predictions_")[[1]][2])))

pooled_dff = matrix(0, 93, length(ll))
for(numl in 1:length(ll)){
  dff = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recall_ENCODE_E2G_PoPS2", "/",
                                ll[numl])))
  pooled_dff[,numl] = dff[,2]
  cat("We are at biosample:", numl, "\n")
}
pooled_dff = as.data.frame(pooled_dff)
rownames(pooled_dff) = dff[,1]
colnames(pooled_dff) = ll2

write.table(pooled_dff, file = "/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Recallcalc_traits_x_celltypes_2024_trial1.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)

