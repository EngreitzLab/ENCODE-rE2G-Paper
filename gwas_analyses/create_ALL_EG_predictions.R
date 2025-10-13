
Ramil_ALL_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/EpiRaction/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/Ramil_ALL_EG_predictions.txt")
{
  library(data.table)
  ramil_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(ramil_dnase_biosamples)

  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
out1 = Ramil_ALL_EG_preds()

ABC2024_ALL_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/ABC2024_ALL_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ABC2024_ALL_EG_preds()


ENCODE_E2G2024_ALL_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/ENCODE_E2G2024_ALL_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ENCODE_E2G2024_ALL_EG_preds()


ENCODE_E2G2024_HiC_ALL_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/ENCODE_E2G2024_HiC_ALL_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ENCODE_E2G2024_HiC_ALL_EG_preds()


Baseline_ALL_EG_preds <- function(
    pred_type = "dist_to_tss",
    preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/Baseline_ALL_dist_to_tss_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = read.table(preds_file, header = F)[,1]
  keep_eids = unique(logreg_dnase_biosamples)

  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0("/lila/data/deyk/ENCODE/Baseline_Preds_2024/",
                          keep_eids[numl], "/", pred_type, ".tsv.gz"))){

      preds_tabb = data.frame(fread(paste0("/lila/data/deyk/ENCODE/Baseline_Preds_2024/",
                                           keep_eids[numl], "/", pred_type, ".tsv.gz")))
      if(pred_type == "dist_to_tss"){
        preds_tabb = preds_tabb[which(preds_tabb$Score < 53966), ]
      }
      if(pred_type == "dist_to_gene"){
        preds_tabb = preds_tabb[which(preds_tabb$Score < 44840), ]
      }
      if(pred_type == "DHS_reads_by_dist_to_tss"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 3.149175673932778e-4), ]
      }
      if(pred_type == "DHS_reads_by_dist_to_tss_norm"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 7.391383676221246e-4), ]
      }
      if(pred_type == "H3K27ac_reads_by_dist_to_tss"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 2.2513369374164912e-4), ]
      }
      if(pred_type == "H3K27ac_reads_by_dist_to_tss_norm"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 8.201190258983458e-4), ]
      }
      if(pred_type == "nearest_expressed_gene"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 0), ]
      }
      if(pred_type == "nearest_expressed_tss"){
        preds_tabb = preds_tabb[which(preds_tabb$Score > 0), ]
      }

      tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end,
                               preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb_list[[numl]] = tabb2
      cat("We are at blood sample:", numl, "\n")
    }
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


out1 = Baseline_ALL_EG_preds(pred_type = "dist_to_tss",
                             preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
                             output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/Baseline_ALL_dist_to_tss_EG_predictions.txt")

out2 = Baseline_ALL_EG_preds(pred_type = "nearest_expressed_gene",
                             preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
                             output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/Baseline_ALL_nearest_expressed_gene_EG_predictions.txt")


opentargets_L2G_file = read.csv(paste0("/data/deyk/opentargets/opentargets_pops/complex_traits_v2_OpenTargets_L2G.Cleaned.study_traits_pops_29Sep24.csv"))
opentargets_set_file = read.csv(paste0("/data/deyk/opentargets/opentargets_pops/complex_traits_v2_OpenTargets_Credible_Sets.study_traits_pops_29Sep24.csv"))

opentargets_L2G_file$key_L2G_file = paste0(opentargets_L2G_file$study_id, "_", opentargets_L2G_file$lead_variant_id)
opentargets_set_file$key_L2G_file = paste0(opentargets_set_file$study_id, "_", opentargets_set_file$lead_variant_id)

matched_keys = intersect(opentargets_L2G_file$key_L2G_file, opentargets_set_file$key_L2G_file)

pooled_tabb_list = list()
for(numl in 1:length(matched_keys)){
  file1 = opentargets_L2G_file[which(opentargets_L2G_file$key_L2G_file == matched_keys[numl]), ]
  file2 = opentargets_set_file[which(opentargets_set_file$key_L2G_file == matched_keys[numl]), ]
  xx = tapply(file1$L2G, file1$gene_name, max)
  chrval = paste0("chr", strsplit(unique(file1$lead_variant_id), ":")[[1]][1])
  pooled_tabb_list[[numl]] = cbind.data.frame(chrval, min(file2$cred_set_min_pos), max(file2$cred_set_max_pos),
                                              rownames(xx), as.numeric(xx))
  cat("We are at variant and trait combination:", numl, " of ", length(matched_keys), " instances \n")
}

pooled_dff = do.call(rbind, pooled_tabb_list)
colnames(pooled_dff) = c("chr", "start", "end", "TargetGene", "Score")
fwrite(pooled_dff, "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/OpenTargets_L2G_predictions.txt",
       sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


ll = list.files("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/All_EG_predictions/")
oo = cbind.data.frame(ll, "ALL", paste0(ll, ".prec_recall"))
write.table(oo, file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/finemap_precision_recall5.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)
