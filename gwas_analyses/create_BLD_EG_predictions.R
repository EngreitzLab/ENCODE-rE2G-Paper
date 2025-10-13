

Ramil_BLD_EG_preds <- function(
                       pred_cell="/data/deyk/ENCODE/EpiRaction/thresholded_predictions/",
                       output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/Ramil_Blood_EG_predictions.txt")
{
  library(data.table)
  ramil_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(ramil_dnase_biosamples[c(grep("Blood", ramil_dnase_biosamples))])

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
out1 = Ramil_BLD_EG_preds()

ABC2024_BLD_EG_preds <- function(
        pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions/",
        output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ABC2024_Blood_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
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

out1 = ABC2024_BLD_EG_preds()


ENCODE_E2G2024_BLD_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G2024_Blood_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
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

out1 = ENCODE_E2G2024_BLD_EG_preds()


ENCODE_E2G2024_HiC_BLD_EG_preds <- function(
    pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions/",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G2024_HiC_Blood_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
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

out1 = ENCODE_E2G2024_HiC_BLD_EG_preds()


Baseline_BLD_EG_preds <- function(
    pred_type = "dist_to_tss",
    preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/Baseline_Blood_dist_to_tss_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = read.table(preds_file, header = F)[,1]
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])

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


out1 = Baseline_BLD_EG_preds(pred_type = "dist_to_tss",
                             preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
                             output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/Baseline_Blood_dist_to_tss_EG_predictions.txt")

out2 = Baseline_BLD_EG_preds(pred_type = "nearest_expressed_gene",
                             preds_file="/lila/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
                             output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/Baseline_Blood_nearest_expressed_gene_EG_predictions.txt")


EpiMap_BLD_EG_preds <- function(scores,
                                pred_cell="/lila/data/deyk/ENCODE/EpiMap/EG_predictions/",
                                key_file="/lila/data/deyk/ENCODE/EpiMap/main_metadata_table.tsv",
                                output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/EpiMap_Blood_EG_predictions.txt"){
  library(data.table)
  biosamples_listed = as.character(sapply(as.character(sapply(list.files(pred_cell, pattern = "linking"),
                                                              function(x) return(strsplit(x, "_")[[1]][3]))),
                                          function(xx) return(strsplit(xx, "[.]")[[1]][1])))
  epimap_biosamples = read.delim(paste0(key_file))
  keep_eids = epimap_biosamples$id[c(grep("Blood & T-cell", epimap_biosamples$GROUP),
                                     grep("HSC & B-cell", epimap_biosamples$GROUP),
                                     grep("Lymphoblastoid", epimap_biosamples$GROUP))]

  keep_eids = intersect(keep_eids, biosamples_listed)

  pooled_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0("/lila/data/deyk/ENCODE/EpiMap/EG_predictions/",
                                         "linking_collated_", keep_eids[numl],".bed.gz")))
    preds_tabb = preds_tabb[which(preds_tabb$Score > 0.0045), ]
    tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_list[[numl]] = tabb2
    cat("Read file", numl, "out of ", length(keep_eids), " files \n")
  }
  pooled_tabb = do.call(rbind, pooled_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = EpiMap_BLD_EG_preds()


pred_cell = "/data/deyk/ENCODE/GraphReg"
tabb_pre1 = data.frame(fread(paste0(pred_cell, "/", "EG_preds_", "K562", "_hg38_GraphReg_FDR_1_L_8_full_predictions.tsv.gz")))
tabb_pre1 = tabb_pre1[which(tabb_pre1$GraphReg.Score > 0.001368199), ]
tabb_pre2 = data.frame(fread(paste0(pred_cell, "/", "EG_preds_", "GM12878", "_hg38_GraphReg_FDR_1_L_8_full_predictions.tsv.gz")))
tabb_pre2 = tabb_pre2[which(tabb_pre2$GraphReg.Score > 0.001368199), ]

tabb2 = rbind.data.frame(tabb_pre1, tabb_pre2)
tabb2 = tabb2[, c("chr", "start", "end", "TargetGene", "GraphReg.Score")]
colnames(tabb2) = c("chr", "start", "end", "TargetGene", "GraphReg.Score")

fwrite(tabb2, "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/GraphReg2023_Blood_EG_predictions.txt",
       sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


ll = list.files("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/")
ll = c("ENCODE_E2G2024_Blood_EG_predictions.txt",
       "ENCODE_E2G2024_HiC_Blood_EG_predictions.txt",
       "Baseline_Blood_nearest_expressed_gene_EG_predictions.txt",
       "GraphReg2023_Blood_EG_predictions.txt",
       "Baseline_Blood_dist_to_tss_EG_predictions.txt",
       "EpiMap_Blood_EG_predictions.txt",
       "Ramil_Blood_EG_predictions.txt",
       "ABC2024_Blood_EG_predictions.txt")
oo = cbind.data.frame(ll, "BLD", paste0(ll, ".prec_recall"))
write.table(oo, file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/finemap_precision_recall2.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################
###################################################. 2023. #########################################################################

ABC_DNaseonly_default_calc_preds <- function(scores,
                                       pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/ABC_2022",
                                       output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/ABC_DNaseonly_default_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell, pattern = "DNase-only")
  biosamples_file = read.delim(paste0(pred_cell, "/", "DNase-only.default_metadata.101322.qc.tsv"))
  biosamples = biosamples_file$Biosample_id[which(biosamples_file$default == "True")]
  biosamples = paste0(biosamples, "_DNase-only")
  keep_eids = biosamples[unique(unlist(lapply(c("T-cell", "T_cell", "T-helper", "CD4", "CD8",  "K562",
                                                "hemato", "GM", "B_cell", "B-cell", "CD14", "CD1c", "macrophage"),
                                              function(x) return(grep(x, biosamples)))))]

  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))){
      if (file.info(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))$size != 0){
        preds_tabb = data.frame(fread(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz")))
        tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$ABC.Score)
        colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
        pooled_tabb_list[[numl]] = tabb2
        cat("We are at biosample:", numl, "\n")
      }
    }
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}



ABC_DNaseH3K27ac_calc_preds <- function(scores,
                                       pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/ABC_2022",
                                       output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/ABC_DNaseH3K27ac_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  biosamples_file = read.delim(paste0(pred_cell, "/", "DNase_H3K27ac.default_metadata_full.101322.qc.tsv"))
  biosamples = paste0(biosamples_file$Biosample)

  keep_eids = biosamples[unique(unlist(lapply(c("T-cell", "T_cell", "T-helper", "CD4", "CD8",  "K562",
                                                "hemato", "GM", "B_cell", "B-cell", "CD14", "CD1c", "macrophage"),
                                              function(x) return(grep(x, biosamples)))))]

  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))){
      if (file.info(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))$size != 0){
        preds_tabb = data.frame(fread(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz")))
        tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$ABC.Score)
        colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
        pooled_tabb_list[[numl]] = tabb2
        cat("We are at biosample:", numl, "\n")
      }
    }
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}



ABC_DNaseH3K27ac_default_calc_preds <- function(scores,
                                  pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/ABC_2022",
                                  output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/ABC_DNaseH3K27ac_default_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  biosamples_file = read.delim(paste0(pred_cell, "/", "DNase_H3K27ac.default_metadata_full.101322.qc.tsv"))
  biosamples = biosamples_file$Biosample[which(biosamples_file$default == "True")]

  keep_eids = biosamples[unique(unlist(lapply(c("T-cell", "T_cell", "T-helper", "CD4", "CD8",  "K562",
                                                "hemato", "GM", "B_cell", "B-cell", "CD14", "CD1c", "macrophage"),
                                              function(x) return(grep(x, biosamples)))))]

  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))){
      if (file.info(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))$size != 0){
        preds_tabb = data.frame(fread(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz")))
        tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$ABC.Score)
        colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
        pooled_tabb = rbind(pooled_tabb, tabb2)
        cat("We are at biosample:", numl, "\n")
      }
    }
  }
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}



EPIraction_default_calc_preds <- function(scores,
                                    pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/EPIraction/preds",
                                    output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/EPIraction_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  keep_eids = ll[unique(unlist(lapply(c("Blood"), function(x) return(grep(x, ll)))))]

  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/", keep_eids[numl]))){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
      cat("We are at biosample:", numl, "\n")
    }
  }

  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


delta_gencode_calc_preds <- function(scores,
                               pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/HiC/delta_gencode",
                               output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/delta_gencode_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  keep_eids = ll[unique(unlist(lapply(c("cd4", "cd8", "mono", "K562", "GM12878", "NK"), function(x) return(grep(x, ll)))))]
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/", keep_eids[numl]))){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }
  }
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


delta_refseq_calc_preds <- function(scores,
                                     pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/HiC/delta_refseq",
                                     output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/delta_refseq_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  keep_eids = ll[unique(unlist(lapply(c("cd4", "cd8", "mono", "K562", "GM12878", "NK"), function(x) return(grep(x, ll)))))]
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/", keep_eids[numl]))){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }
  }
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


hiccups_refseq_calc_preds <- function(scores,
                                    pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/HiC/hiccups_refseq",
                                    output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/hiccups_refseq_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  keep_eids = ll[unique(unlist(lapply(c("cd4", "cd8", "mono", "K562", "GM12878", "NK"), function(x) return(grep(x, ll)))))]
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/", keep_eids[numl]))){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }
  }
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


hiccups_gencode_calc_preds <- function(scores,
                                      pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/HiC/hiccups_gencode",
                                      output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/hiccups_gencode_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  keep_eids = ll[unique(unlist(lapply(c("cd4", "cd8", "mono", "K562", "GM12878", "NK"), function(x) return(grep(x, ll)))))]
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/", keep_eids[numl]))){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }
  }
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}




BingABC_BLD_EG_preds <- function(
                         pred_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/2022_BingRen_scATAC_EnhancerPredictionsFull.txt.gz",
                         key_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/BingABC_biosamples.csv",
                         output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/BingABC_Blood_EG_predictions.txt")
{
  library(data.table)
  bingabc_biosamples = read.csv(paste0(key_file))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb_pre = tabb_pre[which(tabb_pre$Score > 0.025), ]
  keep_eids = bingabc_biosamples[grep("BLD", bingabc_biosamples[,2]), 1]

  tabb = tabb_pre[which(tabb_pre$CellType %in% keep_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene, tabb$Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

EpiMap_BLD_EG_preds <- function(scores,
                        pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap",
                        key_file="/lila/data/deyk/ENCODE/EpiMap/main_metadata_table.tsv",
                        output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/EpiMap_Blood_EG_predictions.txt"){
    library(data.table)
    biosamples_listed = as.character(sapply(as.character(sapply(list.files(pred_cell, pattern = "linking"),
                                                                function(x) return(strsplit(x, "_")[[1]][3]))),
                                            function(xx) return(strsplit(xx, "[.]")[[1]][1])))
    epimap_biosamples = read.delim(paste0(key_file))
    keep_eids = epimap_biosamples$id[c(grep("Blood & T-cell", epimap_biosamples$GROUP),
                                       grep("HSC & B-cell", epimap_biosamples$GROUP),
                                       grep("Lymphoblastoid", epimap_biosamples$GROUP))]

    keep_eids = intersect(keep_eids, biosamples_listed)

    pooled_list = list()
    for(numl in 1:length(keep_eids)){
      preds_tabb = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap/",
                                           "linking_collated_", keep_eids[numl],".bed.gz")))
      preds_tabb = preds_tabb[which(preds_tabb$Score > 0.0045), ]
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_list[[numl]] = tabb2
      cat("Read file", numl, "out of ", length(keep_eids), " files \n")
    }
    pooled_tabb = do.call(rbind, pooled_list)
    fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


Baseline_BLD_EG_preds <- function(
                             pred_type = "dist_to_tss",
                             preds_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Baseline_EG/preds.txt",
                             output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_dist_to_tss_Blood_EG_predictions.txt")
{
  library(data.table)
  preds_biosample = read.table(preds_file, header = F)[,1]
  keep_eids = c("CD14-positive_monocyte_ID_370", "CD4-positive_alpha-beta_T_cell_ID_120", "CD8-positive_alpha-beta_T_cell_ID_148",
                "GM12878_ID_2724", "GM23338_ID_2088", "K562_ID_2644", "natural_killer_cell_ID_2608", "T-cell_ID_61")

  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Baseline_EG/",
                          keep_eids[numl], "/", pred_type, ".tsv.gz"))){

      preds_tabb = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Baseline_EG/",
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

      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end,
                               preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }
  }
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


Random_calc_preds <- function(scores,
                              pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/ABC_2022",
                              output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Random_Blood_EG_predictions.txt"){
  library(data.table)
  ll = list.files(pred_cell)
  biosamples_file = read.delim(paste0(pred_cell, "/", "DNase_H3K27ac.default_metadata_full.101322.qc.tsv"))
  biosamples = biosamples_file$Biosample[which(biosamples_file$default == "True")]

  keep_eids = biosamples[unique(unlist(lapply(c("T-cell", "T_cell", "T-helper", "CD4", "CD8",  "K562","killer",
                                                "hemato", "GM", "B_cell", "B-cell", "CD14", "CD1c", "macrophage"),
                                              function(x) return(grep(x, biosamples)))))]

  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    if(file.exists(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))){
      if (file.info(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz"))$size != 0){
        preds_tabb = data.frame(fread(paste0(pred_cell, "/preds", "/", keep_eids[numl], "/", "EnhancerPredictionsAllPutative.txt.gz")))
        tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$ABC.Score)
        colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
        pooled_tabb_list[[numl]] = tabb2
        cat("We are at biosample:", numl, "\n")
      }
    }
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  colnames(pooled_tabb) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}



###########################################  Output data files   ####################################################


features = c("H3K27ac_reads_by_dist_to_tss_norm",
             "dist_to_gene",
             "dist_to_tss",
             "nearest_gene",
             "nearest_tss",
             "reads_by_dist_to_tss_norm",
             "within_100kb_of_tss")

for(numf in 1:length(features)){
  out1 = Baseline_BLD_EG_preds (pred_type = features[numf],
                                output_file = paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_",
                                features[numf], "_Blood_EG_predictions.txt"))
  cat("We are at feature:", features[numf])
}

out2 = Ramil_BLD_EG_preds()
out3 = BingABC_BLD_EG_preds()
out4 = EpiMap_BLD_EG_preds()

out5 = ABC_DNaseonly_default_calc_preds()
out6 = ABC_DNaseH3K27ac_calc_preds()
out7 = ABC_DNaseH3K27ac_default_calc_preds()
out8 = EPIraction_default_calc_preds()

out9 =  delta_gencode_calc_preds()
out10 = delta_refseq_calc_preds()
out11 = hiccups_gencode_calc_preds()
out12 = hiccups_refseq_calc_preds()


features = c("H3K27ac_reads_by_dist_to_tss_norm",
             "H3K27ac_reads_by_dist_to_tss",
             "DHS_reads_by_dist_to_tss",
             "DHS_reads_by_dist_to_tss_norm",
             "dist_to_gene",
             "dist_to_tss",
             "nearest_expressed_gene",
             "nearest_expressed_tss")

for(numf in 1:length(features)){
  out1 = Baseline_BLD_EG_preds (pred_type = features[numf],
                                output_file = paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_",
                                                     features[numf], "_Blood_EG_predictions.txt"))
  cat("We are at feature:", features[numf])
}




EPIraction_BLD_EG_preds <- function(
    pred_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Ramil_predictions.complete.thresh0_009.tsv.gz",
    key_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Ramil_biosamples_key.csv",
    output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/EPIraction_Blood_EG_predictions.txt")
{
  library(data.table)
  ramil_biosamples = read.csv(paste0(key_file))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  keep_eids = ramil_biosamples[grep("BLD", ramil_biosamples[,2]), 1]

  tabb = tabb_pre[which(tabb_pre$CellType %in% keep_eids == T), ]
  tabb2 = tabb[, c("chr", "start", "end", "TargetGene", "Score")]
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}




ll = list.files("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions")
oo = cbind.data.frame(ll, "BLD", paste0(ll, ".prec_recall"))
write.table(oo, file = "/n/groups/price/kushal/ENCODE/code/Figure4/Primary_onlyEG/finemap_precision_recall.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)

pred_cell = "/n/groups/price/kushal/ENCODE/data/EG_predictions/GraphReg"
tabb_pre1 = data.frame(fread(paste0(pred_cell, "/", "EG_preds_", "K562", "_hg38_GraphReg_FDR_1_L_8_full_predictions.tsv.gz")))
tabb_pre1 = tabb_pre1[which(tabb_pre1$GraphReg.Score > 0.001368199), ]
tabb_pre2 = data.frame(fread(paste0(pred_cell, "/", "EG_preds_", "GM12878", "_hg38_GraphReg_FDR_1_L_8_full_predictions.tsv.gz")))
tabb_pre2 = tabb_pre2[which(tabb_pre2$GraphReg.Score > 0.001368199), ]

tabb2 = rbind.data.frame(tabb_pre1, tabb_pre2)
tabb2 = tabb2[, c("chr", "start", "end", "TargetGene", "GraphReg.Score")]
colnames(tabb2) = c("chr", "start", "end", "TargetGene", "GraphReg.Score")

fwrite(tabb2, "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/GraphReg2023_Blood_EG_predictions.txt",
       sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)







