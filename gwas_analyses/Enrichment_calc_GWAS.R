library(data.table)
library(R.utils)
library(optparse)
library(GenomicRanges)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--biosample", type="character", default = "data/traits.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

biosample = opt$biosample
# biosample = "encode_e2g_predictions_8988T_ENCSR000EID_thresholded_predictions.tsv.gz"

snps_prom = read.table("/data/deyk/extras/Promoter_Coding_variants.txt")[,1]

biosamples = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG_records_top2/")
pops_files = list.files("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG_records_top2")
traitnames2 = as.character(sapply(pops_files, function(y) return(strsplit(y, "_PoPS")[[1]][1])))
traitnames1 = list.files("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets")

common_traits = intersect(traitnames1, traitnames2)

pooled_tabb1 = c()
for(numchr in 1:22){
  tabb1 = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/ENCODE_E2G_1KG_records_top2/",
                                  biosample, "/", "ENCODE_E2G_1KG_records.", numchr, ".txt")))
  pooled_tabb1 = rbind(pooled_tabb1, tabb1)
  cat("We are at chr:", numchr, "\n")
}

enr_traits_vec = c()
for(num_trait in 1:length(common_traits)){
  pooled_tabb2 = c()
  tabb2 = data.frame(fread(paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/GWAS_Variant_PoPS_1KG_records_top2/",
                                  common_traits[num_trait], "_PoPS_1KG_records..txt")))
  pooled_tabb2 = rbind(pooled_tabb2, tabb2)
  baseline_snps = setdiff(pooled_tabb2$SNP, snps_prom)

  xx1 = pooled_tabb1[match(intersect(baseline_snps, pooled_tabb1$SNP), pooled_tabb1$SNP), 2:3]
  xx2 = pooled_tabb2[match(intersect(baseline_snps, pooled_tabb1$SNP), pooled_tabb2$SNP), 2:3]
  pooled_xx = cbind.data.frame(xx1, xx2)
  match_xx = apply(pooled_xx, 1, function(yy) return(max(table(as.character(yy[!is.na(yy)])))))

  prop_snps_all = length(which(match_xx > 1))/9168131

  finemap_tabb = data.frame(fread(paste0("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets/",
                                         common_traits[num_trait], "/", "variant.list.txt")))
  snps_finemapped = intersect(finemap_tabb$rsid[which(finemap_tabb$pip > 0.1)], baseline_snps)
  if(length(snps_finemapped) < 25 | length(intersect(snps_finemapped, pooled_tabb1$SNP)) <= 10){
    enr_value = 1
  }else{
    bxx1 = pooled_tabb1[match(intersect(snps_finemapped, pooled_tabb1$SNP), pooled_tabb1$SNP), 2:3]
    bxx2 = pooled_tabb2[match(intersect(snps_finemapped, pooled_tabb1$SNP), pooled_tabb2$SNP), 2:3]
    pooled_bxx = cbind.data.frame(bxx1, bxx2)
    match_bxx = apply(pooled_bxx, 1, function(yy) return(max(table(as.character(yy[!is.na(yy)])))))

    prop_snps_gwas_celltype = length(which(match_bxx > 1))/length(snps_finemapped)

    enr_value = prop_snps_gwas_celltype/prop_snps_all
    enr_value[enr_value < 1] = 1
  }
  enr_traits_vec = c(enr_traits_vec, enr_value)
  cat("We are at trait:", num_trait)
}

outdff = data.frame("Trait" = common_traits,
                    "Enrichment" = enr_traits_vec)

write.table(outdff, file = paste0("/data/deyk/ENCODE/ENCODE_E2G_2024/processed/Enrichment_ENCODE_E2G_PoPS_top2_trial2/",
                                  biosample, "_", "ENR.txt"),
            row.names = F, col.names = T, sep = "\t", quote=F)

