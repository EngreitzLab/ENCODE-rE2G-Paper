
library(arrow)
library(data.table)
dff = read_parquet("/lila/data/deyk/opentargets/l2g_OpenTargets/part-00099-38980086-a2d9-4855-88f8-4936870e64a0-c000.snappy.parquet")

ll = list.files("/lila/data/deyk/opentargets/study_index_OpenTargets/")
merged_dff = c()
for(numl in 1:length(ll)){
  temp_dff = as.data.frame(read_parquet(paste0("/lila/data/deyk/opentargets/study_index_OpenTargets/", ll[numl])))
  xx = cbind.data.frame(temp_dff$study_id, temp_dff$trait_reported)
  merged_dff = rbind.data.frame(merged_dff, xx)
  cat("We are at file:", numl, "\n")
}

colnames(merged_dff) = c("Study.ID", "Trait")



trait_names = c("Attention deficit hyperactivity disorder", "Alzheimer", "Anorexia",
"Autism", "Bipolar disorder", "Celiac", "Crohn's disease", "Depressive symptoms",
"Diverticulitis", "Inflammatory bowel disease", "Lupus", "Major depressive disorder",
"Multiple sclerosis", "Neuroticism", "Rheumatoid arthritis", "Schizophrenia", "Ulcerative colitis",
"Cardiovascular disease", "Atrial fibrillation", "Albumin-globulin ratio", "at menarche",
"at menopause", "Autoimmune disease", "Albumin", "Alkaline phosphatase",
"Alanine aminotransferase", "Alzheimer", "Apolipoprotein A", "Apolipoprotein B",
"Aspartate aminotransferase", "Asthma", "Balding", "Basophil count", "Body fat",
"Blood clot in the lung", "Body mass index", "Breast cancer", "Body", "Calcium",
"Coronary artery disease", "Cholelithiasis", "Colorectal cancer", "C-reactive protein",
"Diastolic blood pressure", "Seen doctor", "Blood clot in the leg",
"Bone mineral density", "Estimated glomerular filtration rate", "Estimated glomerular filtration rate",
"Eosinophil count", "Fed up feelings", "FEV1/FVC", "Fibroblastic disorders", "Gamma glutamyl transferase",
"Glaucoma", "Glucose", "Guilty feelings", "Hemoglobin", "Hemoglobin A1c", "High density lipoprotein cholesterol",
"Height", "Hematocrit", "Hypothyroidism", "Inflammatory bowel disease", "Insulin", "Inguinal hernia",
"Insomnia", "Irritability", "Low density lipoprotein cholesterol", "Lipoprotein A",
"Loneliness", "loss of chromosome", "Lymphocyte count", "Mean arterial pressure",
"Mean corpuscular hemoglobin", "Mean corpuscular hemoglobin concentration", "Multisite chronic pain",
"Mean corpuscular volume", "Migraine", "Miserableness", "Monocyte count", "Mood swing", "Morning person",
"Nervous feelings", "Neuroticism", "Neutrophil count", "Platelet count", "Pulse pressure", "Prostate cancer",
"Red blood cell count", "Risk taking", "Systolic blood pressure", "Sensitivity",
"Sex hormone binding globulin", "Smoking", "Smoking", "Suffer from nerves", "Type 2 diabetes",
"Type 2 diabetes", "Total bilirubin", "Total cholesterol", "Tense", "Testosterone", "Testosterone",
"Triglyceride", "Total protein", "Uric acid", "Urea", "Vitamin D", "White blood cell count",
"Waist-to-hip ratio", "Worrier", "Worry too long after embarrassment")

study_traits = c()
for(numt in 1:length(trait_names)){
  idx1 = unique(unlist(sapply(c(trait_names[numt]), function(x) return(grep(x, merged_dff$Trait)))))
  study_traits = rbind(study_traits, cbind(merged_dff[idx1, ], trait_names[numt]))
}
colnames(study_traits) = c("Study.ID", "Trait full name", "Trait Category")

write.table(study_traits, file = "/data/deyk/opentargets/opentargets_pops/study_traits_pops.txt",
            row.names = F, col.names = T, sep = "\t", quote=F)


study_traits = read.delim("/data/deyk/opentargets/opentargets_pops/study_traits_pops.txt")
ll=list.files("/lila/data/deyk/opentargets/l2g_OpenTargets/")
names = c()
for(numl in 1:length(ll)){
  dff = as.data.frame(read_parquet(paste0("/lila/data/deyk/opentargets/l2g_OpenTargets/", ll[numl])))
  names = c(names, intersect(dff$study_id, study_traits$Study.ID))
  cat("We are at chunk:", numl, "\n")
}

traits2 = unique(names)
study_traits2 = study_traits[which(traits2 %in% study_traits$Study.ID == T), ]

write.table(study_traits2, file = "/data/deyk/opentargets/opentargets_pops/study_traits_pops_filtered_L2G.txt",
            row.names = F, col.names = T, sep = "\t", quote=F)


