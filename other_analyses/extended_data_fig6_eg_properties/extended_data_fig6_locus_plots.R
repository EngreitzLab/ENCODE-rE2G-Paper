library(Gviz)
library(GenomicInteractions)
library(rtracklayer)
library(dplyr)
library(data.table)

setwd("/Users/ejagoda/Documents/HGRM/EP_benchmarking/ABC/revision/Fig3/")

#Fig 3b
chr = "chr13"
start = 32340000
end =  33914307
dhs_bed_path = "https://www.encodeproject.org/files/ENCFF922AGQ/@@download/ENCFF922AGQ.bed.gz"
dnase_bigwig_path = "https://www.encodeproject.org/files/ENCFF719PPB/@@download/ENCFF719PPB.bigWig"
dnase_bigwig_file = "ENCFF719PPB.bigWig"
e2g_file_path = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_activated_CD4-positive_alpha-beta_T_cell_ENCDO923WRX_ENCFF968IHM_DNaseOnly_thresholded_predictions.bedpe.gz"
gene = "N4BP2L2"
enh = "no"
cell_type = "take2_activated_CD4-positive,_alpha-beta_T_cell"
fig_pannel = "3b"

#Fig3c
chr = "chr12"
start = 12700000
end =  13200000
#dhs_bed_path = "https://www.encodeproject.org/files/ENCFF821YBP/@@download/ENCFF821YBP.bed.gz"
dnase_bigwig_path = "https://www.encodeproject.org/files/ENCFF019DLM/@@download/ENCFF019DLM.bigWig"
dnase_bigwig_file = "ENCFF019DLM.bigWig"
e2g_file_path = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_pulmonary_artery_endothelial_cell_ENCDO332AAA_ENCFF186ZTP_DNaseOnly_thresholded_predictions.bedpe.gz"
gene = "no"
enh = "chr12:13097540-13099593"
cell_type = "pulmonary_artery_endothelial_cell"
fig_pannel = "3c"

#Fig3d
chr = "chr6"
start = 44210000
end =  44320000
dhs_bed_path = "https://www.encodeproject.org/files/ENCFF616WCY/@@download/ENCFF616WCY.bed.gz"
dnase_bigwig_path = "https://www.encodeproject.org/files/ENCFF153OJC/@@download/ENCFF153OJC.bigWig"
dnase_bigwig_file = "ENCFF153OJC.bigWig"
e2g_file_path = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_Right_ventricle_myocardium_inferior_ENCDO520EJG_ENCFF353FPQ_DNaseOnly_thresholded_predictions.bedpe.gz"
gene = "AARS2"
enh = "no"
cell_type = "Right_ventricle_myocardium_inferior"
fig_pannel = "3d"

#Fig3e 
#multi cell
fig_pannel = "3e"
chr = "chr1"
start = 28832906
end =  28938750
gene = "EPB41"
enh = "no"
n_cell_types = 3

#dhs_bed_path1 = "https://www.encodeproject.org/files/ENCFF616WCY/@@download/ENCFF616WCY.bed.gz"
dnase_bigwig_path1 = "https://www.encodeproject.org/files/ENCFF284ARV/@@download/ENCFF284ARV.bigWig"
dnase_bigwig_file1 = "ENCFF284ARV.bigWig"
e2g_file_path1 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_urinary_bladder_ENCDO647EYG_ENCFF351RLC_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type1 = "urinary_bladder"

dnase_bigwig_path2 = "https://www.encodeproject.org/files/ENCFF366YAV/@@download/ENCFF366YAV.bigWig"
dnase_bigwig_file2 = "ENCFF366YAV.bigWig"
e2g_file_path2 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_testis_ENCDO836UAR_ENCFF882ETO_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type2 = "testis"

dnase_bigwig_path3 = "https://www.encodeproject.org/files/ENCFF960LEJ/@@download/ENCFF960LEJ.bigWig"
dnase_bigwig_file3 = "ENCFF960LEJ.bigWig"
e2g_file_path3 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_brain_microvascular_endothelial_cell_ENCDO227AAA_ENCFF203RLS_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type3 = "brain_microvascular_endothelial"



#Fig3f
#multi cell
fig_pannel = "3f"
chr = "chr1"
start = 15337428
end =  15760806
gene = "DNAJC16"
enh = "no"
n_cell_types = 3

#dhs_bed_path1 = "https://www.encodeproject.org/files/ENCFF616WCY/@@download/ENCFF616WCY.bed.gz"
dnase_bigwig_path1 = "https://www.encodeproject.org/files/ENCFF284ARV/@@download/ENCFF284ARV.bigWig"
dnase_bigwig_file1 = "ENCFF284ARV.bigWig"
e2g_file_path1 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_urinary_bladder_ENCDO647EYG_ENCFF351RLC_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type1 = "urinary_bladder"

dnase_bigwig_path2 = "https://www.encodeproject.org/files/ENCFF366YAV/@@download/ENCFF366YAV.bigWig"
dnase_bigwig_file2 = "ENCFF366YAV.bigWig"
e2g_file_path2 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_testis_ENCDO836UAR_ENCFF882ETO_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type2 = "testis"

dnase_bigwig_path3 = "https://www.encodeproject.org/files/ENCFF960LEJ/@@download/ENCFF960LEJ.bigWig"
dnase_bigwig_file3 = "ENCFF960LEJ.bigWig"
e2g_file_path3 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_brain_microvascular_endothelial_cell_ENCDO227AAA_ENCFF203RLS_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type3 = "brain_microvascular_endothelial"

#okay for many
#Fig3f
#multi cell
fig_pannel = "3f_many"
chr = "chrX"
start = 277280
end =  739464
gene = "SHOX"
enh = "no"
n_cell_types = 3

#dhs_bed_path1 = "https://www.encodeproject.org/files/ENCFF616WCY/@@download/ENCFF616WCY.bed.gz"
dnase_bigwig_path1 = "https://www.encodeproject.org/files/ENCFF368ZUB/@@download/ENCFF368ZUB.bigWig"
dnase_bigwig_file1 = "ENCFF368ZUB.bigWig"
e2g_file_path1 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_adipocyte_ENCDO029QEV_ENCFF239CCZ_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type1 = "adipocyte"

dnase_bigwig_path2 = "https://www.encodeproject.org/files/ENCFF515BQL/@@download/ENCFF515BQL.bigWig"
dnase_bigwig_file2 = "ENCFF515BQL.bigWig"
e2g_file_path2 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_fibroblast_of_skin_of_left_biceps_ENCDO576LXR_ENCFF720NZB_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type2 = "fibroblast_of_skin_of_left_biceps"

dnase_bigwig_path3 = "https://www.encodeproject.org/files/ENCFF418OBI/@@download/ENCFF418OBI.bigWig"
dnase_bigwig_file3 = "ENCFF418OBI.bigWig"
e2g_file_path3 = "https://mitra.stanford.edu/engreitz/oak/Projects/Benchmarking/Predictors/ENCODE-E2G/thresholded_bedpe/encode_e2g_predictions_osteoblast_ENCDO271AAA_ENCFF977LFT_DNaseOnly_thresholded_predictions.bedpe.gz"
cell_type3 = "osteoblast"

#SHOX
#encode_e2g_predictions_adipocyte_ENCDO029QEV_ENCFF239CCZ_DNaseOnly_thresholded_predictions.bedpe.gz
#encode_e2g_predictions_fibroblast_of_skin_of_left_biceps_ENCDO576LXR_ENCFF720NZB_DNaseOnly_thresholded_predictions.bedpe.gz
#encode_e2g_predictions_osteoblast_ENCDO271AAA_ENCFF977LFT_DNaseOnly_thresholded_predictions.bedpe.gz

if (n_cell_types != 1){
  name = paste0("locus_plot_",n_cell_types,"cell_types_",gene,"_",chr,"_",start,"_",end)
  outdir = "locus_plots_from_R"
}



#download stuff manually - this doesn't work right
for (i in 1:n_cell_types){
  download.file(url = get(paste0("dnase_bigwig_path",i)), destfile = get(paste0("dnase_bigwig_file",i)))
}



#dhs <- fread(dhs_bed_path)


# import gtf file containing GENCODE v32 genome annotations as a GRanges object

annot <- import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
                format = "gtf")

#annot = import("gencode.v29.annotation.gtf.gz",format = "gtf")

# extract annotations on protein-coding genes and lincRNAs exons
genes <- annot[annot$type == "exon" &
                 annot$gene_type %in% c("protein_coding", "lincRNA") & 
                 annot$transcript_type %in% c("protein_coding", "lincRNA")]

genes <- genes %>% 
  as.data.frame() %>%
  select(chromosome = seqnames, start, end, width, strand, feature = gene_type,
         gene = gene_id, exon = exon_id, transcript = transcript_id,
         symbol = gene_name)

genes_track <- GeneRegionTrack(genes, genome = "hg38", name = "Genes",
                               transcriptAnnotation = "symbol",
                               collapseTranscripts = "meta",
                               col = "#666666", fill = "#666666",)

#genes_track <- GeneRegionTrack(genes, genome = "hg38", name = "Genes",
#                              transcriptAnnotation = "symbol",
#                              collapseTranscripts = "FALSE",
#                              col = "#666666", fill = "#666666",)
                               

# create axis track
axis_track <- GenomeAxisTrack()

## Add DNase-seq hypersensitive sites (DHS)


# create GRanges object
#dhs <- makeGRangesFromDataFrame(dhs, seqnames.field = "V1", start.field = "V2",
#                                end.field = "V3", starts.in.df.are.0based = TRUE)

# create annotation tracks containing DHS in our locus
locus <- GRanges(seqnames = chr, ranges = IRanges(start, end))
#dhs_track <- AnnotationTrack(subsetByOverlaps(dhs, locus), genome = "hg38",
#                             name = "DHS", fill = "steelblue", col = "steelblue",stacking = "pack")

# create DataTrack from DNase-seq bigWig file

for (i in 1:n_cell_types){
  assign(paste0("dnase_bigwig",i),get(paste0("dnase_bigwig_file",i)))
  if (i == 1){
    dnase_track1 = DataTrack(range = get(paste0("dnase_bigwig",i)), genome = "hg38", name = "DNase-seq",
                             type = "polygon", fill.mountain = rep("steelblue", 2),
                             col.mountain = "steelblue")
  }else if (i == 2){
    dnase_track2 = DataTrack(range = get(paste0("dnase_bigwig",i)), genome = "hg38", name = "DNase-seq",
                             type = "polygon", fill.mountain = rep("steelblue", 2),
                             col.mountain = "steelblue")
  } else if (i == 3){
    dnase_track3 = DataTrack(range = get(paste0("dnase_bigwig",i)), genome = "hg38", name = "DNase-seq",
                             type = "polygon", fill.mountain = rep("steelblue", 2),
                             col.mountain = "steelblue")
  }
}




## Adding regulatory interactions
#We will start by loading ENCODE-rE2G predicted regulatory interactions from a .bedpe
#file. The interactions in the file are stored in a `Pairs` object:
  
for (i in 1:n_cell_types){
  e2g_file <- get(paste0("e2g_file_path",i))
  e2g =  import(e2g_file, format = "bedpe")
  mcols(e2g)$distTSS = abs(as.numeric(start(e2g@first))-as.numeric(start(e2g@second)))
  # get gene names from interaction names
  genes <- sub("\\|.+", "", mcols(e2g)[["name"]])
  
  # create GenomicInteractions object
  e2g_int <- GenomicInteractions(anchor1 = S4Vectors::first(e2g),
                                 anchor2 = S4Vectors::second(e2g),
                                 counts = mcols(e2g)[["score"]],
                                 gene = genes,
                                 enh = paste0(e2g@first),
                                 distTSS = mcols(e2g)[["distTSS"]])
  if (gene != "no" & enh == "no"){
    e2g_int_subset <- e2g_int[e2g_int$gene == gene,]
  }else if (gene == "no" & enh != "no"){
    e2g_int_subset <- e2g_int[e2g_int$enh == enh,]
  }else if (gene != "no" & enh != "no"){
    e2g_int_subset <- e2g_int[e2g_int$enh == enh & e2g_int$gene == gene,]
  }else if (gene == "no" & enh == "no"){
    e2g_int_subset ==  e2g_int
  }
  
  mcols(e2g_int_subset)$counts = 0.25
  #assign(paste0("e2g_int_subset",i), e2g_int_subset)
  e2g_track <- InteractionTrack(e2g_int_subset, chromosome = chr, name = "ENCODE-E2G")
  displayPars(e2g_track) <- list(col.interactions = "#DC6464",
                                 col.anchors.line = "#DC6464",
                                 col.anchors.fill = "#DC6464",
                                 col.outside = "#DC6464")
  
  assign(paste0("e2g_track",i),InteractionTrack(e2g_int_subset, chromosome = chr, name = "ENCODE-E2G"))

}


## Adding regulatory interactions

#We can include the ENCODE-rE2G score in the `counts` column of the metadata and also
#include the name of the involved gene for filtering shown later.


#for 3d

#xx = data.frame(table(e2g_int$gene))

## Adding regulatory interactions
#Interactions from predictive models like ENCODE-rE2G can get very busy. Depending on
#the purpose of a plot, it's possible to only show interactions for one gene by
#subsetting the interactions.

# subset interactions for only one gene





# re-create interaction track for subset


if (fig_pannel == "3b"){
  e2g_int_subset = e2g_int_subset[2:3,]
}

if (fig_pannel == "3c"){
  e2g_int_subset = e2g_int_subset[1:5,]
}

#for fig3c
#e2g_int_subset = e2g_int_subset[1:5,]


track_sizes <- c( 0.25, rep(0.5,n_cell_types), rep(0.5,n_cell_types))
# make locus plot showing all tracks



#As you can notice, some of the labels on the plot aren't visible anymore and maybe
#the gray boxes aren't that nice in a figure for publication. We can change the
#display parameters of all tracks further to rotate the labels, adjust the font size
#and color, and remove the gray title background.

# combine all tracks into one list

#genes_track@stacking = "full"
if (n_cell_types == 3){
  all_tracks <- list(genes_track, e2g_track1,dnase_track1,e2g_track2,dnase_track2,e2g_track3,dnase_track3)
}


plotTracks(trackList = all_tracks[2:3], sizes = track_sizes[2:3], chromosome = chr,
           from = start, to = end, background.title = "transparent",
           col.axis = "black", fontcolor.title = "black", fontface.title = 1,
           rotation.title = 0, cex.axis = 1, cex.title = 1, margin = 20)

# make locus plot with custom label settings
pdf_name = paste0(outdir,"/",fig_pannel,"_part1_",name,".pdf")
pdf(pdf_name, height = 5, width = 10)
plotTracks(trackList = all_tracks[1:4], sizes = track_sizes[1:4], chromosome = chr,
           from = start, to = end, background.title = "transparent",
           col.axis = "black", fontcolor.title = "black", fontface.title = 1,
           rotation.title = 0, cex.axis = 1, cex.title = 1, margin = 20)
dev.off()

pdf_name = paste0(outdir,"/",fig_pannel,"_part2_",name,".pdf")
pdf(pdf_name, height = 5, width = 10)
plotTracks(trackList = all_tracks[c(1,5:7)], sizes = track_sizes[c(1,5:7)], chromosome = chr,
           from = start, to = end, background.title = "transparent",
           col.axis = "black", fontcolor.title = "black", fontface.title = 1,
           rotation.title = 0, cex.axis = 1, cex.title = 1, margin = 20)
dev.off()



## Saving the plot to a file

#To export out plot to a .pdf file that can be edited in Adobe Illustrator, we have to
#plot it to an open pdf file connection.

# open pdf file connection
pdf_name = paste0(outdir,"/",fig_pannel,"_",name,".pdf")
pdf(pdf_name, height = 5, width = 10)

# make locus plot with updated display parameters
plotTracks(trackList = all_tracks, sizes = track_sizes, chromosome = chr,
           from = start, to = end, background.title = "transparent",
           col.axis = "black", fontcolor.title = "black", fontface.title = 1,
           rotation.title = 0, cex.axis = 1, cex.title = 1, margin = 20)

# close pdf file
dev.off()

setEPS()
eps_name = paste0(outdir,"/",fig_pannel,"_",name,".eps")
postscript("eps_name")
# make locus plot with updated display parameters
plotTracks(trackList = all_tracks, sizes = track_sizes, chromosome = chr,
           from = start, to = end, background.title = "transparent",
           col.axis = "black", fontcolor.title = "black", fontface.title = 1,
           rotation.title = 0, cex.axis = 1, cex.title = 1, margin = 20)

# close pdf file
dev.off()


