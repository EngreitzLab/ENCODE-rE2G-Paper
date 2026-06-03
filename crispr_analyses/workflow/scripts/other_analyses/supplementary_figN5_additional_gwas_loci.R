# Make additional GWAS locus plots for Supplementary Figure N5.2

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicAlignments)
  library(Gviz)
  library(GenomicInteractions)
})

# Define functions ---------------------------------------------------------------------------------

# make locus plot for one variant
make_locus_plot_variant <- function(variant, pred_tracks, axis_track, genes_track, window_width) {
  
  # highlight GWAS variant in prediction tracks
  highlight_width <- window_width * 0.0075
  pred_tracks_hl <- HighlightTrack(trackList = pred_tracks,
                                   start = start(variant) - highlight_width / 2,
                                   width = highlight_width,
                                   chromosome = as.character(seqnames(variant)),
                                   fill = "#E6E6E6", col = NA)
  
  
  # combine with genes and axis tracks and set track heights
  all_tracks <- c(list(axis_track, genes_track), pred_tracks_hl)
  track_sizes <- c(1, 0.5, rep(1, length(pred_tracks_hl)))
  
  # make locus plot for provided variant
  ext <- window_width / 2  # how far to extend from variant to get to specified window width
  plotTracks(trackList = c(all_tracks), sizes = track_sizes,
             chromosome = as.character(seqnames(variant)), from = start(variant) - ext,
             to = end(variant) + ext, background.title = "transparent", col.axis = "black",
             fontcolor.title = "black", fontface.title = 1, rotation.title = 0, cex.axis = 0.6,
             cex.title = 0.6, margin = 20)
  
}


# function to load E-G predictions in .bedpe format and create an interaction track
make_e2g_track <- function(file, name, color = "gray", gene = NULL, scale_score = FALSE,
                           remove_promoters = TRUE, plot.outside = TRUE) {
  
  # load E2G file
  e2g <- import(file, format = "bedpe")
  
  # remove self-promoters if specified
  if (remove_promoters == TRUE) {
    enhancers <- GenomicAlignments::first(e2g)
    promoters <- GenomicAlignments::second(e2g)
    self_promoters <- width(pintersect(enhancers, promoters)) > 0
    e2g <- e2g[!self_promoters]
  }
  
  # create GenomicInteractions object
  genes <- sub("\\|.+", "", mcols(e2g)[["name"]])
  e2g_int <- GenomicInteractions(anchor1 = S4Vectors::first(e2g),
                                 anchor2 = S4Vectors::second(e2g),
                                 counts = mcols(e2g)[["score"]],
                                 gene = genes)
  
  # only extract interactions of a specific gene (or multiple genes)
  if (!is.null(gene)) {
    e2g_int <- e2g_int[e2g_int$gene %in% gene]
  }
  
  # set counts to 1 if scaling is disabled
  if (scale_score == FALSE) {
    e2g_int$counts <- 1    
  }
  
  # create interaction track to show
  e2g_track <- InteractionTrack(e2g_int, name = name)
  
  # set other display parameters
  displayPars(e2g_track) <- list(col.interactions = color, col.anchors.line = color,
                                 col.anchors.fill = color, plot.outside = plot.outside, 
                                 col.outside = color)
  
  return(e2g_track)
  
}

# helper function to download bigwig file if it doesn't exist
download_bigwig <- function(file, dir) {
  downloaded_file <- file.path(dir, basename(file))
  if (file.exists(downloaded_file) == FALSE) {
    message("Downloading bigWig file: ", basename(file))
    dir.create(dirname(downloaded_file), recursive = TRUE, showWarnings = FALSE)
    download.file(url = file, destfile = downloaded_file)
  }
}

# Make genomic tracks ------------------------------------------------------------------------------

# import gtf file containing genome annotations
annot <- import(snakemake@input$genome_annot, format = "gtf")

# extract annotations on protein-coding genes exons
genes <- annot[annot$type == "exon" & annot$gene_type == "protein_coding" & 
                 annot$transcript_type == "protein_coding" & annot$tag == "CCDS"]

# create data frame containing required gene annotation data
genes <- genes %>% 
  as.data.frame() %>%
  select(chromosome = seqnames, start, end, width, strand, feature = gene_type, gene = gene_id,
         exon = exon_id, transcript = transcript_id, symbol = gene_name)

# make genes track
genes_track <- GeneRegionTrack(genes, genome = "hg38", name = "Genes",
                               transcriptAnnotation = "symbol",
                               collapseTranscripts = "meta",
                               col = "#676767", fill = "#676767")

# make genome coordinates track
axis_track <- GenomeAxisTrack()

# load ENCODE-rE2G metadata
e2g_meta <- read_tsv(snakemake@input$e2g_metadata, show_col_types = FALSE)

# get .bedpe files with thresholded ENCODE-rE2G predictions and .bigWig files with DNase-seq signal
# to show in locus plot
e2g_samples <- c("hematopoietic_multipotent_progenitor_cell_ENCSR845CFB", "K562_ENCSR000EOT", 
                 "activated_T-cell_ENCSR466SUZ", "B_cell_ENCSR891VOV",
                 "common_myeloid_progenitor__CD34-positive_ENCSR122VUW",
                 "adipocyte_ENCSR405TXU")
e2g_meta_samples <- filter(e2g_meta, Sample %in% e2g_samples)
e2g_bedpe <- deframe(select(e2g_meta_samples, `Biosample term name`, thresholded_bedpe_path))
dnase_bigwig <- deframe(select(e2g_meta_samples, `Biosample term name`, dnase_bigwig_url))

# create interaction tracks for all ENCODE-rE2G samples to show
e2g_genes <- c("BCL11A", "TFRC", "KIT")
e2g_tracks <- mapply(FUN = make_e2g_track, file = e2g_bedpe, name = names(e2g_bedpe),
                     MoreArgs = list(gene = e2g_genes, plot.outside = snakemake@params$plot_outside))

# download all DNase-seq bigWig files if needed
invisible(lapply(dnase_bigwig, FUN = download_bigwig, dir = snakemake@params$scratch_dir))

# create DataTracks for DNase-seq signals
downloaded_files <- file.path(snakemake@params$scratch_dir, basename(dnase_bigwig))
names(downloaded_files) <- names(dnase_bigwig)
dnase_tracks <- lapply(downloaded_files, FUN = DataTrack, genome = "hg38", type = "polygon",
                       name = "DNase-seq", fill.mountain = rep("steelblue", 2),
                       col.mountain = "steelblue")

# create list of E2G and DNase-seq tracks
cell_type_order <- c("hematopoietic multipotent progenitor cell", "K562", "activated T-cell", 
                     "B cell", "common myeloid progenitor, CD34-positive", "adipocyte")
pred_tracks <- unlist(Map(list, e2g_tracks[cell_type_order], dnase_tracks[cell_type_order]))


# Make locus plots ---------------------------------------------------------------------------------

# coordinates for GWAS variants shown in plots (based on information from NHGRI-EBI GWAS catalog)
rs4927708_coords <- GRanges("chr3", IRanges(196047774))
rs218265_coords <- GRanges("chr4", IRanges(54542832))
rs7599488_coords <- GRanges("chr2", IRanges(60491212))

message("Making plot for rs4927708 locus...")
pdf(snakemake@output$rs4927708, width = 10, height = 3.5)
make_locus_plot_variant(rs4927708_coords, pred_tracks = pred_tracks, axis_track = axis_track,
                        genes_track = genes_track, window_width = 80000)
dev.off()

message("Making plot for rs218265 locus...")
pdf(snakemake@output$rs218265, width = 10, height = 3.5)
make_locus_plot_variant(rs218265_coords, pred_tracks = pred_tracks, axis_track = axis_track,
                        genes_track = genes_track, window_width = 4e+05)
dev.off()

message("Making plot for rs7599488 locus...")
pdf(snakemake@output$rs7599488, width = 10, height = 3.5)
make_locus_plot_variant(rs7599488_coords, pred_tracks = pred_tracks, axis_track = axis_track,
                        genes_track = genes_track, window_width = 4e+05)
dev.off()

message("All done!")
