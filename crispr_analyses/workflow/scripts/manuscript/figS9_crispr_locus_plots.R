## Make locus plots for supplementary figures

# attach required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(GenomicAlignments)
  library(Gviz)
  library(GenomicInteractions)
})

## Define functions ================================================================================

# main function to create a locus plot for one gene
make_locus_plot <- function(gene, downstream, upstream, e2g_file, e2g_ext_file, crispr_pairs, dhs,
                            gene_annot, dnase_bam, h3k27ac_bam, e2g_threshold = 0,
                            e2g_ext_threshold = 0) {
  
  # process input data -----------------------------------------------------------------------------
  
  # get gene tss coordintates
  gene_tss <- resize(range(gene_annot[gene_annot$gene_name == gene]), width = 1, fix = "start")
  
  # set gene locus coordinates (window to show)
  locus <- GRanges(seqnames = seqnames(gene_tss),
                   IRanges(start(gene_tss) - downstream, end(gene_tss) + upstream))
  
  # create data frame containing required data for gene annotation track
  genes <- gene_annot %>% 
    as.data.frame() %>%
    select(chromosome = seqnames, start, end, width, strand, feature = gene_type, gene = gene_id,
           exon = exon_id, transcript = transcript_id, symbol = gene_name)
  
  # get ENCODE-E2G scores for gene of interest
  e2g_scores <- extract_scores_gene(file = e2g_file, gene = gene, score_col = "Score",
                                    threshold = e2g_threshold)
  e2g_ext_scores <- extract_scores_gene(file = e2g_ext_file, gene = gene,
                                        score_col = "Score",
                                        threshold = e2g_ext_threshold)
  
  # extract significant CRISPR interactions for the given gene
  crispr_gene <- crispr_pairs %>% 
    filter(Regulated == TRUE, measuredGeneSymbol == gene) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # extract unique crispr enhancers and create GRanges object
  crispr_enh <- crispr_pairs %>% 
    select(chr, start, end) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame() %>% 
    reduce()
  
  # make Gviz tracks -------------------------------------------------------------------------------
  
  # make genome coordinates track
  gtrack <- GenomeAxisTrack()
  
  # make genes track
  geneTrack <- GeneRegionTrack(genes, genome = "hg38", name = "Genes", fill = "#666666",
                               col = "#666666", transcriptAnnotation = "symbol",
                               collapseTranscripts = "meta",
                               chromosome = as.character(seqnames(locus)))
  
  # create annotation tracks for DHS and CRISPR enhancers
  dhs_track <- AnnotationTrack(subsetByOverlaps(dhs, locus), genome = "hg38", name = "DHS",
                               fill = "black", col = "black")
  crispr_enh_track <- AnnotationTrack(subsetByOverlaps(crispr_enh, locus), genome = "hg38",
                                      name = "CRISPRi tested elements", fill = "black",
                                      col = "black")
  
  # create DNase-seq and H3K27ac tracks
  dnase_track <- AlignmentsTrack(dnase_bam, from = start(locus), to = end(locus),
                                 chromosome = seqnames(locus), genome = "hg38", type = "coverage",
                                 fill = "steelblue", col = "steelblue", name = "DNase-seq")
  h3k27ac_track <- AlignmentsTrack(h3k27ac_bam, from = start(locus), to = end(locus),
                                   chromosome = seqnames(locus), genome = "hg38", type = "coverage",
                                   fill = "goldenrod", col = "goldenrod", name = "H3K27ac")
  
  # create ENCODE-E2G interaction track
  e2g_pairs <- GenomicInteractions(e2g_scores,
                                   gene_tss[rep_len(1, length(e2g_scores))],
                                   p.value = 1 - e2g_scores$score - e2g_threshold)
  e2g_track <- InteractionTrack(e2g_pairs, chromosome = as.character(seqnames(locus)),
                                name = "ENCODE-E2G")
  displayPars(e2g_track) <- list(col.interactions = "#C5204C", col.anchors.line = "#C5204C",
                                 col.anchors.fill = "#C5204C", plot.outside = FALSE,
                                 interaction.measure = "counts",
                                 interaction.dimension = "width")
  
  # create ENCODE-E2G_extended interaction track
  e2g_ext_pairs <- GenomicInteractions(e2g_ext_scores,
                                       gene_tss[rep_len(1, length(e2g_ext_scores))],
                                       p.value = 1 - e2g_ext_scores$score - e2g_ext_threshold)
  e2g_ext_track <- InteractionTrack(e2g_ext_pairs,
                                    chromosome = as.character(seqnames(locus)),
                                    name = "ENCODE-E2G_extended")
  displayPars(e2g_ext_track) <- list(col.interactions = "#83062E", col.anchors.line = "#83062E",
                                     col.anchors.fill = "#83062E", plot.outside = FALSE,
                                     interaction.measure = "counts",
                                     interaction.dimension = "width")
  
  # create CRISPR interaction track
  crispr_pairs <- GenomicInteractions(crispr_gene, gene_tss[rep_len(1, length(crispr_gene))],
                                      p.value = 1-abs(crispr_gene$EffectSize))
  crispr_track <- InteractionTrack(crispr_pairs, chromosome = as.character(seqnames(locus)),
                                   name = "CRISPRi")
  displayPars(crispr_track) <- list(col.interactions = "black", col.anchors.line = "black",
                                    col.anchors.fill = "black", plot.outside = FALSE,
                                    interaction.measure = "counts",
                                    interaction.dimension = "width")
  
  # Create locus plot ------------------------------------------------------------------------------
  
  # combine all tracks into one list
  all_tracks <- c(gtrack, geneTrack, e2g_track, e2g_ext_track, crispr_track, crispr_enh_track, 
                  h3k27ac_track, dnase_track, dhs_track)
  
  # set common display parameters across all tracks
  all_tracks <- lapply(all_tracks, FUN = function(x) {
    displayPars(x) <- list(background.title = "transparent", col.axis = "black",
                           fontcolor.title = "black", fontface.title = 1,
                           rotation.title = 45, cex.title = 0.5)
    return(x)
  })
  
  # relative heights for all tracks
  track_sizes <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.5, 0.5, 0.1)
  
  # plot tracks  
  plotTracks(trackList = all_tracks, chromosome = as.character(seqnames(locus)),
             from = start(locus), to = end(locus), sizes = track_sizes)
  
}

# function to load predictions for one sample, extract scores for one gene
extract_scores_gene <- function(file, gene, score_col, threshold = 0) {
  
  # load predictions and extract scores for desired gene
  pred <- fread(file)
  pred <- filter(pred, TargetGene == gene)
  
  # create GRanges object with scores for all enhancers for desired gene
  scores <- pred %>% 
    filter(class != "promoter") %>% 
    select(chr = `#chr`, start, end, score = all_of(score_col)) %>% 
    filter(score >= threshold)
  
  # create GRanges track from scores data frame
  if (nrow(scores) > 0) {
    scores <- makeGRangesFromDataFrame(scores, keep.extra.columns = TRUE)
  } else {
    scores <- GRanges()
  }
  
  return(scores)
  
}

## Load data =======================================================================================

# download and import gtf genome annotations
annot <- import(snakemake@input$genome_annot, format = "gtf")

# extract annotations on protein-coding genes and linRNAs exons
genes <- annot[annot$type == "exon" &
                 annot$gene_type %in% c("protein_coding", "lincRNA") & 
                 annot$transcript_type %in% c("protein_coding", "lincRNA")]

# filter out problematic hemoglobin transcripts
genes <- genes[!genes$transcript_id %in% c("ENST00000380252.6", "ENST00000292896.3", "ENST00000380237.5")]

# load output from CRISPR benchmarking
crispr_data <- fread(snakemake@input$crispr_data)

# extract tested CRISPR E-G pairs
crispr_pairs <- crispr_data %>% 
  select(chr = chrom, start = chromStart, end = chromEnd, name, measuredGeneSymbol,
         EffectSize, Regulated) %>% 
  distinct()

# load DNase-seq peaks and create GRanges object
dhs <- fread(snakemake@input$dhs_peaks)
dhs <- makeGRangesFromDataFrame(dhs, seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                starts.in.df.are.0based = TRUE)

# bam files containing mapped DNase-seq and H3K27ac reads
dnase_bam <- snakemake@input$dnase_bam
h3k27ac_bam <- snakemake@input$h3k27ac_bam

# files containing ENCODE-E2G predictions
k562_e2g <- snakemake@input$e2g_links
k562_e2g_ext <- snakemake@input$e2g_ext_links

## Make locus plots ================================================================================

message("Making plot for HBE1 locus...")
pdf(snakemake@output$hbe1, height = 6.5, width = 21)
make_locus_plot(gene = "HBE1", downstream = 75000, upstream = 75000, e2g_file = k562_e2g,
                e2g_ext_file = k562_e2g_ext, crispr_pairs = crispr_pairs, dhs = dhs,
                gene_annot = genes, dnase_bam = dnase_bam, h3k27ac_bam = h3k27ac_bam)
dev.off()

# message("Making plot for MYC locus...")
pdf(snakemake@output$myc, height = 6.5, width = 21)
make_locus_plot(gene = "MYC", downstream = 2100000, upstream = 100000, e2g_file = k562_e2g,
                e2g_ext_file = k562_e2g_ext, crispr_pairs = crispr_pairs, dhs = dhs,
                gene_annot = genes, dnase_bam = dnase_bam, h3k27ac_bam = h3k27ac_bam)
dev.off()

message("Making plot for zoomed-in MYC locus...")
pdf(snakemake@output$myc_zoomed, height = 6.5, width = 21)
make_locus_plot(gene = "MYC", downstream = 100000, upstream = 100000, e2g_file = k562_e2g,
                e2g_ext_file = k562_e2g_ext, crispr_pairs = crispr_pairs, dhs = dhs,
                gene_annot = genes, dnase_bam = dnase_bam, h3k27ac_bam = h3k27ac_bam)
dev.off()

message("Making plot for GATA1 locus...")
pdf(snakemake@output$gata1, height = 6.5, width = 21)
make_locus_plot(gene = "GATA1", downstream = 50000, upstream = 50000, e2g_file = k562_e2g,
                e2g_ext_file = k562_e2g_ext, crispr_pairs = crispr_pairs, dhs = dhs,
                gene_annot = genes, dnase_bam = dnase_bam, h3k27ac_bam = h3k27ac_bam)
dev.off()

message("All done!")
