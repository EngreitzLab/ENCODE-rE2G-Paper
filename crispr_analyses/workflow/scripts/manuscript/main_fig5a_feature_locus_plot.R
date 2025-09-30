## Make feature locus plot for one gene

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(Gviz)
  library(GenomicInteractions)
})

## Define functions ================================================================================

# main function to create feature locus plot
make_feature_locus_plot <- function(features, links, gene_annot, gene, feature_config, link_col,
                                    upstream = 1e+06, downstream = 1e+06) {
  
  # process input ----------------------------------------------------------------------------------
  
  # get gene TSS coordinates
  tss <- resize(range(gene_annot[gene_annot$symbol == gene]), width = 1, fix = "start")
  
  # only retain features and links for selected gene score_cols
  features <- filter(features, TargetGene == !!gene)
  links <- filter(links, grepl(X7, pattern = paste0("^", gene, "\\|")))
  
  # create tracks ----------------------------------------------------------------------------------
  
  # create links track
  link_name <- ifelse(!is.null(names(link_col)), yes = names(link_col), no = "links")
  links_track <- create_links_track(links, tss = tss, link_col = link_col, name = link_name)
  
  # create predictor score tracks (if any in table)
  score_cols <- pull(filter(feat_config, feature_type == "score"), feat_col)
  if (length(score_cols) > 0) {
    score_tracks <- lapply(score_cols, FUN = create_enh_feature_track, features = features,
                           feat_config = feat_config)
  } else {
    score_tracks <- NULL
  }
  
  # create candidate elements track
  element_track <- AnnotationTrack(features, genome = "hg38",
                                   chromosome = as.character(seqnames(tss)),
                                   name = "Candidate elements",
                                   fill = "black", col = "black")
  
  # create enhancer feature tracks
  enh_feat_cols <- pull(filter(feat_config, feature_type == "enh_feature"), feat_col)
  enh_feature_tracks <- lapply(enh_feat_cols, FUN = create_enh_feature_track, features = features,
                               feat_config = feat_config)
  
  # create promoter feature tracks
  prom_feat_cols <- pull(filter(feat_config, feature_type == "prom_feature"), feat_col)
  prom_feature_tracks <- lapply(prom_feat_cols, FUN = create_prom_feature_track, features = features,
                                feat_config = feat_config)
  
  # make chromosome picture track
  ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = as.character(seqnames(tss)))
  
  # make genome coordinates track
  gtrack <- GenomeAxisTrack(genome = "hg38", chromosome = as.character(seqnames(tss)))
  
  # make genes track
  geneTrack <- GeneRegionTrack(gene_annot, genome = "hg38", name = "Genes",
                               transcriptAnnotation = "symbol", collapseTranscripts = "meta",
                               fill = "gray55", col = "gray55")
  
  # Assemble plot ----------------------------------------------------------------------------------
  
  # all tracks to plot
  all_tracks <- c(ideoTrack, gtrack, geneTrack, links_track, score_tracks, element_track,
                  prom_feature_tracks, enh_feature_tracks)
  
  # set common display parameters across all tracks
  all_tracks <- lapply(all_tracks, FUN = function(x) {
      displayPars(x) <- list(background.title = "transparent", col.axis = "black",
                             fontcolor.title = "black", fontface.title = 1,
                             rotation.title = 45, cex.title = 0.5)
      return(x)
    })
  
  # set track heights
  track_heights <- c(0.5, 0.5, 0.5, 0.5, rep(0.5, length(score_tracks)), 0.15,
                     rep(0.5, length(prom_feature_tracks)),
                     rep(0.5, length(enh_feature_tracks)))
  
  # make locus plot
  plotTracks(trackList = all_tracks, sizes = track_heights,
             chromosome = as.character(seqnames(tss)),
             from = c(max(start(tss) - upstream, 0)),
             to = start(tss) + downstream)
  
}

# function to create interaction track for links
create_links_track <- function(links, tss, link_col, name) {
  enh <- makeGRangesFromDataFrame(links, seqnames.field = "X1", start.field = "X2",
                                  end.field = "X3", keep.extra.columns = TRUE)
  inter <- GenomicInteractions(enh, tss[rep_len(1, length(enh))],  p.value = 1 - enh$X8)
  track <- InteractionTrack(inter, name = name)
  displayPars(track) <- list(col.interactions = link_col, col.anchors.line = link_col,
                             col.anchors.fill = link_col, interaction.measure = "p.value")
  return(track)
}

# function to create data track for a given enhancer feature
create_enh_feature_track <- function(feat_col, features, feat_config) {
  config <- filter(feat_config, feat_col == !!feat_col)
  feature <- makeGRangesFromDataFrame(select(features, chr, start, end, data = all_of(feat_col)),
                                      starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  track <- DataTrack(range = feature, data = "data", genome = "hg38",
                     type = config$plot_type, col = config$color, fill.histogram = config$color,
                     col.histogram = config$color, name = config$feat_name)
  return(track)
}

# function to create data track for a given promoter feature
create_prom_feature_track <- function(feat_col, features, feat_config) {
  config <- filter(feat_config, feat_col == !!feat_col)
  feature <- features %>% 
    mutate(start = TargetGeneTSS - 1000, end = TargetGeneTSS + 1000) %>% 
    select(chr, start, end, data = all_of(feat_col)) %>% 
    distinct() %>% 
    mutate(data = if_else(is.logical(data), true = as.integer(data), false = data)) %>% 
    makeGRangesFromDataFrame(., starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  track <- DataTrack(range = feature, data = "data", genome = "hg38",
                     type = config$plot_type, col = config$color, fill.histogram = config$color,
                     col.histogram = config$color, name = config$feat_name)
  return(track)
}


## Create locus plot for a given gene ==============================================================

# load ENCODE-E2G features
e2g_feat <- read_tsv(snakemake@input$e2g_features, show_col_types = FALSE)
e2g_feat <- dplyr::rename(e2g_feat, chr = `#chr`)

# filter out any promoter elements, except 'self promoters' of target genes
e2g_feat <- filter(e2g_feat, !(class == "promoter" & isSelfPromoter == FALSE))

# correct 'numTSSEnhGene.Feature' and 'numCandidateEnhGene.Feature' to not count target gene TSS or
# the element itself
e2g_feat <- e2g_feat %>% 
  mutate(numTSSEnhGene.Feature = pmax(0, numTSSEnhGene.Feature - 1)) %>% 
  mutate(numCandidateEnhGene.Feature = pmax(0, numCandidateEnhGene.Feature - 1))

# load bedpe file with predicted links
e2g_links <- read_tsv(snakemake@input$e2g_links, col_names = FALSE, show_col_types = FALSE)

# colors for links
link_col <- c("ENCODE-E2G links" = "#DC6464")

# load feature config file
feat_config <- read_tsv(snakemake@input$feat_config, show_col_types = FALSE)

# download and import gtf genome annotations
annot <- import(snakemake@input$genome_annot, format = "gtf")

# extract annotations on protein-coding genes and lincRNAs exons
genes <- annot[annot$type == "exon" &
                 annot$gene_type %in% c("protein_coding", "lincRNA") & 
                 annot$transcript_type %in% c("protein_coding", "lincRNA")]

# create GRanges containing required gene annotation data
genes <- genes %>% 
  as.data.frame() %>%
  select(chromosome = seqnames, start, end, width, strand, feature = gene_type, gene = gene_id,
         exon = exon_id, transcript = transcript_id, symbol = gene_name) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# # Re-compute num*EnhGene features ------------------------------------------------------------------
# 
# # function to compute the number of features between enhancers and TSSs
# compute_num_features <- function(pairs, feature) {
#   
#   # get upper and lower window for each pair
#   pairs <- pairs %>% 
#     mutate(pair_id = paste0(name, ";", TargetGene)) %>% 
#     mutate(lower = if_else(start < TargetGeneTSS, true = start, false = TargetGeneTSS),
#            upper = if_else(start < TargetGeneTSS, true = TargetGeneTSS, false = end))
#   
#   # create GRanges object with window for each pair
#   pair_ranges <- pairs %>% 
#     select(chr, start = lower, end = upper, pair_id) %>% 
#     makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
#   
#   # overlap with features
#   overlaps <- countOverlaps(query = pair_ranges, subject = feature) - 1
#   
#   return(overlaps)
#   
# }
# 
# # selected gene locus
# gene <- "PRKAR2B"
# 
# # extract features for selected gene and sort according to element coordinates
# e2g_feat <- e2g_feat %>% 
#   filter(TargetGene == !!gene) %>% 
#   arrange(start, end)
# 
# # get TSS coordinates for all genes
# tss <- annot[annot$type == "gene" & annot$gene_type %in% c("protein_coding", "lincRNA")]
# tss <- split(tss, f = tss$gene_name)
# tss <- unlist(resize(tss, width = 1, fix = "start"))
# 
# # add gencode TSS coordinates to features
# tss <- tibble(chr = as.character(seqnames(tss)), TargetGeneTSS = start(tss), TargetGene = tss$gene_name)
# e2g_feat <- e2g_feat %>% 
#   select(-TargetGeneTSS) %>% 
#   left_join(tss, by = c("chr", "TargetGene")) %>% 
#   relocate(TargetGeneTSS, .after = "TargetGene") %>% 
#   arrange(start, end)
# 
# # compute number of candidates between enhancers and TSS
# candidates <- distinct(select(e2g_feat, chr, start, end))
# candidates <- makeGRangesFromDataFrame(candidates)
# e2g_feat$numCandidateEnhGene.Feature <- compute_num_features(e2g_feat, feature = candidates)
# 
# # compute the number of TSSs between enhancers and TSS
# tss <- makeGRangesFromDataFrame(tss, start.field = "TargetGeneTSS", end.field = "TargetGeneTSS")
# e2g_feat$numTSSEnhGene.Feature <- compute_num_features(e2g_feat, feature = tss)

# Make locus plot ----------------------------------------------------------------------------------

# selected gene locus
gene <- "PRKAR2B"

# make locus plot
pdf("results/manuscript/plots/main_fig5a_prkar2b_features.pdf", height = 6, width = 10)
make_feature_locus_plot(e2g_feat, links = e2g_links, gene_annot = genes, gene = gene,
                        feature_config = feat_config, link_col = link_col, upstream = 2e+05,
                        downstream = 2e+05)
dev.off()
