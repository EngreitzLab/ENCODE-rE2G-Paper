## Make plots correlating ENCODE-rE2G prediction metrics with gene expression levels from Gasperini
## et al., 2019

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cowplot)
  library(ggpubr)
})

# load gene expression levels and set new column names
cpms <- fread(snakemake@input$expr, col.names = c("TargetGeneEnsemblID", "Gasperini_CPM"))

# load ENCODE-rE2G predictions
e2g <- fread(snakemake@input$pred)

# define gene expression metrics
metric_key <- c(n_E2G_links = "# predicted enhancers",
                max_E2G_score = "Max E2G score for predicted enhancers",
                sum_E2G_score = "Sum (E2G scores of predicted enhancers)",
                sum_ABC_score = "Sum (ABC scores of predicted enhancers)",
                normalizedDNase_prom.Feature = "DNase signal at promoter")

# compute ENCODE-rE2G metrics per gene
e2g_metrics <- e2g %>% 
  filter(isSelfPromoter == FALSE) %>% 
  group_by(TargetGene, TargetGeneEnsemblID, normalizedDNase_prom.Feature) %>% 
  summarize(n_E2G_links = n(),
            max_E2G_score = max(Score),
            sum_E2G_score = sum(Score),
            sum_ABC_score = sum(ABC.Score.Feature),
            .groups = "drop") 

# add gene expression levels
e2g_metrics <- left_join(e2g_metrics, cpms, by = "TargetGeneEnsemblID")

message("Genes in Gasperini file: ", nrow(cpms))
message("Genes before filtering: ", nrow(e2g_metrics))
message("Non-NA genes: ", nrow(e2g_metrics %>% filter(!is.na(Gasperini_CPM))))
    
# convert ENCODE-rE2G metrics to long format
e2g_long <- e2g_metrics %>% 
  pivot_longer(cols = all_of(names(metric_key)), names_to = "metric", values_to = "value") %>% 
  mutate(metric = metric_key[metric]) %>% 
  filter(!is.na(Gasperini_CPM))

# function to create plot for correlating one ENCODE-rE2G metric with gene expression
plot_e2g_metric_vs_expr <- function(metric, e2g_metrics, metric_key) {
  
  # define plotting parameters based on metric to plot
  if (metric == "n_E2G_links") {
    x_min <- 0; x_max <- 50
  } else if (metric == "sum_E2G_score") {
    x_min <- 0; x_max <- 16
  } else if (metric == "normalizedDNase_prom.Feature") {
    x_min <- 0.3; x_max <- NA
  } else if (metric == "max_E2G_score") {
    x_min <- 0.2; x_max <- NA
  } else {
    x_min <- 0; x_max <- NA
  }
  
  # make plot for given metric
  ggplot(e2g_metrics, aes(x = !!sym(metric), y = Gasperini_CPM)) +
    geom_point(color = "#730c0d", shape = 16, size = 1, alpha = 0.2) +
    labs(x = metric_key[metric], y = "TPM of target gene") +
    ylim(c(0, 100)) + xlim(c(x_min, x_max)) +
    scale_y_log10() + 
    stat_cor(aes(label = after_stat(paste0("~rho==", signif(..r.., 3)))),
             method = "spearman", size = 3, label.x.npc = "left", label.y.npc = "top", vjust = 1) +
    theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(size = 7, color = "#000000"),
                            axis.title = element_text(size = 9, color = "#000000"),
                            axis.ticks = element_line(color = "#000000"), legend.position = "none") 
  
}

# make plots correlating each ENCODE-rE2G metric with gene expression levels
plots <- lapply(names(metric_key), FUN = plot_e2g_metric_vs_expr, e2g_metrics = e2g_metrics,
                metric_key = metric_key)

# arrange panels in one figure and save to output file
plot_grid(plotlist = plots, nrow = 2, align = "vertical")
ggsave2(filename = snakemake@output[[1]], height = 8, width = 10)
