
# required functions for bootstrapping
source("/oak/stanford/groups/engreitz/Users/kmualim//bootstrap_performance_v0.3/bootstrap_performance_functions.R")

# other functions are needed for dealing with the output of the CRISPR benchmarking pipeline
source("/oak/stanford/groups/engreitz/Users/kmualim/bootstrap_performance_v0.3/crisprComparisonPlotFunctions.R")

# Prepare input data -------------------------------------------------------------------------------

# load merged data from CRISPR benchmarking pipeline
merged <- fread("/oak/stanford/groups/engreitz/Users/kmualim/CRISPR_comparison/results/example/expt_pred_merged_annot.txt")

# load pred_config file that was used in benchmarking
pred_config <- fread("/oak/stanford/groups/engreitz/Users/kmualim/CRISPR_comparison/resources/example/pred_config_dhs.txt",
                     colClasses = c("alpha" = "numeric", "color" = "character"))

# process pred_config file like CRISPR benchmarking pipeline does for analyses
pred_config <- processPredConfig(pred_config, merged = merged)

# process merged data for benchmarking analyses, including filtering for ValidConnection == TRUE
merged <- processMergedData(merged, pred_config = pred_config,
                            filter_valid_connections = TRUE, include_missing_predictions = TRUE,
                            distToTSS_as_kb = FALSE)

# convert merged data for bootstrapping
merged_bs <- convertMergedForBootstrap(merged, pred_config = pred_config)

# Run bootstraps to get confidence intervals -------------------------------------------------------

# bootstrap AUPRC and compute confidence intervals
ci_auprc <- bootstrapPerformanceIntervals(merged_bs, metric = "auprc", R = 10000, conf = 0.95,
                                          ci_type = "perc", ncpus = 2)

delta_auprc <- bootstrapDeltaPerformance(merged_bs, metric = "auprc", R = 10000, conf = 0.95,
                                         ci_type = "perc", ncpus = 4)

# Make AUPRC barplot -------------------------------------------------------------------------------

# add pretty predictor names
ci_auprc <- left_join(ci_auprc, select(pred_config, pred_uid, pred_name_long),
#                      by = c("id" = "pred_uid"))

# get color for each predictor
pred_colors <- deframe(select(pred_config, pred_name_long, color))

# GENERATE FIG5E
pdf("auprc_ci.pdf", width=15, height=5)
 plot AUPRC barplot
ggplot(ci_auprc, aes(y = fct_reorder(pred_name_long, .x = full), x = full, fill = pred_name_long)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(xmin = lower, xmax = upper), position = position_dodge(0.9), width = 0.4,
                color = "black") +
  labs(x = "AUPRC", fill = "Predictor") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = pred_colors) +
  theme_classic() +
  theme(axis.title.y = element_blank())
dev.off()

# OBTAIN COMPARISON VALUES FOR MAIN TEXT
plotBootstrappedIntervals(delta_auprc)
