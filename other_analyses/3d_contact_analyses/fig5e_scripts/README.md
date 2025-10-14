In order to generate Figure 5e, we first ran the CRISPR Benchmarking pipeline located at this github repo: https://github.com/EngreitzLab/CRISPR_comparison

The configuration files to run the pipeline are located in this folder: 
- config_dhs.yml
- pred_config_dhs.txt 

OUTPUT of the CRISPR Comparison code is a merged file labelled: expt_pred_merged_annot.txt

This was later fed into: plot_auprc_example.R located in this folder 
- plot_auprc_example.R relies on other scripts that have been provided in this folder 
- the output of plot_auprc_example.R generates the AUPRC barplot Fig 5e & the numbers used for comparing different measures of contact listed in the MAIN TEXT
