#' @title DSP NGS Analysis
#'
#' Pipeline for DSP NGS Analysis that creates ROC curves, calculates AUC, runs DE analysis and creates downstream visualizations
#' Requires: rstudioapi to be installed prior to running on Rstudio
#'
#' @param config Configuration file, see config.yml example file
#' @param DSPcounts DSP negative normalized counts spreadsheet with annotations
#'
#' @return Resulting visualizations and DE analysis
#'
#' @example R/DSP_NGS.R
#'
#' @export dsp_de_analysis
#'
#'
#' @section TODO
# #### plot GOI's vs all other genes, one color per GOI vs background
# #### label top 20 DE genes
# #### box plots of top 5 groups of genes by median of de significance (-log10P) + GOI (or custom set) on y vs gene groups on x axison Plot ROC
# Plot log normalized counts for genes as a function of gene_group

# main function
dsp_de_analysis <- function(de_results = NULL,
                            norm_counts,
                            genes,
                            grouping_var,
                            base_level,
                            results_dir,
                            de_dir,
                            control_var,
                            residuals_dir,
                            qq_dir,
                            pval_cutoff,
                            fc_cutoff,
                            samp_notes) {

  #' @section Log Params
  # Write input parameters to log file
  flog.info("\n \n ############   Experimental Parameters   ################## \n", name="DSP_NGS_log")
  flog.info(paste0("Gene(s) of interest: ", gene_group), name="DSP_NGS_log")
  flog.info(paste0("Gene grouping variable: ", grouping_var), name="DSP_NGS_log")
  flog.info(paste0("Base Level: ", base_level), name="DSP_NGS_log")
  flog.info(paste0("Control: ", control_var), name="DSP_NGS_log")
  flog.info(paste0("Pvalue Cutoff: ", pval_cutoff), name="DSP_NGS_log")
  flog.info(paste0("Fold Change Cutoff: ", fc_cutoff), name="DSP_NGS_log")
  flog.info(paste0("Gene Set Colors: ", color_genes), name="DSP_NGS_log")
  flog.info(paste0("Sample Colors: ", color_samples), name="DSP_NGS_log")

  # log2 normalize counts and transpose dataframe
  rownames(norm_counts) <- norm_counts$Gene
  tdf <- within(norm_counts, rm("Gene"))
  trans_df <- data.frame(t(tdf))
  gene_log2 <- data.frame(log2(trans_df))
  gene_log2$Sample_ID <- rownames(trans_df)

  #' @section Parsing Dataframe
  # Subset samp_notes for gene_grouping variable
  sub_list <-c("Sample_ID", eval(grouping_var), eval(control_var))
  samps <- samp_notes[,sub_list]

  # add gene_grouping labels
  goi_df <- base::merge(gene_log2, samps, by="Sample_ID")

  # convert labels to binary classifier
  class_label <- ifelse(goi_df[grouping_var] == eval(base_level), 1, 0)
  if (sum(class_label) == 0){
    print("Check that config base level variable matches a group in ROI/AOI annotations sheet")
    flog.error("Check that config base level variable matches a group in ROI/AOI annotations sheet", name="DSP_NGS_log")
  }
  colnames(class_label) <- "level"
  goi_df <- cbind(goi_df, class_label) # add labels column to goi_df
  #goi_df$level <- factor(goi_df$level, levels=c(0,1), ordered="True")

  #' @section ROC curve
  #' @return PNG file of ROC curve for top 5 genes with highest coeff
  paste("Plot ROC Curve...")
  # plot ROC curve for top5 genes
  ROC_plot(goi_df, grouping_var, "poisson", gene_group, de_dir, genes)


  #' @section DE Analysis
  # Run LME4 linear mixed model on log negative normalized counts for genes of interest
  #' @return df of FC and Pvals above threshold specified in config file,
  # as well as residual and qqplots for goi

  # reset goi_df and let them loop through levels
  goi_df <- base::merge(gene_log2, samps, by="Sample_ID")

  # run_de_file <- paste(getwd(), "de_analysis.R", sep="/") # Used for testing

  if(is.null(de_results)) {
    de_res <- run_DE(df = goi_df,
                     grouping_var = grouping_var,
                     control_var = control_var,
                     base_level = base_level,
                     residuals_dir = residuals_dir,
                     qq_dir = qq_dir,
                     pval_cutoff = pval_cutoff)
    results_file = paste(de_dir, "de_results.csv", sep="/")
    write.csv(de_res, results_file)
  } else if(is.character(de_results)) {
    de_res <- read.delim(de_results, sep = ',', header = TRUE, as.is = TRUE)
  } else if(is.data.frame(de_results) & all(c('FC','Pval','gene') %in% colnames(de_results))) {
    de_res <- de_results
  } else {
    stop('There is an issue with provided data frame, please cehck de_results variable\n')
  }
  #### uncomment for testing ###
  # saveRDS(de_res, file="de_res.R")
  # de_results <- readRDS(file = "de_res.R")

  # Generate volcano plot
  # return PNG and Plotly versions of Volcano Plot

  tests <- unique(de_res$test)
  for (test in tests) {
    de_vol <- de_res[de_res$test == test, ]
    paste("Generating volcano plot...")
    volcano_file <- paste0(de_dir, "/", test, "_volcano.", fileType)
    plot_volcano(de = de_vol,
                 negative_label = strsplit(test, split=" vs ")[[1]][1],
                 positive_label = strsplit(test, split=" vs ")[[1]][2],
                 target_group = gene_group,
                 targets = genes,
                 fav_targets = fav_genes,
                 outfile = volcano_file,
                 point_color = default_color,
                 n_targets = n_top,
                 plt_title = test,
                 color_list = colors$gene_groups)
  }
  return(de_res)
 }
