#' @title run_ssGSEA
#'
#' Pathway Analysis helper function to run ssGSEA
#'
#' @param norm_data the normalized log2 transformed data matrix
#' @param annots the annotations data.frame
#' @param pathways are the genesets that we are working with, default is Reactome
#' @param test is the test or contrast that is being compared for pathway analysis, from DE results
#' @param outdir the output directory for plots
#'
#' @return ssGSEA_results The results from ssGSEA analysis
#'
#' @examples
#'  run_ssGSEA(norm_data = norm_data, pathways = pathways, test = test, outdir = outdir)
#'
#' @export run_ssGSEA

run_ssGSEA <- function(norm_data = norm_data, annots, pathways = pathways, test = test, outdir = outdir){
  # Create directory for ssGSEA results
  ssGSEA_dir <- paste(outdir, "ssGSEA", sep = "/")
  dir.create(ssGSEA_dir)

  # Run ssGSEA in GSVA package
  ssGSEA_results <- gsva(norm_data, pathways, method="ssgsea", min.sz = 5, max.sz=500, verbose=FALSE, parallel.sz=1)

  # Write the file
  write.csv(as.data.frame(ssGSEA_results), file=paste0(ssGSEA_dir, "/", test, "_ssGSEA_results.csv"))

  # Format results for DE
  ssGSEA_df <- as.data.frame(ssGSEA_results)
  # dsp_internal auto log2 prior to DE, applying inverse to preserve NES
  ssGSEA_df <- 2 ^ ssGSEA_df
  ssGSEA_df[["Gene"]] <- rownames(ssGSEA_df)
  # recode pathways to numeric prior to running through dsp_de_analysis_parallel
  ssGSEA_df$Gene <- 1:nrow(ssGSEA_df)
  ssGSEA_df_translation <- data.frame(path=row.names(ssGSEA_df), id=1:nrow(ssGSEA_df))
  ssGSEA_df <- ssGSEA_df[, c("Gene", setdiff(colnames(ssGSEA_df), "Gene"))]

  # Run DE in parallel
  de_ssGSEA <- dsp_de_analysis_parallel(de_results = NULL,
                                    norm_counts = ssGSEA_df,
                                    samp_notes = annots,
                                    grouping_var = grouping_var,
                                    base_level = base_level,
                                    control_var = control_var,
                                    de_dir = ssGSEA_dir,
                                    n_processors = numCores)
  # Transfer back the path names
  de_ssGSEA$gene <- ssGSEA_df_translation[match(de_ssGSEA$gene, ssGSEA_df_translation$id), 'path']

  # Reformat and rewrite results to disk
  if(class(de_ssGSEA)=="data.frame"){
    colnames(de_ssGSEA)[colnames(de_ssGSEA) == "gene"] <- "Pathway"
    colnames(de_ssGSEA)[colnames(de_ssGSEA) == "FC"] <- "NES.FC"
    write.csv(de_ssGSEA, file=paste(ssGSEA_dir, "de_results.csv", sep="/"))

    # Generate volcano plot for each test
    for(test in unique(de_ssGSEA$test)) {
      vol_ssGSEA <- de_ssGSEA[de_ssGSEA$test == test, ]
      volcano_file <- paste0(ssGSEA_dir, "/", test, "_ssGSEA_Volcano-AllPathways.", fileType)
      volcano_title <- paste0("ssGSEA Volcano Plot for ", test)
      vol_ssGSEA$Pathway <- str_wrap(vol_ssGSEA$Pathway, 25)
      plot_volcano(de = vol_ssGSEA,
                   negative_label = strsplit(test, split=" vs ")[[1]][1],
                   positive_label = strsplit(test, split=" vs ")[[1]][2],
                   target_group = "pathway",
                   fav_targets = NULL,
                   targets = unlist(vol_ssGSEA$Pathway),
                   outfile = volcano_file,
                   point_color = default_color,
                   top_method = "Significance",
                   n_targets = 10,
                   plt_title = volcano_title,
                   target_ID = "Pathway",
                   FC_ID = "NES.FC",
                   Pval_ID = "Pval",
                   color_list = NULL,
                   save_plot = TRUE)
    }
  } else {
    cat("Statistical analysis on ssGSEA results could not be performed,
          ensure grouping variables are filled out appropriately in the config file.")
  }

  return(ssGSEA_results)
}
