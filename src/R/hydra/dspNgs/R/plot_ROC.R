#' Plot ROC Curve
#'
#' Transposes counts df and creates ROC curve plot with top 5 genes with largest coefficients
#'
#' @param df negative normalized DSP count data transposed to wide format, one column per gene
#' @param grouping_var column from annotations to group by
#' @param family glm distribution, aka poisson
#' @param gene_group gene group of interest
#' @param outdir path to output directory
#' @param genes list of genes of interest
#'
#' @return image file of top 5 genes ROC curve
#' @return dataframe of AUC for each gene
#'
#' @examples
#'  ROC_plot(top5_coeff, gene_group, outdir)
#'
#' @export ROC_plot


# calculate auc table and return
# create ROC plot and save to image
ROC_plot <- function(df, grouping_var, family, gene_group, outdir, genes){

  # drop non-numeric columns
  drop_list <- c("Sample_ID", control_var, grouping_var)
  testing_df <- df[, !colnames(df) %in% drop_list] # features and classes only

  fit_df = NULL

  for (x in colnames(testing_df[, !colnames(testing_df) %in% "level"])){

    form = as.character(paste0("level ~ ", x)) # create formula for pairwise comparison of level and each gene column


    fit <- glm(formula=as.formula(form), data=testing_df, family=family) # fit the model

    # format the output
    out_row = data.frame(coeff=coefficients(fit))
    out_row$gene <- rownames(out_row) # add gene information as rows
    out_row <- out_row[out_row$gene != "(Intercept)",] # remove intercept label row

    # add the row to the fit_df
    fit_df = rbind(fit_df, out_row)

  }

  # sort fit_df by abs(coeff)
  fit_df <- fit_df[order(abs(fit_df$coeff), decreasing=TRUE),]
  top5_all_list <- as.vector((head(fit_df$gene, 5)))
  top5_all <- as.vector(c(head(fit_df$gene, 5), "level")) # pull out top 5 genes by abs(coeff)

  # subset counts for top 5 coeff
  if (tolower(gene_group[1]) == "all"){

    ### plot top 5 of all genes ##

    # subset counts table for top 5's
    top5_all_counts <- df[colnames(df) %in% top5_all]

    #  melt dataframes into long format, one row per gene per count
    top5_all_df <- melt_roc(data=df,
                            d="level",
                            m=top5_all_list)

    # setup plot variables
    all_title ="ROC: Top 5 All Genes"
    # create outfile path from input outdir
    all_outfile = paste0(outdir, "/Top5_all_ROC.", fileType)

    #  plot roc
    p_all <- ggplot(top5_all_df, aes(d = D, m = M, color = name)) +
      geom_roc(show.legend = FALSE) +
      geom_abline(slope = 1, intercept = 0, color = 'gray') +
      geom_roc(n.cuts = 0, lwd = 1) +
      ggtitle(all_title) +
      scale_x_continuous("Specificity") +
      scale_y_continuous("Sensitivity")
    if (is.null(preset_theme)) {
      p_all <- p_all + theme_light(base_size = 14)
    }

    # # pull out AUC for each gene
    all_auc_df <- data.frame(calc_auc(p_all))
    # create list of gene  name + auc for each roi
    all_lvls <- paste(levels(as.factor(top5_all_df$name)), round(all_auc_df$AUC,2), sep=": ")

    # save to file
    match.fun(fileType)(all_outfile)
    print(p_all +
            annotate("text", x = .80, y = .10,
                     label = paste("Avg AUC =", round(mean((calc_auc(p_all))$AUC), 2))) +
            scale_color_hue(labels = c(all_lvls)) +
            labs(color="AUC")

    )

    dev.off()

  } else {
    ### plot top 5 of all genes ##

    # subset counts table for top 5's
    top5_all_counts <- df[colnames(df) %in% top5_all]

    #  melt dataframes into long format, one row per gene per count
    top5_all_df <- melt_roc(data=df,
                            d="level",
                            m=top5_all_list)

    # setup plot variables
    all_title ="ROC: Top 5 All Genes"
    # create outfile path from input outdir
    all_outfile = paste0(outdir, "/Top5_all_ROC.", fileType)

    #  plot roc
    p_all <- ggplot(top5_all_df, aes(d = D, m = M, color = name)) +
      geom_roc(show.legend = FALSE) +
      geom_abline(slope = 1, intercept = 0, color = 'gray') +
      geom_roc(n.cuts = 0, lwd = 1) +
      ggtitle(all_title) +
      scale_x_continuous("Specificity") +
      scale_y_continuous("Sensitivity")
    if (is.null(preset_theme)) {
      p_all <- p_all + theme_light(base_size = 14)
    }

    # # pull out AUC for each gene
    all_auc_df <- data.frame(calc_auc(p_all))
    # create list of gene  name + auc for each roi
    all_lvls <- paste(levels(as.factor(top5_all_df$name)), round(all_auc_df$AUC,2), sep=": ")

    # save to file
    match.fun(fileType)(all_outfile, height = 8, width = 8, units = "in", res = 120)
    print(p_all +
            annotate("text", x = .80, y = .10,
                     label = paste("Avg AUC =", round(mean((calc_auc(p_all))$AUC), 2))) +
            scale_color_hue(labels = c(all_lvls)) +
            labs(color="AUC")

    )

    dev.off()

    ### plot top 5 from gene groups of interest ###
    for(i in seq_along(gene_group)) {

      # pull out top 5 from genes of interest by abs(coeff)
      top5_geneGroup <- fit_df[fit_df$gene %in% unlist(genes[i]),]
      top5_geneGroup <- top5_geneGroup[order(abs(top5_geneGroup$coeff), decreasing=TRUE),] # sort by abs(coeff)
      top5_geneGroup <- head(top5_geneGroup, 5) # pull out the top 5
      top5_geneGroup_lvl <- as.vector(c(top5_geneGroup$gene, "level")) # create list of columns to subset
      top5_geneGroup_list <- as.vector(top5_geneGroup$gene)
      top5_group_counts <- df[colnames(df) %in% top5_geneGroup_lvl] # subset counts dataframe to feed to ROC

      # melt dataframe into long format, one row per gene per count
      top5_geneGroup_df <- melt_roc(data=df,
                                    d="level",
                                    m=top5_geneGroup_list)
      # setup plot variables
      plot_title = str_wrap(paste0("ROC: ", gsub('\\.',' ',gene_group[i])), wrap_num)
      # create outfile path from input outdir
      outfile = paste(outdir, paste0(gene_group[i], ".Top5_ROC.", fileType), sep="/")

      #  plot roc
      p <- ggplot(top5_geneGroup_df, aes(d = D, m = M, color = name)) +
        geom_roc(show.legend = FALSE) +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_roc(n.cuts = 0, lwd = 1) +
        ggtitle(plot_title) +
        scale_x_continuous("Specificity") +
        scale_y_continuous("Sensitivity")
      if (is.null(preset_theme)) {
        p <- p + theme_light(base_size = 14)
      }
      # # pull out AUC for each gene
      auc_df <- data.frame(calc_auc(p))
      # create list of gene  name + auc for each roi
      lvls <- paste(levels(as.factor(top5_geneGroup_df$name)), round(auc_df$AUC,2), sep=": ")


      # save to file
      match.fun(fileType)(outfile, height = 8, width = 8, units = "in", res = 120)
      print(p +
              annotate("text", x = .80, y = .10,
                       label = paste("Avg AUC =", round(mean((calc_auc(p))$AUC), 2))) +
              scale_color_hue(labels = c(lvls)) +
              labs(color="AUC")
      )

      dev.off()
    }
  }
}
