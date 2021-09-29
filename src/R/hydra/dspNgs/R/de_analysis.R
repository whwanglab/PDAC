#' DE Analysis
#'
#' Runs a linear mixed effects model on the dataset to analyze differential expression using the LME4 package
#'
#' @param df negative normalized DSP count data transposed to wide format, one column per gene
#' @param control_var column from annotations to use as experimental control variable
#' @param grouping_var column from annotations to test during DE analysis
#' @param base_level variable to use as base level for testing, all other vars looped through during analysis
#' @param residuals_dir path to output directory image subfolder for residual images
#' @param qq_dir path to output directory image subfolder for qq plots
#' @param pval_cutoff pvalue cutoff for significance from config file
#' @param make_plots logical as to whether qq & residual plots should be written
#' @param random_slope logical to turn on or off random slopes depending on study
#'
#' @return PNG file of top 5 genes ROC curve
#' @return dataframe of AUC for each gene
#'
#' @examples
#'  run_DE(df, control_var, base_level, genes,residuals_dir, qq_dir, pval_cutoff)
#'
#' @export run_DE
#' @export de_genes_by_group

run_DE <- function(df, control_var, grouping_var, base_level, residuals_dir, qq_dir, pval_cutoff,
                   make_plots = TRUE, random_slope = TRUE){

  # create list of genes to loop through
  gene_list <- colnames(df)[!colnames(df) %in% c("Sample_ID", grouping_var, control_var, "level")]

  # create empty matrix to store DE results, preallocating space for efficiency
  de_res <- data.frame(matrix(ncol=10, nrow=0))

  # order levels so that base_level is first, catch errors
  lvls <- levels(as.factor(df[[grouping_var]]))
  if(!base_level %in% lvls) {
    stop('Base level not provided as part of selected column for DE tested.\nPlease confirm that base_level is in the specified column\n')
  } else if(length(lvls) <= 1) {
    stop('The column specified doesn\'t have enough levels to run DE analysis.\nPlease select a different column\n')
  }
  lvls <- lvls[c(which(lvls == base_level),
                 which(lvls != base_level))]
  df[[grouping_var]] <- factor(x = df[[grouping_var]], levels = lvls)


  # loop through levels
  for (lvl in lvls[-1]) {
    ind <- df[[grouping_var]] %in% c(base_level, lvl)
    test_df <- df[ind, ]
    if(is.null(control_var)){
      n_ctrl <- 1
    }else{
      n_ctrl <- length(unique(test_df[[control_var]]))
    }

    qq.plots <- vector(length(gene_list), mode='list')
    res.plots <- vector(length(gene_list), mode='list')
    counter <- 0
    if (n_ctrl > 1) {
      cat(paste0('\nRunning Mixed Effect Models for DE (', grouping_var, '): ', base_level, ' vs ', lvl, '\n'))
    } else {
      cat(paste0('\nRunning Standard Model for DE (', grouping_var, '): ', base_level, ' vs ', lvl, '\n'))
    }
    pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3)
    # loop through each gene in the list
    for (gene in gene_list){
      counter <- counter + 1
      setTxtProgressBar(pb = pb, value = counter)

      # test the number of levels in control_var, if 1 use lm, if > 1 use lmer
      if(n_ctrl > 1) {
        # mod <- lmer(gene ~ grouping_var + (1 + grouping_var|control_var), data=df)
        formula = paste0(gene, " ~ `", grouping_var, "` + (1 + `", grouping_var, "` | `", control_var, "`)")
        mod <- lmer( as.formula(formula), data=test_df,
                     # hard prevent singular boundary warnings:
                     control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

        # summary stats dataframe
        cf <- data.frame(coefficients(summary(mod)))[2,]
        cf$gene <- gene
        colnames(cf)[1] <- "FC"
        colnames(cf)[5] <- "Pval"
        cf$Significance <- -log10(cf$Pval)
        rownames(cf) <- c()
      } else {
        formula = paste0(gene, " ~ `", grouping_var, "`")
        mod <- lm(as.formula(formula), data = test_df)

        # summary table (mimic lmer in case both used)
        cf <- data.frame(coefficients(summary(mod)))[2,]
        colnames(cf)[1] <- "FC"
        colnames(cf)[4] <- "Pval"
        cf$df <- NA
        cf <- cf[,c(1,2,5,3,4)]
        cf$gene <- gene
        cf$Significance <- -log10(cf$Pval)
      }
      # add DE results to de_res dataframe, giving unique row id
      if(gene == gene_list[1]) {
        de_res <- cf
      } else {
        de_res <- rbind(de_res, cf)
      }
      if(draw_DE_plots) {
        # plot residuals versus fit
        res_p <- plot(mod, type = c("p", "smooth"),
                      main=paste0("Residuals vs Fit: ", gene),
                      ylab="Residuals",
                      xlab="Fit")

        res.plots[[counter]] <- res_p

	    # QQ Plots
        qq_p <- qqmath(mod, id = pval_cutoff,
                       main = list(paste0("QQ Plot: ", gene), cex = 1))

        qq.plots[[counter]] <- qq_p
      }
    }
    if(draw_DE_plots) {
      counter <- 0
      pb <- txtProgressBar(min = 0, max = length(qq.plots), style = 3)
      cat("\nDrawing qq plots\n")
      pdf(paste0(models_dir, "/qq_plots.pdf"), onefile = TRUE)
      for (plot in qq.plots){
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter)
        print(plot)
      }
      dev.off()

      pb <- txtProgressBar(min = 0, max = length(res.plots), style = 3)
      cat("\nDrawing residual plots\n")
      counter <- 0
      pdf(paste0(models_dir, "/residuals_plots.pdf"), onefile = TRUE)
      for (plot in res.plots){
        counter <- counter + 1
        setTxtProgressBar(pb = pb, value = counter)
        print(plot)
      }
      dev.off()
    }
    # add column names back to dataframe
    colnames(de_res) <- colnames(cf)
    de_res$test <- paste(base_level, 'vs', lvl)
    # adjust p-values
    de_res$fdr <- p.adjust(de_res$Pval, 'fdr')
    de_res$fwer <- p.adjust(de_res$Pval, 'bonferroni')

    # create output matrix or append new results
    if(lvl == lvls[2]) {
      de_res_out <- de_res
    } else {
      de_res_out <- rbind(de_res_out, de_res)
    }
  }
  return(de_res_out)
}


#' @title de_genes_by_group
#'
#' Function used to determine unique de genes by de group
#'
#' @param de de results dataframe defaults to de_results
#' @param pval pval cutoff defaults to pval_cutoff
#' @param fc fc cutoff defaults to fc_cutoff
#' @param group grouping variable defaults to grouping_var
#' @param SP Segment Properties defaults to dfs[[3]]
#'
#' @return list of de genes unique to each de group
#'
#' @examples
#'
#' de_genes_by_group(de = de_results, pval = pval_cutoff, fc = fc_cutoff, group = grouping_var, SP = dfs[[3]])
#'
de_genes_by_group <- function(de = de_results, pval = pval_cutoff, fc = fc_cutoff, group = grouping_var, SP = dfs[[3]]){
  tmp_de <- de[which(de$Pval < pval),]
  tests <- as.data.frame(str_split(tmp_de$test, pattern = " vs ", simplify = T))
  colnames(tests) <- c("Negative", "Positive")

  genes_de <- list()

  for(i in unique(SP[[group]])){
    neg <- grep(tests$Negative,pattern = i)
    pos <- grep(tests$Positive,pattern = i)
    neg_de <- tmp_de[neg,]
    pos_de <- tmp_de[pos,]

    neg_genes <- subset(neg_de, (neg_de$Pval < pval & neg_de$FC <= -fc))$gene
    pos_genes <- subset(pos_de, (pos_de$Pval < pval & pos_de$FC >= fc))$gene

    genes_de[[i]] <- unique(c(neg_genes, pos_genes))

    for(name in names(genes_de)){
      if(name != i){
        i_remove <- which(genes_de[[i]] %in% genes_de[[name]])
        name_remove <- which(genes_de[[name]] %in% genes_de[[i]])
        if(length(i_remove) > 0){
          genes_de[[i]] <- genes_de[[i]][-i_remove]
          genes_de[[name]] <- genes_de[[name]][-name_remove]
        }
      }
    }
  }

  return(genes_de)
}
