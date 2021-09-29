#' @title QC_plots
#'
#' Creates QC plots of normalization factors and distributions of counts
#'
#' @param data_matrix the data matrix list
#' @param norm_method the selected norm method for output data matrix
#' @param grp_var the selected column name for the grouping analysis
#' @param ctrl_var the selected column name to control by for DE
#' @param HKs the list of housekeeping genes
#' @param outdir the output directory for plots
#' @param colors the list of specified colors
#' @param txt_type boolean specifying if txt files from Azorius were the input
#'
#' @return norm_counts chosen normalization method
#'
#' @examples
#'  qc_plots(data_matrix = dfs, norm_method = "Q3", outdir = od)
#'
#'
#' @export qc_plots

# Create QC plots
qc_plots <- function(data_matrix,
                       norm_method = NULL,
                       grp_var = NULL,
                       ctrl_var = NULL,
                       HKs = NULL,
                       outdir,
                       colors,
                       txt_type=TRUE) {

  #####
  # get the data setup
  # set the gene names to be the rownames
  raw_counts <- as.data.frame(data_matrix[[1]])
  rownames(raw_counts) <- raw_counts[,1]
  raw_counts <- raw_counts[,-1]

  # set the annot from samp_notes matrix
  annot <- as.data.frame(data_matrix[[3]])

  # check that all column names match rownames
  if (all(colnames(raw_counts) %in% annot$Sample_ID) == FALSE){
    stop('Sample_ID in raw counts doesn\'t match Sample_ID in annotations.')
  }

  # reorder to match annot
  raw_counts <- raw_counts[,match(annot$Sample_ID, colnames(raw_counts))]

  # if there are NA's replace the "NA" with NA
  annot <- replace(annot, annot == "NA", NA)

  # do a check for the neg and hk norm factors
  ifelse(!"NormFactorHK" %in% colnames(annot) || any(is.na(annot$NormFactorHK)), do_hk <- FALSE, do_hk <- TRUE)
  ifelse(any(grep(pattern = "NormFactorNeg_", colnames(annot))), do_neg <- TRUE, do_neg <- FALSE)

  # if no HKs in config overlap with genes in dataset, reset do_hk to FALSE and warn user
  # not applicable for DSP-DA output
  if(length(which(rownames(raw_counts) %in% HKs)) == 0 && txt_type){
    do_hk = FALSE
    cat('Warning: No housekeepers found in the data. HK normalization will not be performed.\n',
        'If you desire HK normalization, adjust hk_genes in the config file to include genes\n',
        'in your data set.')
  }

  # check how many pools in the dataset, so you know how many norm factors to work with
  if(do_neg){
    col_indx <- grep(pattern = "NormFactorNeg_", colnames(annot))
    neg_norm_factor <- c()
    for(pool in 1:length(col_indx)){
      neg_norm_factor <- append(neg_norm_factor, colnames(annot)[col_indx[pool]])
    }
  }

  # if DSP-DA data, calculate the q3 normalization factor from raw counts
  if(!txt_type && !("NormFactorQ3" %in% colnames(annot))) {
    negs_to_remove <-
      data_matrix[[5]][data_matrix[[5]][["CodeClass"]] == "Negative",
                         "TargetName"]
    annot$NormFactorQ3 <-
      nanoNormFactors(
        apply(
          data_matrix[[1]][!data_matrix[[1]][["TargetName"]] %in%
                             negs_to_remove, 2:ncol(data_matrix[[1]])],
          2,
          quantile,
          probs = 0.75,
          na.rm=T))
  }

  # make a list of the names of normalization factors that are available - including the data frames needed (always require Q3)
  # then calculate the normalized counts with factors, except negs - pull from data matrix
  q3 <- data.frame(t(t(raw_counts)/annot$NormFactorQ3))
  # Preserve DSP-DA sample ID format
  colnames(q3) <- colnames(raw_counts)
  ifelse(do_hk, print("Do HK"), print("No HK"))
  ifelse(do_neg, print("Do Neg"), print("No Neg"))
  if(!do_hk && !do_neg){

    norm_factor_list <- c("NormFactorQ3")
    distrib <- list("Q3"= q3)

  }else if(do_neg && !do_hk){

    norm_factor_list <- c("NormFactorQ3", neg_norm_factor)
    nc <- as.data.frame(data_matrix[[2]])
    rownames(nc) <- nc[,1]
    nc <- nc[,-1]
    distrib <- list("NegCtrl"= nc, "Q3"= q3)

  }else if(do_hk && !do_neg){

    norm_factor_list <- c("NormFactorHK", "NormFactorQ3")
    if (txt_type) {
      hk <- data.frame(t(t(raw_counts)/annot$NormFactorHK))
    # Use DSP-DA normalized data as-is
    } else {
      hk <- as.data.frame(data_matrix[[2]])
      rownames(hk) <- hk[, 1]
      hk <- hk[, -1]
    }
    distrib <- list("HK" = hk, "Q3"= q3)

  }else {

    norm_factor_list <- c("NormFactorHK", "NormFactorQ3", neg_norm_factor)
    nc <- as.data.frame(data_matrix[[2]])
    rownames(nc) <- nc[,1]
    nc <- nc[,-1]
    hk <- data.frame(t(t(raw_counts)/annot$NormFactorHK))
    distrib <- list("NegCtrl"= nc,"HK" = hk, "Q3"= q3)

  }

  # Set default plot size
  width <- 12
  height <- 8

  #####
  # plotting HKgeomeans
  ## if(do_hk && txt_type){
  if(FALSE){

    # calculate the HK geomeans based off the list of HKs
    hk_geo <- geometric.mean(raw_counts[rownames(raw_counts) %in% HKs,])

    # just append to annot, we already sorted to match sample id sort
    annot$HK_Geo <- hk_geo

    df_hkgeo <- subset(annot, select = c("Sample_ID", "HK_Geo"))
    p1 <- ggplot(df_hkgeo, aes(x=reorder(Sample_ID, HK_Geo), y = HK_Geo)) +
            geom_point(stat = "identity") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(title = "HK Geomeans for each Sample", x = "Sample ID", y = "HK Geomeans")

    #####
    # plotting norm factors
    # subset for just the hks
    hk_sub <- subset(hk, rownames(hk) %in% HKs)
    HKnormalized <- t(hk_sub)
    norm_factors <- t(annot$NormFactorHK)

    # for each sample, calculate the MSE of its HKs from the average profile
    avg_hk_profile <- apply(HKnormalized,2,mean)
    HK_mse <- rowMeans((sweep(HKnormalized,2,avg_hk_profile))^2)

    # identify high HK_mse outliers:
    outlier_HKMSE <- HK_mse > mean(HK_mse)+2*sd(HK_mse)

    # identify hkgeomean outliers:
    outlier_normfactor <- abs(norm_factors-mean(norm_factors))>2*sd(norm_factors)

    # create a data frame for plotting ease
    df <- as.data.frame(t(rbind(HK_mse, norm_factors)))

    # plot the normalization factor
    p2 <- ggplot(df, aes(x=V2, y=HK_mse)) +
            geom_point() +
            geom_text_repel(data=subset(df, outlier_HKMSE|outlier_normfactor), label = rownames(subset(df, outlier_HKMSE|outlier_normfactor))) +
            labs(title = "HK Norm Factors vs Deviance from Typical HK Profile", x = "Normalization Factor", y = "MSE of sample from mean profile")

    #####
    # plotting the heatmaps of HKs
    hk_raw <- subset(raw_counts, rownames(hk) %in% HKs)
    hk_q3 <- subset(q3, rownames(q3) %in% HKs)
    p3 <- pheatmap(cor(t(hk_raw)), main = "HK Counts in Raw Data")
    p4 <- pheatmap(cor(t(hk_q3)), main = "HK Counts in Q3 Data")

    if (fileType == "pdf" | fileType == "svg"){
      width <- width*0.85
      height <- height*0.85
    }

    # save to file
    match.fun(fileType)(paste0(outdir, "/HKgeomean.", fileType), width = width, height = height, units = "in", res = 120)
    print(p1)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/Norm_factors.", fileType), width = width, height = height, units = "in", res = 120)
    print(p2)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/HK_Heatmap_raw.", fileType), width = width, height = height, units = "in", res = 120)
    print(p3)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/HK_Heatmap_Q3.", fileType), width = width, height = height, units = "in", res = 120)
    print(p4)
    dev.off()
  }

  # plotting the heatmap of HKs and negs
  ## if(do_hk && do_neg){
  if(FALSE){

    hk_neg <- subset(nc, rownames(nc) %in% HKs)
    p5 <- pheatmap(cor(t(hk_neg)), main = "HK Counts in Neg Data")

    # save to file
    match.fun(fileType)(paste0(outdir, "/HK_Heatmap_Neg.", fileType), width = width, height = height, units = "in", res = 120)
    print(p5)
    dev.off()

  }

  ## #####
  ## # plotting pairs plots of norm factors
  ## p6 <- ggpairs(subset(annot, select = norm_factor_list)) +
  ##         labs(title = "Pairs Plots of Norm Factors")

  ## #####
  ## # plotting distributions
  ## # create a data frame with all the distributions together to plot easier
  ## many_distros <- data.frame()
  ## counter <- 1

  ## # for loop to add them together
  ## for( dist in distrib){

  ##   temp_dist <- melt(log2(t(dist)))
  ##   temp_dist <- cbind(temp_dist, NormMethod = rep(names(distrib)[counter], dim(dist)[1]))

  ##   # bind to make a large df
  ##   many_distros <- rbind(many_distros, temp_dist)
  ##   counter <- counter + 1

  ## }

  ## # plot the distributions
  ## p7 <- ggplot(many_distros, aes(x=value, fill = NormMethod)) +
  ##         facet_wrap(~NormMethod, nrow = 3) +
  ##         geom_density(alpha=.5) +
  ##         labs(title = "Distributions of Normalized Counts by Normalization Methods", x = "Log2 Counts" , y = "Density")

  ## # plot count distributions to check how the normalization affects each DE grouping variable
  ## # merge all nrom tables together, add DE grouping variable annotations
  ## all_norm <- rbind(
  ##   reshape2::melt(raw_counts, id.vars=NULL) %>% mutate(norm="Raw"),
  ##   if(do_neg) reshape2::melt(nc, id.vars=NULL) %>% mutate(norm="Neg"),
  ##   if(do_hk) reshape2::melt(hk, id.vars=NULL) %>% mutate(norm="HK"),
  ##   reshape2::melt(q3, id.vars=NULL) %>% mutate(norm="Q3")
  ## )
  ## colnames(all_norm) = c('Sample_ID', 'Count', 'Normalization')
  ## all_norm = merge(all_norm, annot[,c('Sample_ID', grouping_var)])
  ## all_norm$Normalization = factor(all_norm$Normalization, levels = c('Raw', 'Neg', 'HK', 'Q3'))

  ## dodge <- position_dodge(width = 1)
  ## p8 <- ggplot(all_norm, aes(x = Normalization, y = log2(Count), fill = all_norm[,grouping_var])) +
  ##   geom_violin(position = dodge) +
  ##   facet_wrap(~Normalization, nrow = 1, scales = 'free_x') +
  ##   geom_boxplot(position = dodge, width = 0.15, outlier.colour=NA, show.legend = FALSE) +
  ##   labs(fill = 'DE grouping variable')
  ## if (is.null(preset_theme)) {
  ##   p8 <- p8 + theme_minimal(base_size = 16)
  ## }
  ## p8 <- p8+ theme(panel.grid = element_blank(),
  ##         axis.line.x = element_line(),
  ##         axis.line.y = element_line(),
  ##         strip.text = element_blank() )

  ## if (length(unique(all_norm[,grouping_var])) <= 4){
  ##   p8 = p8 + scale_fill_manual(values = c('#61c237', '#376cc2', '#7c7c7c', '#9b65c4'))
  ## }

  ## #####
  ## # plotting the heatmaps of the normalized counts

  ## #restart counter and setup groups
  ## counter <- 1
  ## annotations <- subset(annot, select = c("Sample_ID", ctrl_var, grp_var))
  ## rownames(annotations) <- annotations[,1]
  ## annotations <- annotations[,-1]

  ## # start clean
  ## if (!is.null(dev.list())) {
  ##   dev.off()
  ## }

  ## # for loop to print heatmaps
  ## for( dist in distrib){

  ##   match.fun(fileType)(paste0(outdir, paste0("/", names(distrib)[counter], " Heatmap."), fileType),
  ##                       width = width, height = width, units = "in", res = 120)
  ##   log_dist <- log2(dist)
  ##   # Preserve DSP-DA name format
  ##   colnames(log_dist) <- colnames(dist)
  ##   ph <- pheatmap(log_dist, show_rownames = FALSE, annotation_col = annotations, annotation_colors = colors,
  ##                  main = paste0(names(distrib)[counter], " Normalized Counts"))
  ##   print(ph)
  ##   dev.off()

  ##   counter <- counter + 1

  ## }

  ## # save to file
  ## match.fun(fileType)(paste0(outdir, "/PairsPlots.", fileType), width = width, height = height, units = "in", res = 120)
  ## print(p6)
  ## dev.off()

  ## match.fun(fileType)(paste0(outdir, "/Distributions.", fileType), width = width, height = height, units = "in", res = 120)
  ## print(p7)
  ## dev.off()

  ## match.fun(fileType)(paste0(outdir, "/DistributionsByDEgrouping.", fileType), width = width*1.3, height = height*1.3, units = "in", res = 120)
  ## print(p8)
  ## dev.off()

  #####
  # return the selected norm method
  if(is.null(norm_method)){
    return()
  } else if (!txt_type) {
    # Bypass renormalization for DSP-DA data
    return(data_matrix[[2]])
  }else if(norm_method == "Q3"){
    norm_counts <- cbind("Gene" = rownames(q3), q3)
    rownames(norm_counts) <- NULL

  }else if(norm_method == "HK"){
    norm_counts <- cbind("Gene" = rownames(hk), hk)
    rownames(norm_counts) <- NULL

  }else if(norm_method == "Neg"){
    norm_counts <- cbind("Gene" = rownames(nc), nc)
    rownames(norm_counts) <- NULL
  }

  return(norm_counts)

}
