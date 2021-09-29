#' Volcano Plot
#'
#' Takes results of differential expression analysis and generates a volcano plot
#'
#' @param de_results output of dsp_de_analysis function with list of up and downregulated targets, e.g. genes or gene sets
#' @param negative_label label for the negative x axis
#' @param positive_label label for the positive x axis
#' @param target_group name of target groups, or , used for adding color contrast to graph
#' @param targets a list of targets, e.g. genes, within the target_group, used for subsetting graph
#' @param outfile name of file for output
#' @param point_color a specific color to be used for plotting points, should be updated to a gradient method
#' @param top_method column within the de_results that should be used to select the top genes, selecting highest values from the column. Use Significance rather than Pval for this purpose
#' @param n_targets number of genes to show in the cloud
#' @param target_ID column within the de_results data frame that contains the target names to be used for labeling graphs
#' @param FC_ID column within the de_results data frame that contains the FC x-axis values to be used for plotting
#' @param Pval_ID column within the de_results data frame that contains the P-value y-axis values to be used for plotting
#' @param color_list the targetset color values if multiple targetsets are being used, can be NULL if only one targetset passed
#' @param save_plot boolean of whether to save the plot as a file or just output to the console, defaults to TRUE and uses config information for file format
#'
#' @return image file of volcano plot
#' @return plotly interactive volcano plot
#'
#' @examples
#'  plot_volcano(de_results, base_level)
#'
#' @export plot_volcano
#' @export change_axis_revlog_trans

# recode to read in config as one function, then pass args to rest of script
plot_volcano <- function(de = NULL,
                         negative_label = NULL,
                         positive_label = NULL,
                         target_group = "All Probes",
                         fav_targets = NULL,
                         targets = NULL,
                         outfile = 'VolcanoePlot-AllGenes',
                         point_color = 'RdPu',
                         top_method = 'Significance',
                         n_targets = 15,
                         plt_title = 'Volcano Plot for All Genes',
                         target_ID = 'gene',
                         FC_ID = 'FC',
                         Pval_ID = 'Pval',
                         color_list = NULL,
                         save_plot = TRUE,
                         show_legend = TRUE) {
  if(is.null(de)) {
    stop("No dataset provided to plot Volcano Plot, please check inputs and workspace\n")
  } else {

    maxFC <- max(abs(de[[FC_ID]]))
    maxPval <- min(de[[Pval_ID]])

    # find closest FDR value to fdr_cutoff and use that pvalue to add y axis cutoff line
    if(!is.null(fdr_cutoff) & "fdr" %in% colnames(de)){
      fdr_pval <- mean(de[[Pval_ID]][which(abs(de$fdr - fdr_cutoff) ==
                                         min(abs(de$fdr - fdr_cutoff)))])
    }else{
      fdr_cutoff <- NULL
      fdr_pval <- NULL
    }

    # create basic volcano plot with correct formatting
    ggFigure <- ggplot(de, aes(x=eval(parse(text=FC_ID)), y=eval(parse(text=Pval_ID))))+
      geom_point(color=point_color)+
      labs(y="Pvalue",
           x=paste(negative_label, "<-", "FC(FC)", "->", positive_label, sep=" "),
           title = plt_title)+
      scale_x_continuous(limits=c(-maxFC, maxFC))

    # this makes for easier testing if not running in DSPDA, will flip yaxis of graph
    # scale_y_continuous(trans=change_axis_revlog_trans(base=10),
    #                    labels=function(x) format(x, trim=TRUE, digits=4,
    #                                              scientific=ifelse(maxPval < 0.0001, TRUE, FALSE),
    #                                              drop0trailing=TRUE))

    if(show_legend == FALSE){
      ggFigure <- ggFigure + theme(legend.position="none")
    }

    # subset de to only include genes either in specified target groups or above pval/fdr threshold
    if(volcano_color == "Gene Group"){
      gene_coloring <- de[de[[target_ID]] %in% unlist(genes),]
      gene_coloring$Target_coloring <- "Not Specified"
      for(t in target_group){
        gene_coloring$Target_coloring[gene_coloring[[target_ID]] %in% genes[[t]] &
                                        gene_coloring$Target_coloring != "Not Specified"] <- "Multiple"
        gene_coloring$Target_coloring[gene_coloring[[target_ID]] %in% genes[[t]] &
                                        gene_coloring$Target_coloring == "Not Specified"] <- t
      }
      gene_coloring$Target_coloring <- str_wrap(gene_coloring$Target_coloring, width=45)

      color_label <- "Target Group\nMembership"

      color_options <- colors$gene_groups
    }else{
      color_options <- colors[[grouping_var]]
      color_options <- color_options[which(names(color_options) %in% c(negative_label, positive_label))]
      if(!is.null(fdr_cutoff) & !is.null(pval_cutoff)){
        if(fdr_pval < pval_cutoff){
          gene_coloring <- de[which(de[[Pval_ID]] < pval_cutoff),]
          label_thresh_low <- paste("pval <", pval_cutoff)
          label_thresh_high <- paste("FDR <", fdr_cutoff)
          high_thresh <- fdr_pval
        }else{
          gene_coloring <- de[which(de[[Pval_ID]] < fdr_pval),]
          label_thresh_low <- paste("FDR <", fdr_cutoff)
          label_thresh_high <- paste("pval <", pval_cutoff)
          high_thresh <- pval_cutoff
        }

        # label points as positive or negative FC
        gene_coloring$Target_coloring <- ifelse(test=gene_coloring[[FC_ID]] < 0,
                                                yes=negative_label,
                                                no=positive_label)

        # label points as above pval or FDR threshold
        gene_coloring$Target_coloring <- ifelse(test=gene_coloring[[Pval_ID]] >= high_thresh,
                                                yes=paste(label_thresh_low, gene_coloring$Target_coloring),
                                                no=paste(label_thresh_high, gene_coloring$Target_coloring))

        # make color options have muted colors for higher threshold
        color_options <- c(color_options, muted(color_options, l=80))
        names(color_options) <- c(paste(label_thresh_high, negative_label),
                                  paste(label_thresh_high, positive_label),
                                  paste(label_thresh_low, negative_label),
                                  paste(label_thresh_low, positive_label))
      }else{
        if(is.null(fdr_cutoff)){
          gene_coloring <- de[which(de[[Pval_ID]] < pval_cutoff),]
          label_thresh <- paste("pval <", pval_cutoff)
        }else{
          gene_coloring <- de[which(de$fdr < fdr_cutoff),]
          label_thresh <- paste("FDR <", fdr_cutoff)
        }
        gene_coloring$Target_coloring <- ifelse(test=gene_coloring[[FC_ID]] < 0,
                                                yes=paste(label_thresh, negative_label),
                                                no=paste(label_thresh, positive_label))

        names(color_options) <- c(paste(label_thresh, negative_label),
                                  paste(label_thresh, positive_label))
      }

      color_label <- "Significance:"

      # color by fc_cutoff if fc_cutoff is not NULL
      if(!is.null(fc_cutoff)){
        gene_coloring$Target_coloring[abs(gene_coloring[[FC_ID]]) < fc_cutoff] <- paste("FC <", fc_cutoff)
        color_options <- c(color_options, fc_color)
        names(color_options)[length(color_options)] <- paste("FC <", fc_cutoff)
      }
    }

    color_options <- c(color_options, default_color)
    names(color_options)[length(color_options)] <- "Not Specified"

    # add coloring to ggplot
    ggFigure <- ggFigure + geom_point(data=gene_coloring, aes(x=eval(parse(text=FC_ID)), y=eval(parse(text=Pval_ID)), color=Target_coloring))+
      labs(color=color_label)+
      scale_color_manual(values=color_options)

    # add threshold line values to y axis
    yaxis <- data.frame(brk=as.numeric(pretty_breaks(n=4)(0:max(-log10(de[[Pval_ID]])))))
    yaxis$brk <- 10^-(yaxis$brk)

    # keep scientific notation if small enough pvalues when changing to character
    yaxis$label <- format(yaxis$brk, trim=TRUE, digits=4,
                          scientific=ifelse(maxPval < 0.0001, TRUE, FALSE),
                          drop0trailing=TRUE)

    # add threshold lines if thresholds are not NULL
    if(!is.null(fc_cutoff)){
      ggFigure <- ggFigure + geom_vline(xintercept=fc_cutoff, linetype="dotted")+
        geom_vline(xintercept=-fc_cutoff, linetype="dotted")
        # annotate("text", x=fc_cutoff+0.35, y=1,
        #          label=paste0("FC=", round(fc_cutoff, digits=2)))
    }
    if(!is.null(pval_cutoff)){
      ggFigure <- ggFigure + geom_hline(yintercept=pval_cutoff, linetype="dotted")
      yaxis <- rbind(yaxis, c(pval_cutoff, paste0('pval=',pval_cutoff)))
    }
    if(!is.null(fdr_cutoff)){
      ggFigure <- ggFigure + geom_hline(yintercept=fdr_pval, linetype="dotted")
      yaxis <- rbind(yaxis, c(fdr_pval, paste0('FDR=',fdr_cutoff)))
    }

    # order yaxis in increasing value
    yaxis$brk <- as.numeric(yaxis$brk)
    yaxis <- yaxis[order(yaxis$brk, decreasing=F),]

    # subset de to only contain genes to label on plot, either by user specified genes or top n_targets by pval
    if(volcano_label == "Fav Genes"){
      gene_labels <- subset(de, subset=eval(parse(text=target_ID)) %in% fav_targets)
    }else{
      # remove gene labels for genes below lowest pvalue cutoff
      gene_labels <- de[which(de[[Pval_ID]] < min(fdr_pval, pval_cutoff, na.rm = T)),]

      # only label genes above fc_cutoff if set by user else only look at pvalue
      if(!is.null(fc_cutoff) & label_fc == FALSE){
        gene_labels <- gene_labels[which(abs(gene_labels[[FC_ID]]) > fc_cutoff),]
      }

      # only keep top # of genes by pvalue
      gene_labels <- gene_labels[head(order(gene_labels[[Pval_ID]], decreasing=FALSE), n=n_targets),]
    }

    # add gene labels to ggplot
    ggFigure <- ggFigure + geom_text_repel(data=gene_labels, aes(x=eval(parse(text=FC_ID)),
                                                                 y=eval(parse(text=Pval_ID)),
                                                                 label=eval(parse(text=target_ID))),
                                           force=5, fontface="bold", min.segment.length=0.1,
                                           size=label_size)

    # add y axis labels
    ggFigure <- ggFigure + scale_y_continuous(trans=change_axis_revlog_trans(base=10), breaks=as.numeric(yaxis$brk),
                                              labels=yaxis$label)

    if("leadingEdge" %in% colnames(gene_labels)){
      gene_labels <- gene_labels[,-which(colnames(gene_labels) == "leadingEdge")]
    }

    file_name <- str_split(outfile, pattern = "_", simplify = T)
    file_name <- file_name[-length(file_name)]
    file_name <- paste(file_name, collapse = "_")
    write.table(x = gene_labels, file = paste(file_name, "volcano_labeled_points.tsv"),
                sep = "\t", quote = F, col.names = NA, row.names = T)

    # save to file
    if(save_plot) {
      if(fileType %in% c('tiff','png','jpeg','bmp')) {
        match.fun(fileType)(outfile, width = ifelse(maxFC > 2, 1600, 1300),
                            height = 1000, res = 150) # it's not print to png!!!!
      } else {
        match.fun(fileType)(outfile)
      }
      print(ggFigure)
      dev.off()
    }
  }
}

#' change_axis_revlog_trans
#'
#' reverse log transform axis; used to return pvalue rather than -log10(pvalue) on yaxis
#' @param base base in which logs are computed
#' @return revlog_trans reverse log transformation
#' @export
change_axis_revlog_trans <- function(base=exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  revlog_trans <- trans_new(name=paste0("revlog-", base),
                            transform=trans,
                            inverse=inv,
                            breaks=log_breaks(base=base),
                            domain=c(1e-100, Inf))

  return(revlog_trans)
}
