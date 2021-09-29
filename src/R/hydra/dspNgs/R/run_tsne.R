#' @title run_tsne
#'
#' Runs t-Distributed Stochastic Neighbor Embedding (tsne) on dsp data
#'
#' @param df read in dsp excel workbook dataset
#' @param outdir output directory for figures
#' @param outlier_cutoff multiplier for outlier cutoff
#' @param color category to color the graph by
#' @param symbol category to change symbols in the graph
#' @param size should the size of the points change by gene_count
#'
#' @return none
#'
#' @examples
#' run_tsne(df, "file/path", outlier_cutoff, color = "slide_id", symbol = "TIS", size = TRUE)
#'
#' @export run_tsne
#' @export run_tsne_genes

run_tsne <- function(df, pca = NULL, outdir, color = "cluster", symbol = "dsp_slide",
                     size = TRUE, perplexity = 30, colors = colors){

  #read in counts without outliers and notes with slide name split
  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]
  set.seed(100)

  if ((!any(names(samp_notes) == color) & color != "cluster") |
      (!any(names(samp_notes) == symbol) & symbol != "cluster")){
    stop(str_wrap(paste("Given coloring or symbol factor was not allowed. OPTIONS: ",
                        paste(names(samp_notes), collapse = ", ")), width = 80), call.=FALSE)
  }


  if (is.null(pca)){
    #run pca
    flog.info(paste0("Running PCA"), name="DSP_NGS_log")
    pca <- run_PCA(df, outdir)
  }

  flog.info(paste0("Running tSNE"), name="DSP_NGS_log")
  #calculate perplexity based on number of samples
  maxPerplex <- (nrow(pca$call$X) - 1) %/% 3

  if (perplexity > maxPerplex){
    perplexity <- maxPerplex
    flog.warn(paste0("Given perplexity is larger than data allows, continuing with max perplexity allowed: ",
                     maxPerplex, "\n"), name="DSP_NGS_log")
    print(paste0("WARNING: Given perplexity is larger than data allows, continuing with max perplexity allowed: ",
                     maxPerplex))
  }

  tsne <- Rtsne(t(norm_counts), pca = TRUE, check_duplicates = FALSE, perplexity = perplexity)


  #add tsne info to samp_notes for ease of plotting
  samp_notes$tsne <- data.frame(tsne$Y)
  samp_notes$gene_count <- colSums(norm_counts[,match(samp_notes$Sample_ID, colnames(norm_counts))])
  samp_notes$cluster <- pca$call$X$clust

  #check options for function
  color_by <- as.character(eval(parse(text=paste0("samp_notes[[\"", color, "\"]]"))))
  if (is.null(color_by)){
    stop("Given coloring factor was not allowed, default: microsat
         options: gene_count, cluster, TIS, slide_id", call.=FALSE)
  }
  symbol_by <- as.character(eval(parse(text=paste0("samp_notes[[\"", symbol, "\"]]"))))
  if (is.null(symbol_by)){
    stop("Given symbol factor was not allowed, default: TIS
         options: gene_count, cluster, microsat, slide_id", call.=FALSE)
  }
  if (size == TRUE){
    size <- samp_notes$gene_count / 50000
  }else {
    size <- 6
  }

  clrs <- unlist(colors[color], use.names = FALSE)
  names(clrs) <- names(colors[[color]])

  tsne <- ggplot(data = samp_notes, aes(x = tsne$X1, y = tsne$X2, colour = str_wrap(color_by, wrap_num),
                                       text = paste("ROI: ", Sample_ID,"\nGene Count: ", round(gene_count),
                                                    "\nColor: ", color_by, "\nSymbol:", symbol_by))) +
    geom_point(size = size, aes(shape = symbol_by)) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "tSNE1", y = "tSNE2", title = "tSNE", colour = color, shape = symbol) +
    scale_color_manual(values = clrs)

  ggsave(paste0(outdir, "/tSNE.", fileType), tsne, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot tSNE
  print(ggplotly(tsne, tooltip = "text"))

  return(samp_notes)
}

run_tsne_genes <- function(df, pca = NULL, outdir, size = TRUE, perplexity = 90, de = de_results, expressed_genes = de_genes){

  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]
  norm_counts <- t(norm_counts)
  norm_counts <- norm_counts[,match(row.names(pca$call$X), colnames(norm_counts))]
  set.seed(100)

  if (is.null(pca)){
    #run pca
    flog.info(paste0("Running PCA"), name="DSP_NGS_log")
    pca <- run_PCA_genes(df, outdir, size = size)
  }

  flog.info(paste0("Running tSNE for genes"), name="DSP_NGS_log")
  #calculate perplexity based on number of samples
  maxPerplex <- (nrow(pca$call$X) - 1) %/% 3

  if (perplexity > maxPerplex){
    perplexity <- maxPerplex
    flog.warn(paste0("Given perplexity is larger than data allows, continuing with max perplexity allowed: ",
                     maxPerplex, "\n"), name="DSP_NGS_log")
    print(paste0("WARNING: Given perplexity is larger than data allows, continuing with max perplexity allowed: ",
                     maxPerplex))
  }

  tsne <- Rtsne(t(norm_counts), pca = TRUE, check_duplicates = FALSE, perplexity = perplexity)


  #add tsne info to samp_notes for ease of plotting
  info <- data.frame(tsne$Y)
  info$cluster <- pca$call$X$clust
  info$gene <- colnames(norm_counts)
  info$gene_count <- colSums(norm_counts)
  info$color <- "Not Specified"

  for (i in 1:length(gene_group)){
    w2kp <- which(info$gene %in% unlist(genes[i]))
    w2kp2 <- which(info$color[w2kp] != "Not Specified")

    info$color[w2kp] <- names(genes[i])
    info$color[w2kp2] <- "Multiple"
  }

  w2kp <- which(info$gene %in% fav_genes)
  info$color[w2kp] <- "Fav.Genes"

  #check options for function
  if (size == TRUE){
    sizes <- info$gene_count / 40000
    w2kp <- which(sizes < 0.15)
    sizes[w2kp] <- 0.15
  }else {
    sizes <- ifelse(info$color == "Not Specified", 0.8, 2)
  }

  info$size <- sizes

  clrs <- unlist(colors["gene_groups"], use.names = FALSE)
  names(clrs) <- names(colors[["gene_groups"]])

  tsne <- ggplot(data = info, aes(x = X1, y = X2, colour = str_wrap(color, wrap_num),
                                  text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                               "\n", "Gene Group",": ", color, sep = ""))) +
    geom_point(size = info$size) +
    geom_point(data = subset(info, color != "Not Specified" & color != "All Probes"), size =  info$size[info$color!="Not Specified" & info$color != "All Probes"]) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "tSNE1", y = "tSNE2", title = "tSNE Genes", colour = "Gene Group", shape = "Cluster") +
    scale_color_manual(values = clrs)+
    theme(legend.position="right", legend.direction = "vertical")

  ggsave(paste0(outdir, "/tSNE_genes.", fileType), tsne, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot tSNE
  print(ggplotly(tsne, tooltip = "text"))

  zero <- ifelse(test = sum(lengths(expressed_genes)) > 0, yes = FALSE, no = TRUE)

  info$group_color <- "Not Specified"
  if(zero == FALSE){
    for (i in 1:length(expressed_genes)){
      w2kp <- which(info$gene %in% unlist(expressed_genes[i]))
      w2kp2 <- which(info$group_color[w2kp] != "Not Specified")

      info$group_color[w2kp] <- names(expressed_genes[i])
      info$group_color[w2kp2] <- "Multiple"
    }

    w2kp <- which(info$group_color == "Multiple")
    info$group_color[w2kp] <- "Not Specified"

    if (size == TRUE){
      sizes <- info$gene_count / 40000
      w2kp <- which(sizes < 0.15)
      sizes[w2kp] <- 0.15
    }else {
      sizes <- ifelse(info$group_color == "Not Specified", 0.8, 2)
    }

    info$size <- sizes

    clrs <- unlist(colors[grouping_var], use.names = FALSE)
    names(clrs) <- names(colors[[grouping_var]])

    clrs[['Not Specified']] <- '#ADADAD'

    tsne_plot <- ggplot(data = info, aes(x = X1, y = X2, colour = group_color,
                                       text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                                    "\n", grouping_var,": ", group_color, sep = ""))) +
      geom_point(size = info$size) +
      geom_point(data = subset(info, group_color != "Not Specified"), size = info$size[info$group_color!="Not Specified"]) +
      scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
      labs(x = "tSNE1", y = "tSNE2", title = "tSNE Genes", colour = grouping_var) +
      scale_color_manual(values = clrs)

    ggsave(paste0(outdir, "/tSNE_genes_", grouping_var, ".", fileType), tsne_plot, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

    #plot info
    print(ggplotly(tsne_plot, tooltip = "text"))
  }

  return(info)
}
