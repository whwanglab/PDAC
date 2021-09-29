#' @title run_umap
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
#' run_umap(df, "file/path", outlier_cutoff, color = "slide_id", symbol = "TIS", size = TRUE)
#'
#' @export run_umap
#' @export run_umap_genes

run_umap <- function(df, pca = NULL, outdir, color = "cluster", symbol = "dsp_slide",
                     size = TRUE, colors = colors){

  #read in counts without outliers and notes with slide name split
  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]

  if ((!any(names(samp_notes) == color) & color != "cluster") |
      (!any(names(samp_notes) == symbol) & symbol != "cluster")){
    stop(str_wrap(paste("Given coloring or symbol factor was not allowed. OPTIONS: ",
                        paste(names(samp_notes), collapse = ", ")), width = 80), call.=FALSE)
  }

  if (is.null(pca)){
    flog.info(paste0("Running PCA"), name="DSP_NGS_log")
    #run pca
    pca <- run_PCA(df, outdir)
  }

  #calculate UMAP
  flog.info(paste0("Calculating UMAP"), name="DSP_NGS_log")

  if (ncol(norm_counts) < 15){
    map <- umap(as.matrix(t(norm_counts)), n_neighbors= ncol(norm_counts)-2)
  }else{
    map <- umap(as.matrix(t(norm_counts)))
  }


  #add umap info to samp_notes for ease of plotting
  samp_notes$umap <- data.frame(map$layout)
  samp_notes$gene_count <- colSums(norm_counts)
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
    size = samp_notes$gene_count / 50000
  }else {
    size = 6
  }

  clrs <- unlist(colors[color], use.names = FALSE)
  names(clrs) <- names(colors[[color]])

  umap <- ggplot(data = samp_notes, aes(x = umap$X1, y = umap$X2, colour = str_wrap(color_by, wrap_num),
                                        text = paste("ROI: ", Sample_ID,"\nGene Count: ",
                                                     round(gene_count), "\nColor: ", color_by,
                                                     "\nSymbol:", symbol_by))) +
    geom_point(size = size, aes(shape = symbol_by)) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "UMAP1", y = "UMAP2", title = "UMAP", colour = color, shape = symbol) +
    scale_color_manual(values = clrs)

  ggsave(paste0(outdir, "/UMAP.", fileType), umap, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot tSNE
  print(ggplotly(umap, tooltip = "text"))

  return(samp_notes)

}

run_umap_genes <- function(df, pca = NULL, outdir, size = TRUE, de = de_results, expressed_genes = de_genes){

  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]
  norm_counts <- t(norm_counts)
  norm_counts <- norm_counts[,match(row.names(pca$call$X), colnames(norm_counts))]

  if (is.null(pca)){
    flog.info(paste0("Running PCA"), name="DSP_NGS_log")
    #run pca
    pca <- run_PCA_genes(df, outdir)
  }

  #calculate UMAP
  flog.info(paste0("Calculating UMAP for genes"), name="DSP_NGS_log")
  map <- umap(as.matrix(t(norm_counts)))


  #add umap info to samp_notes for ease of plotting
  info <- data.frame(map$layout)
  info$cluster <- pca$call$X$clust
  info$gene <- rownames(info)
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

  umap <- ggplot(data = info, aes(x = X1, y = X2, colour = str_wrap(color, wrap_num),
                                  text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                               "\n", "Gene Group",": ", color, sep = ""))) +
    geom_point(size = info$size) +
    geom_point(data = subset(info, color != "Not Specified" & color != "All Probes"), size =  info$size[info$color!="Not Specified" & info$color != "All Probes"]) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "UMAP1", y = "UMAP2", title = "UMAP Genes", colour = "Gene Group") +
    scale_color_manual(values = clrs)+
    theme(legend.position="right", legend.direction = "vertical")

  ggsave(paste0(outdir, "/UMAP_genes.", fileType), umap, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot UMAP
  print(ggplotly(umap, tooltip = "text"))

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

    umap_plot <- ggplot(data = info, aes(x = X1, y = X2, colour = group_color,
                                         text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                                      "\n", grouping_var,": ", group_color, sep = ""))) +
      geom_point(size = info$size) +
      geom_point(data = subset(info, group_color != "Not Specified"), size = info$size[info$group_color!="Not Specified"]) +
      scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
      labs(x = "UMAP1", y = "UMAP2", title = "UMAP Genes", colour = grouping_var) +
      scale_color_manual(values = clrs)

    ggsave(paste0(outdir, "/UMAP_genes_", grouping_var, ".", fileType), umap_plot, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

    #plot info
    print(ggplotly(umap_plot, tooltip = "text"))
  }

  return(info)
}
