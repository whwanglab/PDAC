#' @title run_PCA
#'
#' Runs principle component analysis (PCA) on dsp data
#'
#' @param df read in dsp excel workbook dataset
#' @param outdir output directory for figures
#' @param outlier_cutoff multiplier for outlier cutoff
#' @param color category to color the graph by
#' @param symbol category to change symbols in the graph
#' @param size should the size of the points change by gene_count
#'
#' @return clustering information
#'
#' @examples
#' pca <- run_PCA(df, "file/path", outlier_cutoff)
#'
#' @export run_PCA
#' @export run_PCA_genes

run_PCA <- function(df, outdir, color = "cluster", symbol = "dsp_slide", size = TRUE, colors = colors){

  #read in counts without outliers and notes with slide name split
  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]

  if ((!any(names(samp_notes) == color) & color != "cluster") |
      (!any(names(samp_notes) == symbol) & symbol != "cluster")){
    stop(str_wrap(paste("Given coloring or symbol factor was not allowed. OPTIONS: ",
                        paste(names(samp_notes), collapse = ", ")), width = 80), call.=FALSE)
  }

  if (class(samp_notes$Sample_ID) == "integer"){
    samp_notes$Sample_ID <- paste0("X", samp_notes$Sample_ID)
  }

  flog.info(paste0("Running PCA"), name="DSP_NGS_log")
  # Compute PCA with ncp = 5
  res.pca <- PCA(t(norm_counts), ncp = 5, graph = FALSE, scale.unit = TRUE)
  # Compute hierarchical clustering on principal components
  res.hcpc <- HCPC(res.pca, graph = FALSE)

  #reorder to match samp_notes
  res.hcpc$call$X <- res.hcpc$call$X[match(samp_notes$Sample_ID, row.names(res.hcpc$call$X)),]

  norm_counts <- norm_counts[,match(samp_notes$Sample_ID, colnames(norm_counts))]

  # eignenvalue bar graph
  outfile = paste0(outdir, "/Eigenvalues.", fileType)
  res.pca$eig <- as.data.frame(res.pca$eig)

  match.fun(fileType)(outfile)
  barplot(res.pca$eig$`percentage of variance`[1:10], ylim=c(0,100), xlab = "Component",
          ylab = "Percent of Variance", main = "Eigenvalues", names.arg = 1:10)
  dev.off()

  outfile = paste(outdir, "AOIvsComponents.csv", sep="/")
  write.table(res.pca$ind$contrib, file = outfile, sep = ",",col.names = NA)


  #create and save hierarchical clustering dendrogram
  outfile = paste0(outdir, "/Cluster_dendrogram.", fileType)

  dend <- fviz_dend(res.hcpc,
            cex = 0.7,                     # Label size
            palette = "jco",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "jco",           # Rectangle color
            labels_track_height = max(nchar(samp_notes$Sample_ID) * 10),   # Augment the room for labels
            horiz = T)

  match.fun(fileType)(outfile, width = max(nrow(samp_notes)*6, 250), height = max((nrow(samp_notes))*25, 500))
  print(dend)
  dev.off()

  #add pca info to samp_notes for ease of plotting
  samp_notes$pca  <- data.frame(res.hcpc$call$X[,1:2])
  samp_notes$gene_count <- colSums(norm_counts)
  samp_notes$cluster <- res.hcpc$call$X$clust

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
  }else if (length(unique(symbol_by)) > 14){
    stop("Too many shapes requested. Use this column as color_by instead.")
  }
  if (size == TRUE){
    size = samp_notes$gene_count / 50000
  }else {
    size = 6
  }

  clrs <- unlist(colors[color], use.names = FALSE)
  names(clrs) <- names(colors[[color]])

  pca <- ggplot(data = samp_notes, aes(x = pca$Dim.1, y = pca$Dim.2, colour = str_wrap(color_by, wrap_num),
                                       text = paste("ROI: ", Sample_ID,"\nGene Count: ", round(gene_count),
                                                    "\nColor: ", color_by, "\nSymbol:", symbol_by))) +
    geom_point(size = size, aes(shape = symbol_by)) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "PCA1", y = "PCA2", title = "PCA", colour = color, shape = symbol) +
    scale_color_manual(values = clrs)

  ggsave(paste0(outdir, "/PCA.", fileType), pca, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot PCA
  print(ggplotly(pca, tooltip = "text"))


  #return clustering
  return(res.hcpc)
}

run_PCA_genes <- function(df, outdir, size = TRUE, de = de_results, expressed_genes = de_genes){

  samp_notes <- df[[3]]
  norm_counts <- df[[2]][, !colnames(df[[2]]) %in% c("Gene","TargetName")]
  norm_counts <- norm_counts[,samp_notes$Sample_ID %in% colnames(norm_counts)]
  norm_counts <- t(norm_counts)

  flog.info(paste0("Running PCA on genes"), name="DSP_NGS_log")
  # Compute PCA with ncp = 5
  res.pca <- PCA(t(norm_counts), ncp = 5, graph = FALSE, scale.unit = TRUE)
  # Compute hierarchical clustering on principal components
  res.hcpc <- HCPC(res.pca, graph = FALSE)

  norm_counts <- norm_counts[,match(row.names(res.hcpc$call$X), colnames(norm_counts))]

  # eignenvalue bar graph
  outfile = paste0(outdir, "/Eigenvalues.", fileType)
  res.pca$eig <- as.data.frame(res.pca$eig)

  match.fun(fileType)(outfile)
  barplot(res.pca$eig$`percentage of variance`[1:10], ylim=c(0,100), xlab = "Component",
          ylab = "Percent of Variance", main = "Eigenvalues", names.arg = 1:10)
  dev.off()

  outfile = paste(outdir, "GenevsComponents.csv", sep="/")
  write.table(res.pca$ind$contrib, file = outfile, sep = ",",col.names = NA)

  #add pca info to samp_notes for ease of plotting
  pca <- data.frame(res.hcpc$call$X[,1:2])
  pca$cluster <- res.hcpc$call$X$clust
  pca$gene <- rownames(pca)
  pca$gene_count <- colSums(norm_counts)
  pca$color <- "Not Specified"

  for (i in 1:length(gene_group)){
    w2kp <- which(pca$gene %in% unlist(genes[i]))
    w2kp2 <- which(pca$color[w2kp] != "Not Specified")

    pca$color[w2kp] <- names(genes[i])
    pca$color[w2kp2] <- "Multiple"
  }

  w2kp <- which(pca$gene %in% fav_genes)
  pca$color[w2kp] <- "Fav.Genes"

  #check options for function
  if (size == TRUE){
    sizes <- pca$gene_count / 30000
    w2kp <- which(sizes < 0.2)
    sizes[w2kp] <- 0.2
  }else {
    sizes <- ifelse(pca$color == "Not Specified", 0.8, 2)
  }

  pca$size <- sizes

  clrs <- unlist(colors["gene_groups"], use.names = FALSE)
  names(clrs) <- names(colors[["gene_groups"]])

  pca_plot <- ggplot(data = pca, aes(x = Dim.1, y = Dim.2, colour = str_wrap(color, wrap_num),
                                     text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                                  "\n", "Gene Group",": ", color, sep = ""))) +
    geom_point(size = pca$size) +
    geom_point(data = subset(pca, color != "Not Specified" & color != "All Probes"), size = pca$size[pca$color!="Not Specified" & pca$color != "All Probes"]) +
    scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
    labs(x = "PCA1", y = "PCA2", title = "PCA Genes", colour = "Gene Group") +
    scale_color_manual(values = clrs)+
    theme(legend.position="right", legend.direction = "vertical")

  ggsave(paste0(outdir, "/PCA_genes.", fileType), pca_plot, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

  #plot PCA
  print(ggplotly(pca_plot, tooltip = "text"))

  zero <- ifelse(test = sum(lengths(expressed_genes)) > 0, yes = FALSE, no = TRUE)

  pca$group_color <- "Not Specified"

  if(zero == FALSE){
    for (i in 1:length(expressed_genes)){
      w2kp <- which(pca$gene %in% unlist(expressed_genes[i]))
      w2kp2 <- which(pca$group_color[w2kp] != "Not Specified")

      pca$group_color[w2kp] <- names(expressed_genes[i])
      pca$group_color[w2kp2] <- "Multiple"
    }


    w2kp <- which(pca$group_color == "Multiple")
    pca$group_color[w2kp] <- "Not Specified"

    if (size == TRUE){
      sizes <- pca$gene_count / 40000
      w2kp <- which(sizes < 0.15)
      sizes[w2kp] <- 0.15
    }else {
      sizes <- ifelse(pca$group_color == "Not Specified", 0.8, 2)
    }

    pca$size <- sizes

    clrs <- unlist(colors[grouping_var], use.names = FALSE)
    names(clrs) <- names(colors[[grouping_var]])

    clrs[['Not Specified']] <- '#ADADAD'

    pca_plot <- ggplot(data = pca, aes(x = Dim.1, y = Dim.2, colour = group_color,
                                       text = paste("Gene: ", gene, "\nGene Count: ", round(gene_count),
                                                    "\n", grouping_var,": ", group_color, sep = ""))) +
      geom_point(size = pca$size) +
      geom_point(data = subset(pca, group_color != "Not Specified"), size = pca$size[pca$group_color!="Not Specified"]) +
      scale_shape_manual(values = c(16,17,15,18,22,23,25,5,3,9,8,2,1,7)) +
      labs(x = "PCA1", y = "PCA2", title = "PCA Genes", colour = grouping_var, shape = "Cluster") +
      scale_color_manual(values = clrs)

    ggsave(paste0(outdir, "/PCA_genes_", grouping_var, ".", fileType), pca_plot, device = fileType, width = 7.5, height = 5, units = 'in', scale = 1.5)

    #plot PCA
    print(ggplotly(pca_plot, tooltip = "text"))
  }else{
    print("No figure generated due to no differentially expressed genes, lower the fc_cutoff and/or pval_cutoff to create figure")
  }

  #return clustering
  return(res.hcpc)
}
