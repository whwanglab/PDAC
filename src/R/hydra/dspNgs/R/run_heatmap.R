#' @title run_heatmap
#'
#' Generic call to heatmap functionality
#'
#' @param data the data to plot, assumed targets * samples
#' @param annot sample annotations
#' @param width file width
#' @param height file height
#' @param samp_vars annotations to include on heatmap
#' @param fav_genes favorite genes to highlight
#' @param ann_colors list of lists of colors
#' @param logData T/F - log2 transform data
#' @param targets used to subset dataset to only genes of interest or NULL for all
#' @param name title of graph and file
#' @param scale passed to pheatmap
#' @param outdir where to print file
#' @param fileType function to call for file printing
#' @param sortBy annotation variable to use for sorting the heatmap
#'
#' @return
#'
#' @examples
#'
#'
#' @export run_heatmap

run_heatmap <- function(data = NULL,
                        data_TargetID = 'Gene',
                        annot = NULL,
                        samp_vars = NULL,
                        fav_genes = NULL,
                        ann_colors = NULL,
                        logData = TRUE,
                        targets = NULL,
                        name = NULL,
                        scale = 'row',
                        outdir = '.',
                        fileType = 'tiff',
                        sortBy = NULL,
                        sortBaseline = NULL,
                        ...) {
  # add sorting variable to samp_vars
  if(!is.null(sortBy)){
    if (!sortBy %in% samp_vars){
      samp_vars <- c(samp_vars, sortBy)
    }
  }

  # factorize sorting
  sortOrder <- NULL
  if(!is.null(sortBaseline)) {
    lvls <- levels(as.factor(annot[[sortBy]]))
    sortByList <- factor(annot[[sortBy]],
                         levels = c(sortBaseline,
                                    lvls[lvls != sortBaseline]))
    sortOrder <- order(sortByList)
    rm(sortByList)
  } #else {
    #sortOrder <- order(sortByList)
  #}

  # subset to specific annotations
  if(!is.null(samp_vars)) {
    annot <- subset(annot, select = c("Sample_ID", samp_vars))
    rownames(annot) <- annot[,1]
    annot <- annot[,-1]
  }

  # subset targets to targets in gene column
  targets <- targets[targets %in% data[, data_TargetID]]

  # pair down and log transform data
  rowlabs <- data[, data_TargetID]
  rownames(data) <- make.names(data[, data_TargetID])
  data <- data[, colnames(data) != data_TargetID]
  if(is.null(targets)) {
    targets <- rownames(data)
  } else {
    rowlabs <- targets
    targets <- make.names(targets)
  }
  targets <- targets[targets %in% rownames(data)]
  data <- data[targets,]
  if(logData) {
    pre_trans_cols <- colnames(data)
    data <- data.frame(log2(data))
    colnames(data) <- pre_trans_cols
  }

  # favorite genes
  gene_ann <- NULL
  if(!is.null(fav_genes)) {
    gene_ann <- data.frame(favGenes = ifelse(rowlabs %in% fav_genes, 'True','False'),
                           row.names = row.names(data))
  }

  # set width / height
  doc_scale = ifelse(fileType %in% c('svg','pdf'), 1, 15)
  width <- ncol(data) * 1.6 * doc_scale
  if(ncol(data) > 50) {
    width <- 1000
  }else if (ncol(data) < 15){
    width <- width * 2.5
  }
  height <- max(nrow(data) * 1.3 * doc_scale,
                ifelse(fileType %in% c('svg','pdf'), 10, 500))
  if(ncol(data) > 63) {
    height <- 800
  }

  # sort graphs
  if(!is.null(sortBy)) {
    data <- data[, sortOrder]
    annot <- annot[sortOrder, ]
  }

  # plot pheatmap
  ph <- pheatmap(data,
                 cluster_cols = is.null(sortBy),
                 fontsize = 10,
                 labels_row = rowlabs,
                 cellheight = ifelse(nrow(data) < 63, 11, NA),
                 cellwidth = ifelse(ncol(data) < 51, 11, NA),
                 border_color = NA,
                 show_colnames = ncol(data) < 51,
                 show_rownames = nrow(data) < 63,
                 annotation_col = annot,
                 annotation_row = gene_ann,
                 annotation_colors = c(ann_colors,
                                       list(favGenes = c(True = 'green3',
                                                         False = 'white'))),
                 scale = scale,
                 main = name)

  # save heatmap
  match.fun(fileType)(paste0(outdir,
                             paste0("/",
                                    name, " Heatmap."),
                             fileType),
                      width = width,
                      height = height)
  print(ph)
  dev.off()
}
