#' @title run_cell_decon
#'
#' Runs our internal immune cell deconvolution algorithm, "InSituSort"
#'
#' @param df dfs data frame containing azorius data
#' @param norm_method Raw data matric
#' @param outdir Directory to output figures
#' @param tumor_high_ids Names of AOIs with nearly 100% tumor content. Must be elements of the colnames of norm_counts.
#'   If NULL, the decon model will ignore tumo-intrinsic expression, which usually works just fine.
#' @param x Horizontal coordinates on which to plot decon results  (e.g. from t-SNE) (optional)
#' @param y Vertical coordinates on which to plot decon results (e.g. from t-SNE) (optional)
#' @param xlab Name of the x-dimension on which results are plotted
#' @param ylab Name of the y-dimension on which results are plotted
#' @param segment Column with segmentation information
#' @param TME label for TME segmentations
#' @return clustering information
#'
#' @examples
#'
#' deconresults <- run_cell_decon(df = dfs, matrix = "Human_Cell_Landscape", norm_method = "Q3",
#'                                outdir = cellType_dir, tumor_high_ids = NULL, x = c(1,2,3,4),
#'                                y = c(4,3,2,1), xlab = "X coordinates", ylab = "Y coordinates",
#'                                segment = "aoi_type", skip = "Tumor")
#' @export run_cell_decon
run_cell_decon <- function(df = dfs, matrix = NULL, norm_method = norm_method, outdir = cellType_dir,
                           tumor_high_ids = NULL, x = NULL, y = NULL, xlab = "", ylab = "",
                           segment = segmentation, skip = NULL){

  counts <- df[[1]]
  norm_counts <- df[[2]]
  samp_notes <- df[[3]]
  probe_notes <- df[[5]]

  probe_notes$Pooling <- as.integer(probe_notes$Pooling)

  indx <- which(!counts$TargetName %in% norm_counts$Gene)
  if (length(indx) > 0){
    counts <- counts[-indx,]
  }

  tumor_high_ids <- str_replace_all(string = tumor_high_ids, pattern = "[- ]", replacement = ".")

  #### prepare data for decon ------------------------------

  # format raw and normalized counts
  raw <- as.matrix(counts[, -1])
  rownames(raw) = counts[, 1]
  norm <- as.matrix(norm_counts[, -1])
  rownames(norm) = norm_counts[, 1]

  #Payman
  rem_gene_indx <- which(!probe_notes$TargetName %in% norm_counts$Gene)
  if (length(rem_gene_indx) > 0) {
    probepool <- probe_notes$Pooling[-rem_gene_indx]
  } else {
    probepool <- probe_notes$Pooling
  }

  ### calculate estimated background on scale of normalized data
  bg = derive_GeoMx_background(norm = norm,
                               probepool = probepool,
                               negnames = "Neg Probe") #"NegProbe-WTX"
  bg = bg[rownames(norm), colnames(norm)]



  if(exists("user_profile_matrix")){
    # read in user defined cell matrix and binned cell groups
    profile_matrix <- read.table(file = user_profile_matrix, header = T,
                                 sep = ",", row.names = 1)
    matrix <- "custom"

    if(!is.null(binned_cellTypes)){
      cellGroups <- read.table(binned_cellTypes, fill = TRUE,
                               stringsAsFactors = FALSE)

      group_names <- cellGroups[,1]
      cellGroups <- cellGroups[,-1]

      cellGroups <- split(cellGroups, seq_len(nrow(cellGroups)))
      cellGroups <- lapply(cellGroups, function(x) x[x != ""])
      names(cellGroups) <- group_names

      if(sum(lengths(cellGroups)) > ncol(profile_matrix) | !all(colnames(profile_matrix) %in% cellGroups)){
        if(!all(colnames(profile_matrix) %in% cellGroups)){
          missingCT <- colnames(profile_matrix)[which(!colnames(profile_matrix) %in% cellGroups)]
          print(paste("WARNING:", paste(missingCT, collapse = ", "), "are not included in given cell type bins. Continuing with no cell type bins"))
        }else{
          print("WARNING: binned_cellTypes has more cell types than the given profile matrix, check if spaces in cell types are changed to . like cell.type.1. Continuing with no cell type bins")
        }

        cellGroups <- as.list(colnames(profile_matrix))
        names(cellGroups) <- colnames(profile_matrix)
      }
    }else{
      cellGroups <- as.list(colnames(profile_matrix))
      names(cellGroups) <- colnames(profile_matrix)
    }

    sink(paste0(outdir, "custom_cellGroups.txt"))
    print(cellGroups)
    sink()

    profile_matrix <- as.matrix(profile_matrix)
  }else if (!is.null(matrix)){
    ### download cell profile matrix if not default NanoString Immune
    profile_matrix <- download_profile_matrix(matrixname = matrix)

    cellGroups <- list()
    done <- NULL
    notCellTypes <- c("stratified", "fetal", "progenitor", "progenitors", "high", "cell",
                      "cells", ".", "call", "calls", "decidual", "glandular", "alveolar", "like")
    cellTypes <- c("doublet","CD34", "muscle", "myofibroblast", "pancreatic", "phagocyte",
                   names(safeTME.matches), "T", "antigen.presenting", "not.available", "epithelial",
                   "interneuron", "neuron", "ureteric.bud", "loop.of.henle", "plasmocyte", "cone", "rod")
    for(i in cellTypes){
      if(i == "NK"){
        group <- "natural\\.killer"
      }else if(i == "pDC"){
        group <- "plasmacytoid\\.dendritic"
      }else if(i == "mDCs"){
        group <- "conventional\\.dendritic\\."
      }else if(i == "Treg"){
        group <- "regulatory\\.T"
      }else if(i == "CD4.T.cells"){
        group <- "cd4.*t\\."
      }else if(i == "CD8.T.cells"){
        group <- "cd8.*t\\."
      }else if(i == "T"){
        group <- c("^(?:(?!cd).)*\\.t\\.", "^t\\..*")
      }else if(i == "doublet"){
        group <- "doublet"
      }else if(i == "not.available"){
        group <- "unknown"
      }else if(i == "endothelial.cells"){
        group <- "endothelial"
      }else if(i == "plasma"){
        group <- "plasma\\."
      }else if(i == "pancreatic"){
        group <- c("\\..\\.pancreatic", "pancreatic\\..\\.")
      }else if(i == "B"){
        group <- c("\\.b\\.", "^b\\..*")
      }else{
        if(endsWith(i, "s")){i <- substr(i,1,nchar(i)-1)}
        group <- c(tolower(paste0(gsub(pattern = "\\.", replacement = "\\\\.", i),"\\.")),tolower(gsub(pattern = "\\.", replacement = "\\\\.", i)))
      }

      w2kp <- NULL
      for(g in group){
        w2kp <- c(w2kp, grep(perl = T, tolower(colnames(profile_matrix)), pattern = g, fixed = F))
      }
      double <- unique(w2kp[which(w2kp %in% done)])
      w2kp <- unique(w2kp[which(!w2kp %in% done)])
      done <- c(done, w2kp)

      if(length(w2kp) > 0){
        if(i %in% names(cellGroups)){
          cellGroups[[i]] <- c(cellGroups[[i]], colnames(profile_matrix)[w2kp])
        }else{
          cellGroups[[i]] <-  colnames(profile_matrix)[w2kp]
        }
      }
    }


    w2kp <- which(grepl(tolower(colnames(profile_matrix)), pattern = "[0-9]$") | grepl(tolower(colnames(profile_matrix)), pattern = "type\\.") |
                    grepl(tolower(colnames(profile_matrix)), pattern = "high") | grepl(tolower(colnames(profile_matrix)), pattern = "intermediated"))
    w2kp <- w2kp[which(!w2kp %in% done)]
    if(length(w2kp) > 0){
      done <- c(done, w2kp)
      for(i in w2kp){
        group <- tolower(unlist(strsplit(colnames(profile_matrix)[i], ".", fixed = T)))
        if(group[1] %in% notCellTypes){
          group <- group[2]
        }else{
          group <- group[1]
        }

        if(endsWith(group, "s") & group != "mucous"){group <- substr(group,1,nchar(group)-1)}

        if(group %in% names(cellGroups)){
          cellGroups[[group]] <- c(cellGroups[[group]], colnames(profile_matrix)[[i]])
        }else{
          cellGroups[[group]] <-  colnames(profile_matrix)[[i]]
        }
      }
    }

    w2kp <- NULL
    for(i in notCellTypes){
      w2kp <- c(w2kp, which(endsWith(tolower(colnames(profile_matrix)), suffix = i)))
    }
    w2kp <- w2kp[which(!w2kp %in% done)]
    done <- c(done, w2kp)
    if(length(w2kp) > 0){
      for(i in w2kp){
        group <- tolower(unlist(strsplit(colnames(profile_matrix)[i], ".", fixed = T)))
        if(group[length(group)-1] %in% notCellTypes){
          group <- group[length(group)-2]
        }else{
          group <- group[length(group)-1]
        }
        if(endsWith(group, "s") & group != "mucous"){group <- substr(group,1,nchar(group)-1)}

        if(group %in% names(cellGroups)){
          cellGroups[[group]] <- c(cellGroups[[group]], colnames(profile_matrix)[[i]])
        }else{
          cellGroups[[group]] <-  colnames(profile_matrix)[[i]]
        }
      }
    }

    w2kp <- which(!1:length(colnames(profile_matrix)) %in% done)
    done <- c(done, w2kp)
    if(length(w2kp) > 0){
      for(i in w2kp){
        group <- unlist(strsplit(colnames(profile_matrix)[i], ".", fixed = T))
        group <- tolower(group[length(group)])
        if(endsWith(group, "s") & group != "mucous"){group <- substr(group,1,nchar(group)-1)}

        if(group %in% names(cellGroups)){
          cellGroups[[group]] <- c(cellGroups[[group]], colnames(profile_matrix)[[i]])
        }else{
          cellGroups[[group]] <-  colnames(profile_matrix)[[i]]
        }
      }
    }

    sum <- 0
    for(i in cellGroups){sum <- sum + as.numeric(length(cellGroups[i]))}

    sink(paste0(outdir, matrix, "_cellGroups.txt"))
    print(cellGroups)
    sink()

    shared_genes <- intersect(rownames(norm), rownames(profile_matrix))
    if(length(shared_genes) < 5){
      stop("Chosen profile matrix is intended for mouse, choose another matrix. In later release these genes will be replaced with human orthologs")
    }

    profile_matrix <- as.matrix(profile_matrix)
  }else{
    cellGroups <- NULL
    w2kp <- NULL
    matrix <- "NanoStringImmune"
    sink(paste0(outdir, matrix, "_cellGroups.txt"))
    print(safeTME.matches)
    sink()
    profile_matrix <- NULL
  }

  if(!is.null(cell_count)){
    cell_counts <- samp_notes[[cell_count]]
  }else{
    cell_counts <- NULL
  }


  #### run main decon function: ----------------------------------
  ## res <- insitusortTILs(norm = norm,
  ##                      raw = raw,
  ##                      bg = bg,
  ##                      X = profile_matrix,
  ##                      is_pure_tumor = is.element(colnames(norm), high_tumor_AOIs),
  ##                      cell_counts = as.numeric(cell_counts),
  ##                      matchmajor = cellGroups)
  res <- spatialdecon(norm = norm,
                      bg = bg,
                      raw = raw,
                      X = profile_matrix,
                      is_pure_tumor = is.element(colnames(norm), high_tumor_AOIs),
                      cell_counts = as.numeric(cell_counts) )

  #### save decon results: ----------------------------------
  save(res, file = paste0(outdir, "cell_decon_results_", matrix, ".RData"))
  write.csv(res$beta, file = paste0(outdir, "cell_decon_abundance_estimates--", matrix, ".csv"))
  write.csv(res$beta.granular, file = paste0(outdir, "cell_decon_abundance_estimates_granular--", matrix, ".csv"))
  if(!is.null(cell_counts)){
    write.csv(res$cell.counts$cell.counts, file = paste0(outdir, "cell_decon_cell_counts--", matrix, ".csv"))
    write.csv(res$cell.counts.granular$cell.counts, file = paste0(outdir, "cell_decon_cell_counts_granular--", matrix, ".csv"))
  }
  #### plot decon results: ----------------------------------

  barplots = function(mat, draw_legend = FALSE, colors = NULL,...) {
    ret <- FALSE
    usecells = rownames(mat)[which(rowSums(mat) > 0)]

    if(is.null(colors)){
      if(matrix == "NanoStringImmune"){
        colors <- cellcols[which(names(cellcols) %in% usecells)]
      }else{
        manycols = c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69',
                     '#FCCDE5','#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C',
                     '#FDBF6F','#FF7F00','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E',
                     '#E6AB02','#A6761D','#666666', sample(grDevices::colors(), 150))
        colors <- manycols[1:length(usecells)]
        names(colors) <- usecells
        w2kp <- match(names(colors), names(cellcols), nomatch = 0)
        if(sum(w2kp>0) > 0){
          colors[1:sum(w2kp>0)] <- cellcols[w2kp]
        }
      }
      ret <- TRUE
    }


    # draw barplot:
    graphics::barplot(mat[usecells, ], cex.lab = 1.5, main = "", border = NA,
                      las = 2, col = colors)

    # draw a legend:
    if (draw_legend) {
      graphics::frame()
      graphics::legend("center",
                       fill = rev(colors),
                       legend = rev(usecells))
    }

    if(ret){return(colors)}
  }


  # stacked plot: barplots of abundances and proportions
  pdf(paste0(outdir, "cell_abundance_barplots--", matrix, ".pdf"), width = 12, height = 10)
  graphics::layout(mat = matrix(c(1,3,2,2), 2), widths = c(11, 3, 11, 3), heights = c(6.5, 3.5, 6.5, 3.5))
  par(mar = c(0.1, 5.5, 0, 0))
  colors <- barplots(mat = res$beta, draw_legend = TRUE, names.arg = rep("", ncol(res$beta)))
  par(mar = c(10, 5.5, 0.1, 0))
  w2kp <- which(apply(res$prop_of_all, 2, function(x) any(!is.na(x))))
  barplots(mat = res$prop_of_all[, w2kp], cex.names = 0.5, colors = colors)
  dev.off()

  if(!is.null(cell_counts)){
    # stacked plot: barplots of cell counts
    pdf(paste0(outdir, "cell_counts_barplots--", matrix, ".pdf"), width = 12, height = 10)
    graphics::layout(mat = matrix(c(1,3,2,2), 2), widths = c(11, 3, 11, 3), heights = c(9.5, 3.5, 6.5, 0.5))
    par(mar = c(0.1, 5.5, 0, 0))
    barplots(mat = res$cell.counts$cell.counts, draw_legend = TRUE, names.arg = rep("", ncol(res$cell.counts$cell.counts)))
    dev.off()
  }



  # heatmaps of results
  if(sum(apply(res$prop_of_all, 2, function(x) any(is.na(x)))) < ncol(res$prop_of_all)*.1){

    # plot decon results in x-y space:
    # *** hydra team: it could be cool to replace the below with plotly output, but not required *** <-------------------
    if ((length(x) > 0) & (length(y) > 0)){
      if(!is.null(segment) & xy_axes == "xy"){
        w2kp <- NULL
        for(i in 1:length(x)){
          temp <- which((x < (x[i] + x[i]*0.0005) & x > x[i] - x[i]*0.0005) & (y < y[i] + y[i]*0.0005 & y > y[i] - y[i]*0.0005))
          w2kp <- c(w2kp, temp[which(temp != i)])
        }
        w2kp <- w2kp[which(samp_notes[w2kp,segment] %in% skip)]
      } else {
        w2kp <- NULL
      }

      pdf(paste0(outdir, "cell_abundance_in_", xy_axes, "_space--", matrix, ".pdf"), width = 12, height = 9)
      graphics::layout(mat = matrix(c(1, 2), 1), widths = c(9, 3), heights = c(9, 9))
      if(length(w2kp) > 0){
        print(paste("Creating florets in space figures without", paste(skip, collapse = ", "), "segments"))
        florets(x = x[-w2kp] / diff(range(x[-w2kp])),
                y = y[-w2kp] / diff(range(y[-w2kp])),
                b = res$beta[names(colors),-w2kp] / quantile(res$beta[names(colors),-w2kp], 0.95) * 0.001,
                legendwindow = FALSE,
                xaxt = "n", yaxt = "n",
                xlab = xlab, ylab = ylab, col = colors)
      }else{
        florets(x = x / diff(range(x)),
                y = y / diff(range(y)),
                b = res$beta[names(colors),] / quantile(res$beta[names(colors),], 0.95) * 0.001,
                legendwindow = FALSE,
                xaxt = "n", yaxt = "n",
                xlab = xlab, ylab = ylab, col = colors)
      }
      par(mar = c(0,0,0,0))
      frame()
      legend("center", fill = rev(colors), legend = rev(names(colors)), cex = 0.8)
      dev.off()
    }

    w2kp <- which(apply(res$prop_of_all, 2, function(x) any(!is.na(x))))
    res$beta <- res$beta[,w2kp]
    res$prop_of_all <- res$prop_of_all[,w2kp]

    frame();dev.off()
    pdf(paste0(outdir, "cell_abundance_heatmaps--", matrix, ".pdf"), width = 12, height = 10)
    pheatmap(t(res$beta), col = viridis_pal(option = "C", begin = 0, end = 1)(100),
             show_rownames = T, main = "abundance")
    pheatmap(t(res$beta), col = colorRampPalette(c("white", "dodgerblue4"))(100),
             show_rownames = T, main = "abundance")
    pheatmap(t(res$prop_of_all), col = viridis_pal(option = "C", begin = 0, end = 1)(100),
             show_rownames = T, main = "Proportion")
    pheatmap(t(res$prop_of_all), col = colorRampPalette(c("white", "dodgerblue4"))(100),
             show_rownames = T, main = "Proportion")
    if(!is.null(cell_counts)){
      pheatmap(t(res$cell.counts$cell.counts), col = viridis_pal(option = "C", begin = 0, end = 1)(100),
               show_rownames = T, main = "Cell Counts")
      pheatmap(t(res$cell.counts$cell.counts), col = colorRampPalette(c("white", "dodgerblue4"))(100),
               show_rownames = T, main = "Cell Counts")
    }
    dev.off()
  }
  return(res)
}
