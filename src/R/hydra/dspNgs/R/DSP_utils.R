#' @title DSP_utils
#'
#' Collection of DSP_hydra functions
#'
#' @return
#'
#' @examples
#'
#' @export read_config
#' @export log_stop
#' @export create_dir
#' @export create_dirs
#' @export make_name_valid
#' @export genes_of_interest
#' @export set_colors
#' @export parse_GeneSets
#' @export write_dim_reduct_files
#' @export ngeoMean
#' @export nanoNormFactors

read_config <- function(){
  config <- config::get(config='user_options')

  # general
  assign("species", tolower(config$species), .GlobalEnv)
  speciesAllowed <- c('mouse', "human")
  if(!species %in% speciesAllowed){
    stop(paste("ERROR:", species, "not valid"))
  }

  # output and inputs directories
  assign("outdir", config$output_dir, .GlobalEnv)
  assign("outname", config$output_name, .GlobalEnv)
  assign("data_path", config$data_path, .GlobalEnv)
  assign("txt", config$txt, .GlobalEnv)
  assign("norm_method", config$norm_method, .GlobalEnv)
  assign("fileType", config$fileType, .GlobalEnv)
  assign("DA_raw", config$DA_raw, .GlobalEnv)
  assign("DA_norm", config$DA_norm, .GlobalEnv)

  ftOptions <- c("png", "jpeg", "pdf", "tiff", "bmp", "svg")
  if(!fileType %in% ftOptions){
    print("WARNING: image file type not valid, continuing with png")
    assign("fileType", "png", .GlobalEnv)
  }

  assign("RData", config$RData_path, .GlobalEnv)

  if(is.null(config$rm_genes_list)){
    if(is.null(config$rm_genes_txt)){
      assign("rm_genes", NULL, .GlobalEnv)
    }else{
      genes <- as.character(read.table(config$rm_genes_txt, header = FALSE)$V1)
      assign("rm_genes", genes, .GlobalEnv)
    }
  }else if(is.null(config$rm_genes_list)){

  }else{
    assign("rm_genes", eval(parse(text=config$rm_genes_list)), .GlobalEnv)
  }

  # cutoffs
  assign("thresh", config$thresh, .GlobalEnv)
  assign("pval_cutoff", config$pval_cutoff, .GlobalEnv)
  assign("fc_cutoff", config$fc_cutoff, .GlobalEnv)
  assign("fdr_cutoff", config$fdr_cutoff, .GlobalEnv)

  if(is.null(pval_cutoff) & is.null(fdr_cutoff)){
    stop(paste("ERROR: pval_cutoff and fdr_cutoff cannot both be NULL"))
  }

  assign("outlier_cutoff", config$outlier_cutoff, .GlobalEnv)
  assign("LOQ_cutoff", config$LOQ_cutoff, .GlobalEnv)
  assign("LOQ_level", config$LOQ_level, .GlobalEnv)
  assign("LOQ_floor", config$LOQ_floor, .GlobalEnv)
  assign("sat_cutoff", config$sat_cutoff, .GlobalEnv)
  assign("sat_factcrit", config$sat_factcrit, .GlobalEnv)
  assign("sat_numcrit", config$sat_numcrit, .GlobalEnv)

  # sample removal
  assign("remove_samples_col", config$remove_samples_col, .GlobalEnv)
  if(!is.null(remove_samples_col)){
    if(startsWith(remove_samples_col, "c(")){
      assign("remove_samples_col", eval(parse(text=remove_samples_col)), .GlobalEnv)
    }
  }
  assign("remove_samples_condition", config$remove_samples_condition, .GlobalEnv)
  if(!is.null(remove_samples_condition)){
    if(startsWith(remove_samples_condition, "c(")){
      assign("remove_samples_condition", eval(parse(text=remove_samples_condition)), .GlobalEnv)
    }
  }
  assign("remove_samples_operators", config$remove_samples_operators, .GlobalEnv)
  if(!is.null(remove_samples_operators)){
    if(startsWith(remove_samples_operators, "c(")){
      assign("remove_samples_operators", eval(parse(text=remove_samples_operators)), .GlobalEnv)
    }
  }
  operator_Options <- c("==", "<", "<=", ">", ">=", "!=")
  if(any(!remove_samples_operators %in% operator_Options)){
    w2kp <- which(!remove_samples_operators %in% operator_Options)
    stop(paste("ERROR:", remove_samples_operators[w2kp], "not valid"))
    # print(paste("WARNING:", remove_samples_operators[w2kp], "not valid; continuing with all =="))
    # remove_samples_operators[w2kp] <- "=="
    # assign("remove_samples_operators", remove_samples_operators, .GlobalEnv)
  }

  if(length(remove_samples_col) != length(remove_samples_operators) |
     length(remove_samples_col) != length(remove_samples_condition)){
    stop("remove samples parameters are not the same length")
  }

  assign("remove_samples_and", config$remove_samples_and, .GlobalEnv)
  assign("remove_samples_keep", config$remove_samples_keep, .GlobalEnv)
  assign("remove_samples_custom", config$remove_samples_custom, .GlobalEnv)

  # load dimensional reduction variables
  assign("perplexity", config$tsne_perplexity, .GlobalEnv)
  assign("color_by", config$color_by, .GlobalEnv)
  assign("shape_by", config$shape_by, .GlobalEnv)
  assign("size", config$size, .GlobalEnv)

  # de variables
  assign("hk_genes", config$hk_genes, .GlobalEnv)
  if(!is.null(hk_genes)){
    assign("hk_genes", eval(parse(text=hk_genes)), .GlobalEnv)
  }


  assign("grouping_var", config$grouping_var, .GlobalEnv)
  assign("n_top", config$n_top, .GlobalEnv)
  assign("n_top_heatmap", config$n_top_heatmap, .GlobalEnv)
  assign("base_level", config$base_level, .GlobalEnv)
  assign("control_var", config$control_var, .GlobalEnv)
  assign("de_results", config$de_results, .GlobalEnv)
  if(!is.null(de_results)){
    de_results <- read.delim(de_results, header = T, sep = ",")
    assign("de_results", de_results, .GlobalEnv)
  }
  assign("draw_DE_plots", config$draw_DE_plots, .GlobalEnv)

  assign("fav_genes", config$fav_genes, .GlobalEnv)
  if(!is.null(fav_genes)){
    if(startsWith(fav_genes, "c(")){
      assign("fav_genes", eval(parse(text=fav_genes)), .GlobalEnv)
    }
  }

  assign("gene_group", config$gene_group, .GlobalEnv)
  if(!is.null(gene_group)){
    if(startsWith(gene_group, "c(")){
    assign("gene_group", eval(parse(text=gene_group)), .GlobalEnv)
    }
  }

  assign("default_color", config$default_color, .GlobalEnv)
  assign("fc_color", config$fc_color, .GlobalEnv)
  assign("label_size", config$label_size, .GlobalEnv)
  assign("label_fc", config$label_fc, .GlobalEnv)

  assign("volcano_color", config$volcano_color, .GlobalEnv)
  if(!volcano_color %in% c("Significance", "Gene Group")){
    print("WARNING: given volcano color not valid, continuing with Significance")
    assign("volcano_color", "Significance", .GlobalEnv)
  }
  assign("volcano_label", config$volcano_label, .GlobalEnv)
  if(!volcano_label %in% c("Significance", "Fav Genes")){
    print("WARNING: given volcano label not valid, continuing with Significance")
    assign("volcano_label", "Significance", .GlobalEnv)
  }

  # ROI variables
  assign("ROI_plot_type", config$ROI_plot_type, .GlobalEnv)
  assign("save_ROI_plots", config$save_ROI_plots, .GlobalEnv)
  assign("ROI_ID", config$ROI_ID, .GlobalEnv)
  if(!is.null(ROI_ID)){
    if(startsWith(ROI_ID, "c(")){
      assign("ROI_ID", eval(parse(text=ROI_ID)), .GlobalEnv)
    }
  }
  assign("AOI_ID", config$AOI_ID, .GlobalEnv)
  assign("segment_Annotation", config$segment_Annotation, .GlobalEnv)
  assign("AOI_names", config$AOI_names, .GlobalEnv)
  if(!is.null(AOI_names)){
    if(startsWith(AOI_names, "c(")){
      assign("AOI_names", eval(parse(text=AOI_names)), .GlobalEnv)
    }
  }
  assign("AOI_colors", config$AOI_colors, .GlobalEnv)
  if(!is.null(AOI_colors)){
    if(startsWith(AOI_colors, "c(")){
      assign("AOI_colors", eval(parse(text=AOI_colors)), .GlobalEnv)
    }
  }
  assign("targets_ROI_plot", config$targets_ROI_plot, .GlobalEnv)
  if(!is.null(targets_ROI_plot)){
    if(startsWith(targets_ROI_plot, "c(")){
      assign("targets_ROI_plot", eval(parse(text=targets_ROI_plot)), .GlobalEnv)
    }
  }
  assign("target_thresh", config$target_thresh, .GlobalEnv)
  assign("target_transformation", config$target_transformation, .GlobalEnv)
  assign("cluster_n", config$cluster_n, .GlobalEnv)

  # cell typing variables
  assign("high_tumor_AOIs", config$high_tumor_AOIs, .GlobalEnv)
  if(!is.null(high_tumor_AOIs)){
    if(startsWith(high_tumor_AOIs, "c(")){
    assign("high_tumor_AOIs", eval(parse(text=high_tumor_AOIs)), .GlobalEnv)
    }
  }
  matrices <- c("NULL","Airway_Epithelium","Atlas_Adult_Retina_10x","Census_Adult_Immune_10x",
                "Census_Newborn_Blood_10x","Diff_Fetal_Neuron_SS2","FetalMaternal_Adult_Blood_10x",
                "FetalMaternal_Adult_Blood_SS2","FetalMaternal_Adult_Decidua_10x",
                "FetalMaternal_Adult_Decidua_SS2","FetalMaternal_Fetal_Placenta_10x",
                "Human_brain","Human_Cell_Landscape","IBD_Adult_Colon_10x",
                "Landscape_Adult_Liver_10x","Lung_plus_neutrophils","Mouse_Brain",
                "Profiling_Adult_BoneMarrow_10x","Reprogram_Embryo_Dendritic_10x",
                "Sensitivity_Adult_Esophagus_10x","Sensitivity_Adult_Lung_10x",
                "Sensitivity_Adult_Spleen_10x","Somatic_Adult_Pancreas_SS2",
                "SpatioTemporal_Adult_Kidney_10x","SpatioTemporal_Fetal_Kidney_10x",
                "Tcell_Adult_Blood_10x","Tcell_Adult_BoneMarrow_10x","Tcell_Adult_Lung_10x",
                "Tcell_Adult_LymphNode_10x")
  if(suppressWarnings(is.na(as.integer(config$profile_matrix)))){
    assign("user_profile_matrix", config$profile_matrix, .GlobalEnv)
    if(!file.exists(user_profile_matrix)){
      print("WARNING: custom profile matrix path not valid, continuing with NanoString Immune")
      assign("profile_matrix", NULL, .GlobalEnv)
      rm(user_profile_matrix, envir = .GlobalEnv)
    }
  }else if(!as.integer(config$profile_matrix) %in% 1:length(matrices)){
    print("WARNING: profile matrix not valid, continuing with NanoString Immune")
    assign("profile_matrix", NULL, .GlobalEnv)
  }else if(as.integer(config$profile_matrix) == 1){
    assign("profile_matrix", NULL, .GlobalEnv)
  }else{
    assign("profile_matrix", matrices[as.integer(config$profile_matrix)], .GlobalEnv)
  }

  assign("binned_cellTypes",config$binned_cellTypes, .GlobalEnv)
  assign("segmentation",config$segmentation, .GlobalEnv)
  assign("skipped_seg",config$skipped_seg, .GlobalEnv)
  assign("cell_count",config$cell_count, .GlobalEnv)
  assign("xy_axes",tolower(config$xy_axes), .GlobalEnv)
  if(!xy_axes %in% c("pca", "tsne", "umap", "xy")){
    print("WARNING: xy_axes not valid, continuing with tsne")
    assign("profile_matrix", "tsne", .GlobalEnv)
  }
  assign("x_position_col",config$x_position, .GlobalEnv)
  assign("y_position_col",config$y_position, .GlobalEnv)
  assign("x_label",config$x_label, .GlobalEnv)
  assign("y_label",config$y_label, .GlobalEnv)

  # pathway analysis variables
  assign("geneListrank", config$geneListrank, .GlobalEnv)
  assign("path_fc", config$path_fc, .GlobalEnv)
  assign("path_pval", config$path_pval, .GlobalEnv)
  assign("enrichment", config$enrichment, .GlobalEnv)
  assign("custom_pathways_path", config$custom_pathways_path, .GlobalEnv)
  assign("exclude_reactome", config$exclude_reactome, .GlobalEnv)

  # load color schema
  assign("color_genes", config$color_genes, .GlobalEnv)
  assign("color_samples", config$color_samples, .GlobalEnv)
  assign("color_grouping", config$color_grouping, .GlobalEnv)

  # load themes
  assign("preset_theme", config$preset_theme, .GlobalEnv)
  assign("theme_adjustments", config$theme_adjustments, .GlobalEnv)
}

# function to stop execution on error and log error message
log_stop <- function(error_message){
  stop(call=TRUE, geterrmessage())
  flog.error(error_message, name="DSP_NGS_log")
  flog.error(geterrmessage(), name="DSP_NGS_log")
}

# create directories
create_dir <- function (path, name = ""){
  if (!dir.exists(path)){
    dir.create(path)
    flog.info(paste(name,"directory created.", sep = " "), name="DSP_NGS_log")
  } else {
    flog.info(paste(name,"directory exists.", sep = " "), name="DSP_NGS_log")
  }
}

create_dirs <- function(outdir){
  create_dir(outdir, "Output")

  # check if RESULTS dir exists, if not create
  assign("results_dir", paste0(outdir, "results"), .GlobalEnv)
  create_dir(results_dir, "Results")

  # check if QC dir exists, if not create
  assign("qc_dir", paste(results_dir, "qc", sep="/"), .GlobalEnv)
  create_dir(qc_dir, "QC")

  # check if DE dir exists, if not create
  assign("de_dir", paste(results_dir, "de", sep = "/"), .GlobalEnv)
  create_dir(de_dir, "DE")

  # check if MODEL Fit dir exists, if not create
  assign("models_dir", paste(de_dir, "modelFit", sep = "/"), .GlobalEnv)
  create_dir(models_dir, "modelFit")

  # check if CLUSTERING dir exists, if not create
  assign("clustering_dir", paste(results_dir, "clustering", sep="/"), .GlobalEnv)
  create_dir(clustering_dir, "Clustering")

  # check if PATHWAY dir exists, if not create
  assign("path_dir", paste(results_dir, "pathway", sep="/"), .GlobalEnv)
  create_dir(path_dir, "Pathway")

  # check if SAMPLES dir exists, if not create
  assign("samples_dir", paste(clustering_dir, "samples", sep="/"), .GlobalEnv)
  create_dir(samples_dir, "Samples")

  # check if GENES dir exists, if not create
  assign("genes_dir", paste(clustering_dir, "genes", sep="/"), .GlobalEnv)
  create_dir(genes_dir, "Genes")

  # check if RData dir exists, if not create
  assign("RData_dir", paste0(outdir, "RData"), .GlobalEnv)
  create_dir(path = RData_dir, name = "RData")

  # check if spatialPlots dir exists, if not create
  assign("spatial_dir", paste(results_dir, "spatialPlots/", sep = "/"), .GlobalEnv)
  create_dir(path = spatial_dir, name = "spatialPlots")

  # check if cellType dir exists, if not create
  assign("cellType_dir", paste(results_dir, "cellTyping/", sep = "/"), .GlobalEnv)
  create_dir(path = cellType_dir, name = "cellTyping")
}

#' @title function to coerce DSP-DA sample names to a valid format
#' @param name_test string indicating sample identification
#' @return string with \code{name_test} contained within backquote
#' @examples make_name_valid("CPA 1 | Slide 1 | Seg-A")
make_name_valid <- function(name_test) {
  if (make.names(name_test) != name_test) {
    name_test <- gsub("`", "", name_test)
    name_test <- paste0("`", name_test, "`")
  }
  return(name_test)
}

#get genes of interest out of string
genes_of_interest <- function (targets,gene_group){
  if (gene_group[1] == "all"){
    genes <- as.character(targets$HUGOSymbol)
    w2kp <- grep(x=genes, pattern = "NegProbe")
    if (length(w2kp) > 0){
      genes <- genes[-w2kp]
    }
  }else{
    w2kp <- NULL
    for (i in 1:length(gene_group)){
      w2kp <- c(w2kp, grep(x = targets$TargetGroup, pattern = gene_group[i]))
    }
    genes <- unique(as.character(targets$HUGOSymbol[w2kp]))
  }

  w2kp <- grep(x=genes, pattern = ";")
  if (length(w2kp) > 0){
    geneSplit <- unlist(strsplit(as.character(genes[w2kp]), split = ";"))
    genes <- c(genes[-w2kp],geneSplit)
  }

  return(genes)
}

#set color scheme for analysis
set_colors <- function(){
  if(color_grouping %in% rownames(brewer.pal.info)) {
    col_gr <- brewer.pal.info[color_grouping, "maxcolors"]
    getPaletteGroup = colorRampPalette(brewer.pal(col_gr, color_grouping))
  } else if(exists(paste0("pal_", color_grouping))) {
    col_gr <- suppressWarnings(get(paste0("pal_", color_grouping))()(51))
    col_gr <- col_gr[!is.na(col_gr)]
    getPaletteGroup = colorRampPalette(col_gr)
    col_gr <- length(col_gr)
  } else {
    stop('Could not identify color_grouping palette, please check against brewer.pal.info or ggsci palettes\n')
  }

  if(color_samples %in% rownames(brewer.pal.info)) {
    col_sm <- brewer.pal.info[color_samples, "maxcolors"]
    getPaletteSamples = colorRampPalette(brewer.pal(col_sm, color_samples))
  } else if(exists(paste0("pal_", color_samples))) {
    col_sm <- suppressWarnings(get(paste0("pal_", color_samples))()(51))
    col_sm <- col_sm[!is.na(col_sm)]
    getPaletteSamples = colorRampPalette(col_sm)
    col_sm <- length(col_sm)
  } else {
    stop('Could not identify color_samples palette, please check against brewer.pal.info or ggsci palettes\n')
  }

  if(color_genes %in% rownames(brewer.pal.info)) {
    col_tg <- brewer.pal.info[color_genes, "maxcolors"]
    getPaletteGenes = colorRampPalette(brewer.pal(col_tg, color_genes))
  } else if(exists(paste0("pal_", color_genes))) {
    col_tg <- suppressWarnings(get(paste0("pal_", color_genes))()(51))
    col_tg <- col_tg[!is.na(col_tg)]
    getPaletteGenes = colorRampPalette(col_tg)
    col_tg <- length(col_tg)
  } else {
    stop('Could not identify color_genes palette, please check against brewer.pal.info or ggsci palettes\n')
  }


  groups <- unique(dfs[[3]][[grouping_var]])
  gr_ln <- length(groups)
  groups.colors <- getPaletteGroup(max(gr_ln, col_gr))[1:gr_ln]
  samples <- unique(dfs[[3]][[control_var]])
  sm_ln <- length(samples)
  samples.colors <- getPaletteSamples(max(sm_ln, col_sm))[1:sm_ln]
  samples.dr <- unique(dfs[[3]][[color_by]])
  smdr_ln <- length(samples.dr)
  samples.dr.colors <- getPaletteSamples(max(smdr_ln, col_sm))[1:smdr_ln]
  samples.concat <- unique(paste(dfs[[3]][[grouping_var]], dfs[[3]][[control_var]]))
  smcc_ln <- length(samples.concat)
  samples.concat.colors <- getPaletteSamples(max(smcc_ln, col_sm))[1:smcc_ln]
  targets <- c(gene_group, "Fav.Genes", "Not Specified", "Multiple")
  tg_ln <- length(targets)
  targets.colors <- getPaletteGenes(max(tg_ln, col_tg))[1:tg_ln]

  colors <- list(
    c(groups = groups.colors),
    c(samples = samples.colors),
    c(samples.dr = samples.dr.colors),
    c(samples.concat = samples.concat.colors),
    c(targets = targets.colors)
  )

  names(colors) <- c(grouping_var, control_var, color_by, "AOIs", "gene_groups")
  w2kp <- which(!duplicated(names(colors)))
  if (length(w2kp) > 0){
    colors <- colors[w2kp]
  }

  names(colors[[grouping_var]]) <- as.character(groups)
  w2kp <- which(is.na(names(colors[[grouping_var]])))
  if (length(w2kp) > 0){
    colors[[grouping_var]] <- colors[[grouping_var]][-w2kp]
  }

  names(colors[[control_var]]) <- as.character(samples)
  w2kp <- which(is.na(names(colors[[control_var]])))
  if (length(w2kp) > 0){
    colors[[control_var]] <- colors[[control_var]][-w2kp]
  }

  names(colors[[color_by]]) <- as.character(samples.dr)
  w2kp <- which(is.na(names(colors[[color_by]])))
  if (length(w2kp) > 0){
    colors[[color_by]] <- colors[[color_by]][-w2kp]
  }

  names(colors[["AOIs"]]) <- as.character(samples.concat)
  w2kp <- which(is.na(names(colors[["AOIs"]])))
  if (length(w2kp) > 0){
    colors[["AOIs"]] <- colors[["AOIs"]][-w2kp]
  }

  names(colors[["gene_groups"]]) <- as.character(targets)

  if ("All Probes" %in% gene_group){
    colors[['gene_groups']][["All Probes"]] <- '#ADADAD'
  }

  colors[['gene_groups']][['Not Specified']] <- '#ADADAD'

  flog.info("\n \n ############   Color Schemes   ################## \n", name="DSP_NGS_log")
  for (i in 1:length(colors)){
    flog.info(paste0(names(colors[i]), ":"), name="DSP_NGS_log")
    line = ""
    for (j in 1:length(colors[[i]])){
      line <- paste(line, paste(names(colors[[i]][j]), colors[[i]][j], sep = " = "), ", ", sep = "")
    }
    flog.info(str_wrap(substr(line, 1, nchar(line)-2), width = 120), name="DSP_NGS_log")
  }



  return(colors)
}

# parse targetNotes to create geneSet lists
parse_GeneSets <- function(notes = dfs[[5]], gene_group) {

  out <- list()
  for (i in 1:nrow(notes)){
    groups <- unlist(strsplit(as.character(notes$TargetGroup[i]), ";", fixed = TRUE))
    for (group in groups){
      if (!is.null(out[[group]])){
        out[group] <- list(c(out[[group]], as.character(notes$HUGOSymbol[i])))
      }else{
        out[group] <- list(as.character(notes$HUGOSymbol[i]))
      }
    }
  }
  if (!"All Probes" %in% gene_group){
    out[["All Probes"]] <- NULL
  }

  return(out)
}

write_dim_reduct_files <- function(pca = pca, umap = umap, tsne = tsne){
  sample_dims <- data.frame(pca[["samples"]]$call$X[,1:2])

  #ensure same order
  order <- match(rownames(sample_dims), tsne[["samples"]]$Sample_ID)
  tsne[["samples"]] <- tsne[["samples"]][order,]
  order <- match(rownames(sample_dims), umap[["samples"]]$Sample_ID)
  umap[["samples"]] <- umap[["samples"]][order,]

  sample_dims <- cbind(sample_dims,tsne[["samples"]]$tsne, umap[["samples"]]$umap, tsne[["samples"]]$cluster, tsne[["samples"]]$gene_count)
  names(sample_dims) <- c("PCA1", "PCA2","tSNE1", "tSNE2", "UMAP1", "UMAP2", "cluster", "gene_count")
  samples_dims <- cbind(sample_dims, tsne[["samples"]][,2:(which(colnames(tsne[["samples"]])== "RawReads")-1)])

  write.csv(x = samples_dims, file = paste0(samples_dir, "/coordinates.csv"))

  gene_dims <- data.frame(pca[["genes"]]$call$X[,1:2])

  #ensure same order
  order <- match(rownames(gene_dims), tsne[["genes"]]$gene)
  tsne[["genes"]] <- tsne[["genes"]][order,]
  order <- match(rownames(gene_dims), umap[["genes"]]$gene)
  umap[["genes"]] <- umap[["genes"]][order,]


  gene_dims <- cbind(gene_dims,tsne[["genes"]]$X1, tsne[["genes"]]$X2, umap[["genes"]]$X1,
                     umap[["genes"]]$X2, tsne[["genes"]]$cluster, tsne[["genes"]]$gene_count,
                     tsne[["genes"]]$color, tsne[["genes"]]$group_color)
  names(gene_dims) <- c("PCA1", "PCA2","tSNE1", "tSNE2", "UMAP1", "UMAP2", "cluster", "gene_count", "Gene_group", grouping_var)

  write.csv(x = gene_dims, file = paste0(genes_dir, "/coordinates.csv"))
}

#' @title function to calculate the geometric mean of a vector
#' @param v numeric vector
#' @return geometric mean of \code{v}
#' @examples ngeoMean(c(0, 1, 5, 8, 3))
ngeoMean <- function(v) {
  v[v == 0] <- 1
  return(geoMean(v, na.rm = T))
}

#' @title function to calculate normalization factors based on geometric mean
#' @param v numeric vector
#' @return vector of nnumeric normalization factors for \code{v}
#' @examples nanoNormFactors(c(0, 1, 5, 8, 3))
nanoNormFactors <- function(v) {
  return(v / ngeoMean(v))
}
