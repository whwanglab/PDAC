#' DSP_hydra VERSION 0.3
#'
#' Pipeline for DSP NGS analysis
#'
#' @param data_path path to dsp dataset excel workbook
#' @param config Configuration file, see config.yml example file
#'
#' @return Resulting visualizations
#'
#' @examples
#' read_dataset(data_path)


#' @section Load Libraries

library(dspNgs)

#set working directory to location of DSP_hydra.R
path <- getActiveDocumentContext()$path
setwd(dirname(path))

#install Cell Decon Package
if(!"InSituSort" %in% installed.packages()[,"Package"]){
  install.packages("logNormReg")
  if(sessionInfo()[2]$platform == "x86_64-pc-linux-gnu (64-bit)"){
    install.packages("../pkgs/InSituSort_0.0.0.9001/InSituSort_0.0.0.9000.tar.gz", repos = NULL, type = "source")
  }else{
    install.packages("../pkgs/InSituSort_0.0.0.9001/InSituSort_0.0.0.9000.zip", repos = NULL, type = "win.binary")
  }
}else{
  platform <- sessionInfo()[2]$platform
  remove.packages("InSituSort")
  if(platform == "x86_64-pc-linux-gnu (64-bit)"){
    install.packages("../pkgs/InSituSort_0.0.0.9001/InSituSort_0.0.0.9000.tar.gz", repos = NULL, type = "source")
  }else{
    install.packages("../pkgs/InSituSort_0.0.0.9001/InSituSort_0.0.0.9000.zip", repos = NULL, type = "win.binary")
  }
}

library(InSituSort)

#' @section Configurtaion
#' Read in options from configuration file
read_config()

if(!is.null(RData)){
  if(file.exists(RData)){
    load(RData)
    flog.logger(name="DSP_NGS_log", TRACE, appender=appender.file(log_file))
    # format output to include -l log level, -n namespace, -f function, -m message
    layout <- layout.format('[~l] [~n.~f] ~m')
    # apply layout
    flog.layout(layout)
  }else{
    print("WARNING: RData file path does not exist. Check file path in config.yml if you want to reload a workspace and rerun read_config()")
    stop()
  }
}

#' @section Setup Logging

# capture today's date
today = format(Sys.time(), "%Y%m%d_%H%M")

outdir <- paste(outdir, outname, sep = "/")
dir.create(outdir)

log_dir <- paste(outdir, "log/", sep = "/")
dir.create(log_dir)

#copy config file into log folder
file.copy(from = "config.yml", to = paste0(log_dir, "/config.yml"), overwrite = TRUE, copy.mode = TRUE)

# create log file path
log_file = paste0(log_dir, today, "_DSP_NGS.log")

# setup a logger to handle TRACE level and greater stdout (ERRORS, WARNINGS, etc.)
# tee outputs to both file and stdout
flog.logger(name="DSP_NGS_log", TRACE, appender=appender.file(log_file))
# format output to include -l log level, -n namespace, -f function, -m message
layout <- layout.format('[~l] [~n.~f] ~m')
# apply layout
flog.layout(layout)

flog.info("\n \n ############   Loaded Packages  ################## \n", name="DSP_NGS_log")

# output attached packages and versions
packinfo <- installed.packages(fields = c("Package", "Version"))
flog.info(capture.output(sessionInfo()), name="DSP_NGS_log")

#' @section Multiprocessing

# use 65% of the cores
numCores  = round(detectCores() * .65)
flog.info("Using %i cores.", numCores, name="DSP_NGS_log")

#' @section Create Directories
if(!endsWith(outdir, "/")){
  outdir <- paste(outdir, "/", sep = "")
}

create_dirs(outdir)

#' @section Read in Dataset
print("Reading in dataset")

# DSP-DA data wrangling
if (txt == FALSE) {
  dfs <- read_dataset(xl_path=data_path, raw_xl=DA_raw,
                        norm_xl=DA_norm, norm_type=norm_method)
  #dfs[[1]] raw collapsed counts
  #dfs[[2]] normalized counts
  #dfs[[3]] ROI/AOI annotations
  #dfs[[4]] summary of experiment
  #dfs[[5]] target notes

  # Get name of normalized data analysis
  pjid <- gsub(".xlsx", "", DA_norm)

# Azorius data wrangling
} else {
  counts <- read.delim(dir(data_path, full.names=T, pattern="_TargetCountMatrix.txt$")[2], sep = "\t",
                       header = TRUE)
  counts[is.na(counts)] <- 1
  norm_counts <- read.delim(dir(data_path, full.names=T, pattern="_NegNorm"),
                            sep = "\t", header = TRUE)
  samp_notes <- read.delim(dir(data_path, full.names=T, pattern="_SegmentProperties.txt$"), sep = "\t",
                           header = TRUE)
  samp_notes$Sample_ID <- str_replace_all(samp_notes$Sample_ID, "-", ".")
  if(any(grepl(x = samp_notes$Sample_ID, pattern = "^\\d"))){
    samp_notes$Sample_ID <- sub("^", "X", samp_notes$Sample_ID)
  }
  proc_sum <- read.delim(dir(data_path, full.names=T, pattern="_DatasetHistory.txt$"), sep = "\t",
                         header = TRUE)
  probe_notes <- read.delim(dir(data_path, full.names=T, pattern="_TargetProperties.txt$"), sep = "\t",
                            header = TRUE, colClasses = c(Pooling='character'))
  if("GlobalOutliers" %in% colnames(probe_notes)) {
    probe_notes <- droplevels(probe_notes[!is.na(probe_notes[["GlobalOutliers"]]),])
  } else {
    warning(paste("GlobalOutliers column missing from target properties file.",
              "No outlier information has been passed from Azorius."))
  }
  dfs <- list(counts, norm_counts, samp_notes, proc_sum, probe_notes)
  rm(counts, norm_counts, samp_notes, proc_sum, probe_notes)

  pjid <- gsub(".*/(.*)_DatasetHistory.txt$", "\\1", dir(data_path, full.names=T, pattern="_DatasetHistory.txt$"))
}

# For Azorius generated data only
if (txt) {
  if(length(grep(x = colnames(dfs[[3]]), pattern = "GeoLOD")) > 0){
    print("ERROR: data generated using old Azorius. Please rerun new Azorius before continuing")
    stop()
  }

  #' @section Change LOQ if needed
  dfs[[3]] <- change_LOQ(dfs, LOQ_level)
}

#sample removal if needed
if(!is.null(remove_samples_col)){
  dfs <- drop_samples(df=dfs,
                      cols = remove_samples_col,
                      conditions = remove_samples_condition,
                      keep = remove_samples_keep,
                      operators = remove_samples_operators,
                      and = remove_samples_and,
                      custom = remove_samples_custom)
}

#' @section Color Selection
colors <- set_colors()

#' @section Set Theme
if(is.null(preset_theme)) {
  # Override any pre-existing themes with ggplot default
  theme_set(theme_gray())
} else {
  theme_set(eval(parse(text=preset_theme)))
}
# Add individual updates to theme
if(!is.null(theme_adjustments)) {
  do.call(theme_update, eval(parse(text=theme_adjustments)))
  preset_theme <- theme_get()
}

#' @section Subset Gene Group
#' Create list of genes in gene_group pathway of interest
targetNotes <- parse_GeneSets(notes = dfs[[5]], gene_group = gene_group)
genes <- targetNotes[gene_group]


if (0 %in% lapply(genes, length)){
  print("Error: The given gene_group(s) in the config file do not match the probe groups in the TargetProperties file. Please check that these are spelled and/or capitalized correctly and restart the analysis")
}

wrap_num <- 45
names(colors$gene_groups) <- str_wrap(names(colors$gene_groups), wrap_num)

#' @section Squencing QC
run_seqQC(dataframelist=dfs,
          loq=LOQ_level,
          qc_dir=qc_dir,
          facet_column='',
          fileType=fileType)

#' @section Normalization
# Normalize data using chosen factor

# normalization plots to determine best one
print("Different normalization plots to determine best one to use")

#run QC plots
qc_plots(data_matrix = dfs,
           grp_var = grouping_var,
           ctrl_var = control_var,
           HKs = hk_genes,
           outdir = qc_dir,
           colors = colors,
           txt_type = txt)

#***************************************************************************************************
#***IF YOU WANT TO CHANGE THE NORMALIZATION FACTOR OR ANYTHING ELSE IN THE CONFIG FILE DO IT NOW****
#***************************************************************************************************

# For Azorius generated data only
if (txt) {
  read_config()

  if(!is.null(norm_method)) {
    # renormalize if needed
    print("Re-normalizing data")

    dfs[[2]] <- qc_plots(data_matrix = dfs,
                         norm_method = norm_method,
                         grp_var = grouping_var,
                         ctrl_var = control_var,
                         HKs = hk_genes,
                         outdir = qc_dir,
                         colors = colors)

    colnames(dfs[[2]])[1] <- "TargetName"
  }
}

#***************************************************************************************************
#***Dropping blacklist genes, outliers, and genes below LOQ cutoffs*********************************
#***************************************************************************************************

# For Azorius generated data only
if (txt) {
  # Drop black list genes
  dfs[[2]] <- remove_blacklist_genes(dfs, rm_genes)

  # Remove outliers and reassociate the counts
  dfs[[2]] <- remove_gene_outliers(df = dfs, outlier_cutoff = outlier_cutoff,
                                     thresh = thresh)

  dfs[[2]] <- remove_LOQ(dfs = dfs, LOQ_level = LOQ_level, LOQ_cutoff = LOQ_cutoff)

  # If Q3 was the normalization method chosen, renormalize Q3 with counts after dropping LOQ genes, if less than 1000, stick with original Q3 factors
  if(norm_method == "Q3"){
    dfs <- run_Q3(dataset = dfs)
  }
} else {
  colnames(dfs[[2]])[1] <- "Gene"
}

#outputs to output folder
write.table(x = dfs[[2]], file = paste0(results_dir, "/", pjid, "_", norm_method, "Norm_TargetCountMatrix.txt"), row.names = F, sep = "\t", quote = FALSE)



#' @section Dimensional Reduction
#' Run PCA, tSNE, UMAP on ROIs and Genes

print("Performing dimensional reduction on samples")
pca <- list()
tsne <- list()
umap <- list()

pca[["samples"]] <- run_PCA(df = dfs,
                            outdir = samples_dir,
                            color = color_by,
                            symbol = shape_by,
                            size = size,
                            colors = colors)

tsne[["samples"]] <- run_tsne(df = dfs,
                              pca = pca[["samples"]],
                              outdir = samples_dir,
                              color = color_by,
                              size = size,
                              symbol = shape_by,
                              perplexity = perplexity,
                              colors = colors)

umap[["samples"]] <- run_umap(df = dfs,
                              pca = pca[["samples"]],
                              outdir = samples_dir,
                              color = color_by,
                              size = size,
                              symbol = shape_by,
                              colors = colors)


#' @section ROI plots: expression
if(length(ROI_ID) > 1) {
  dfs[[3]]$ROI_ID <- make_ID(IDs = ROI_ID, df = dfs[[3]])
  ROI_ID <- 'ROI_ID'
}

ROI_plots <- plot_ROI2D(data_df = dfs[[2]],              # input data expression set (normalized)
                        annot_df = dfs[[3]],             # input annotations
                        ROI_ID = ROI_ID,                 # ROI column
                        AOI_ID = AOI_ID,                 # AOI column
                        plot_type = ROI_plot_type,       # ROI plot type: 'joy','heatmap','bar'
                        title = 'Gene Expression',       # Additional title to be plotted on ROI graphs for future reference
                        skip_ROIs = FALSE,               # whether ROIs should be graphed
                        targets = targets_ROI_plot,      # list of targets or method for identifying top genes
                        target_thresh = target_thresh,   # threshold value. %ile of data to plot. 0 (all) - 1 (none)
                        plot_clusters = TRUE,            # whether to plot clusters from dendrogram on ROI graph
                        cluster_pal = color_grouping,    # color for palette plotting clusters onto ROI graph
                        cluster_n = cluster_n,           # number of clusters to make
                        cluster_min = 5,                 # minimum size at which to draw breaks on top of graph
                        draw_bgd = TRUE,                 # whether the background should be drawn TRUE or FALSE
                        split_groups = TRUE,             # whether lines should be drawn indicating clusters
                        bgd_color = 'white',             # line around plot (50% transparent black) behind bars or joy plot
                        bgd_fill = alpha('black',0.5),   # background color for plot (50% transparent black) behind bars or joy plot
                        segment_ann = segment_Annotation,# segment annotation column
                        segment_names = AOI_names,       # segment names to be plotted based on annotation column
                        segment_colors = AOI_colors,     # segment colors to be used in graphing
                        transf = target_transformation,  # transformation to be applied to all expression values
                        hole_diameter = 5,               # relative value for how large the window into the graph should be (higher = larger)
                        width = 1)                       # relative value for how wide the data plotting area should be (higher = larger
if(save_ROI_plots) {
  counter <- 0
  pb <- txtProgressBar(min = 0, max = length(ROI_plots)-3, style = 3)
  pdf(file = paste0(spatial_dir,
                    'ROI_ExpressionPlots_',
                    paste(AOI_names, collapse = "_"),
                    '.pdf'))

  cat('\n  Writing plots to file, this will take a few minutes\n')

  print(ROI_plots$legend + labs(title = 'Dendrogram Clustering Key for Graphs'))
  print(ROI_plots[[4]]$interactive + labs(title = 'Linear Representation of ROI Graphs'))
  for(roi in 4:length(ROI_plots)) {
    print(ROI_plots[[roi]]$plot)
    counter <- counter + 1
    setTxtProgressBar(pb = pb, value = counter)
  }
  dev.off()
}

# plot an interactive graph (change 4 to the ROI ID):
ggplotly(ROI_plots[[4]]$interactive)

if(xy_axes == "pca"){
  x_position <- pca[["samples"]]$call$X$Dim.1
  y_position <- pca[["samples"]]$call$X$Dim.2
}else if(xy_axes == "tsne"){
  x_position <- tsne[["samples"]]$tsne$X1
  y_position <- tsne[["samples"]]$tsne$X2
}else if(xy_axes == "umap"){
  x_position <- umap[["samples"]]$umap$X1
  y_position <- umap[["samples"]]$umap$X2
}else{
  x_position <- dfs[[3]][[x_position_col]]
  y_position <- dfs[[3]][[y_position_col]]

  if(is.null(x_position) | is.null(y_position)){
    print("Error: The x or y position columns in the config file to not match the column names in the Segment Properties file. Please check that these are in the file and spelled and/or capitalized correctly and restart the analysis")
  }
}


#' @section Cell type deconvolution
print("Cell Type Deconvolution")
decon <- run_cell_decon(df = dfs,
                        matrix = profile_matrix,
                        norm_method = norm_method,
                        tumor_high_ids = high_tumor_AOIs,
                        x = x_position,
                        y = y_position,
                        xlab = x_label,
                        ylab = y_label,
                        segment = segmentation,
                        skip = skipped_seg,
                        outdir = cellType_dir)


#' @section Differential Expression
print("Differential Expression")
# Runs the DE analysis in parallel
de_results <- dsp_de_analysis_parallel(de_results = de_results,
                                  norm_counts = dfs[[2]],
                                  samp_notes = dfs[[3]],
                                  grouping_var = grouping_var,
                                  base_level = base_level,
                                  control_var = control_var,
                                  de_dir = de_dir,
                                  n_processors = numCores)

# Feeds the de_results from above to process additional plots/analyses
de_results <- dsp_de_analysis(de_results = de_results,
                              norm_counts = dfs[[2]],
                              genes = genes,
                              grouping_var = grouping_var,
                              base_level = base_level,
                              results_dir = results_dir,
                              de_dir = de_dir,
                              control_var = control_var,
                              residuals_dir = residuals_dir,
                              qq_dir = qq_dir,
                              pval_cutoff = pval_cutoff,
                              fc_cutoff = fc_cutoff,
                              samp_notes = dfs[[3]])

#' @section Top GeneSet Summary
flog.info(paste0("Plotting top gene sets from annotations"), name="DSP_NGS_log")
for(comp in unique(de_results$test)) {
  GS_hits <- find_topGS(de_results = subset(de_results, test == comp),
                        targetNotes = genes,
                        fileType = fileType,
                        outdir = de_dir,
                        name = comp)
  GS_hits$test <- comp
  if(comp == unique(de_results$test)[1]) {
    GS_hits_table <- GS_hits
  } else {
    GS_hits_table <- rbind(GS_hits_table, GS_hits)
  }
}

GS_hits_all <- find_topGS(targetNotes = targetNotes,
                          de_results = de_results,
                          fileType = fileType,
                          outdir = de_dir,
                          name = 'All DE Tests')

#' @section Heatmaps for Selected Genes
flog.info(paste0("Creating Heatmaps for Gene Sets and Favorite Genes"), name="DSP_NGS_log")

# geneset genes
for(geneset in gene_group) {
  if (geneset != "All Probes"){
    run_heatmap(data = dfs[[2]],
                annot = dfs[[3]],
                samp_vars = c(control_var, grouping_var),
                fav_genes = fav_genes,
                ann_colors = colors,
                targets = genes[[geneset]],
                name = geneset,
                outdir = clustering_dir,
                fileType = fileType)
  }
}

if(!is.null(fav_genes)){
  # favorite genes
  run_heatmap(data = dfs[[2]],
              annot = dfs[[3]],
              samp_vars = c(control_var, grouping_var),
              fav_genes = NULL,
              ann_colors = colors,
              targets = fav_genes,
              name = 'Fav.Genes',
              outdir = clustering_dir,
              fileType = fileType)
}


# top DE genes (all contrasts)
tmp_de <- de_results[order(de_results$Pval, decreasing = FALSE), ]
sig_genes <- unique(unlist(tmp_de$gene)[1:n_top_heatmap])

run_heatmap(data = dfs[[2]],
            annot = dfs[[3]],
            samp_vars = c(control_var, grouping_var),
            fav_genes = NULL,
            ann_colors = colors,
            targets = sig_genes,
            name = 'Top DE Genes',
            outdir = clustering_dir,
            fileType = fileType,
            sortBy = grouping_var,
            sortBaseline = base_level)

#dimensional reduction on genes
print("Performing dimensional reduction on genes")

de_genes <- de_genes_by_group()

pca[["genes"]] <- run_PCA_genes(df = dfs,
                                outdir = genes_dir,
                                size = size,
                                expressed_genes = de_genes)

tsne[["genes"]] <- run_tsne_genes(df = dfs,
                                  pca = pca[["genes"]],
                                  outdir = genes_dir,
                                  size = size,
                                  perplexity = perplexity * 10,
                                  expressed_genes = de_genes)

umap[["genes"]] <- run_umap_genes(df = dfs,
                                  pca = pca[["genes"]],
                                  outdir = genes_dir,
                                  size = size,
                                  expressed_genes = de_genes)

#' @section Pathway analysis
#' run pathway analysis on given ROIs

flog.info(paste0("Performing Pathway Analysis"), name="DSP_NGS_log")
enrich_res <- run_pathways(df = dfs,
                         de_results = de_results,
                         geneListrank = geneListrank,
                         path_fc = path_fc,
                         path_pval = path_pval,
                         enrichment = enrichment,
                         outdir = path_dir,
                         custom_pathways_path = custom_pathways_path,
                         exclude_reactome = exclude_reactome)

flog.info(paste0("Creating Volcano Plot for Pathway Analysis"), name="DSP_NGS_log")
for(comp in names(enrich_res)) {

  # pull each enrichment data
  temp_res <- enrich_res[[comp]]

  # add significance and create volcano plot
  temp_res$Significance <- -log10(temp_res$pval)
  temp_res$pathway_wrap <- str_wrap(temp_res$pathway, 25)
  plot_volcano(de = temp_res,
               negative_label = strsplit(comp, split=" vs ")[[1]][1],
               positive_label = strsplit(comp, split=" vs ")[[1]][2],
               target_group = 'Pathway',
               fav_targets = NULL,
               targets = unlist(temp_res$pathway_wrap),
               outfile = paste0(path_dir, "/", comp, "/", comp, '_GSEA_Volcano-AllPathways.', fileType),
               point_color = default_color,
               top_method = 'Significance',
               n_targets = 10,
               plt_title = paste0('GSEA Volcano Plot for ', comp),
               target_ID = 'pathway_wrap',
               FC_ID = 'NES',
               Pval_ID = 'pval',
               color_list = NULL,
               save_plot = TRUE)

  pathways <- distinct(temp_res, pathway, .keep_all = TRUE)
  pathways <- head(pathways[order(pathways$NES),], n = 10)

  overlapOutput <- pathway_overlap(pathways,
                                   all = TRUE,
                                   species_data = species)


  write.csv(apply(overlapOutput,2,as.character),
            paste0(path_dir, "/", comp, "/", comp, '_pathways_instersect.csv'),
            row.names = FALSE)

  pathway_intersect(pathways = pathways,
                    species_data = species,
                    outfile = paste0(path_dir, "/", comp, "/", comp, '_GSEA_Pathways_Intersect.', fileType),
                    save_plot = TRUE)

  pathway_chord(pathways = pathways,
                species_data = species,
                outfile = paste0(path_dir, "/", comp, "/", comp, '_GSEA_Pathways_Chord.', fileType),
                save_plot = TRUE)


}

write_dim_reduct_files(pca = pca,
                       umap = umap,
                       tsne = tsne)

save.image(file = paste0(RData_dir, "/", today, "_RData"))
