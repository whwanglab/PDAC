#' Reads in DSP dataset
#'
#' @export parallel_read
#' @export read_dataset
#' @export reformat_dataset
#' @export remove_gene_outliers
#' @export remove_LOQ
#' @export remove_blacklist_genes
#' @export change_LOQ
#' @export drop_expected_negs
#' @export drop_samples

#' @title Function that parallelizes spreadsheet read
#' @param file string value with file name
#' @param get_cols boolean indicating if column names should be preserved
#' @return list of data.frame objects containing data from /code{file} tabs
#' @examples
#' parallel_read(file="./testData/DA_data/Q3_normalization/Q3 norm.xlsx")
parallel_read <- function(file, get_cols=TRUE){

  # detect available cores
  numCores  = round(parallel::detectCores() * .65)

  # get sheet names
  sheets <- readxl::excel_sheets(file)

  if(.Platform$OS.type == "windows") {
    print("OS = Windows")
    cl <- makePSOCKcluster(numCores)
    tmp_dfs <-
      parLapplyLB(cl,
                    sheets,
                    function(sheet, file) {
                      readxl::read_excel(file, sheet=sheet, col_names=get_cols)
                    },
                    file)
    stopCluster(cl)
  } else {
    tmp_dfs <- parallel::mclapply(excel_sheets(file),
                                    readxl::read_excel,
                                    path = file,
                                    col_names = get_cols,
                                    mc.cores=numCores)
  }
  names(tmp_dfs) <- sheets
  return(tmp_dfs)
}

#' @title function to reformat DSP-DA data to be compatible with Hydra
#' @param curr_dfs list of current data frames resulting
#'   from \code{\link{read_dataset}}
#' @param norm_type string identifying the type of normalization used
#' @return list of data frames containing reformatted /code{curr_dfs}
#' @examples
#' Q3_dfs <- read_dataset(xl_path="./testData/DA_data/Q3_normalization",
#'                          raw_xl="ProbeQC - Default.xlsx",
#'                          norm_xl="Q3 Norm.xlsx",
#'                          norm_type="Q3")
#' reformat_dataset(curr_dfs=Q3_dfs, norm_type="Q3")
reformat_dataset <- function(curr_dfs, norm_type) {

  # Remove appended sheet names for column names read in
  for (idx in names(curr_dfs)) {
    if (nrow(curr_dfs[[idx]]) < 1) {
      stop(paste0("Missing excel tab with ", idx, " data."))
    }
    # Associate column names with appropriate format
    colnames(curr_dfs[[idx]]) <- curr_dfs[[idx]][1, ]
    curr_dfs[[idx]] <- curr_dfs[[idx]][-1, ]
  }

  # dfs[[1]] Raw collapsed counts
  # Change counts to numeric and at least 1
  curr_dfs[[1]][,2:ncol(curr_dfs[[1]])] <-
    sapply(curr_dfs[[1]][,2:ncol(curr_dfs[[1]])], as.numeric)
  curr_dfs[[1]][is.na(curr_dfs[[1]])] <- 1
  colnames(curr_dfs[[1]]) <-
    sapply(colnames(curr_dfs[[1]]), make_name_valid)
  rownames(curr_dfs[[1]]) <- curr_dfs[[1]][["TargetName"]]

  # dfs[[2]] Normalized counts
  # Change counts to numeric
  curr_dfs[[2]][,2:ncol(curr_dfs[[2]])] <-
    sapply(curr_dfs[[2]][,2:ncol(curr_dfs[[2]])], as.numeric)
  colnames(curr_dfs[[2]]) <-
    sapply(colnames(curr_dfs[[2]]), make_name_valid)
  rownames(curr_dfs[[2]]) <- curr_dfs[[2]][["TargetName"]]

  # dfs[[3]] Segment properties
  # Replace SegmentDisplayName column name with Sample_ID and reorder
  curr_dfs[[3]][["Sample_ID"]] <- curr_dfs[[3]][["SegmentDisplayName"]]
  curr_dfs[[3]] <-
    select(curr_dfs[[3]], Sample_ID, everything(), -SegmentDisplayName)
  # Ensure naming format compatible
  curr_dfs[[3]][["Sample_ID"]] <-
    sapply(curr_dfs[[3]][["Sample_ID"]], make_name_valid)
  # Identify numeric columns
  cols_num <- c("RawReads", "AlignedReads", "DeduplicatedReads", "TrimmedReads", "StitchedReads")
  curr_dfs[[3]][cols_num] <- sapply(curr_dfs[[3]][cols_num], as.numeric)

  # dfs[[4]] Experiment summary
  # No changes to dataset summary data frame
  curr_dfs[[4]] <-
    curr_dfs[[4]]

  # dfs[[5]] Target properties
  # Make sure target notes has a Pooling column of character type
  colnames(curr_dfs[[5]])[colnames(curr_dfs[[5]]) == "ProbePool"] <- "Pooling"
  curr_dfs[[5]] <-
    transform(curr_dfs[[5]], Pooling=as.character(Pooling))
  # Replace DSP-DA naming convention
  curr_dfs[[5]][["TargetGroup"]] <- gsub("All Targets", "All Probes", curr_dfs[[5]][["TargetGroup"]])
  curr_dfs[[5]][["TargetGroup"]] <- gsub(",", ";", curr_dfs[[5]][["TargetGroup"]])

  # Relabel segment property LOQ, negative geomean, and normalization factor columns
  pool_ids <- unique(curr_dfs[[5]][["Pooling"]])
  LOQ_key <- "Standard deviation amount for the LOQ"
  for (pool in pool_ids) {
    # Get LOD used in DSP-DA to override user input
    DA_LOQ <-
      as.numeric(curr_dfs[[4]][curr_dfs[[4]][, 1] == LOQ_key &
                                 !is.na(curr_dfs[[4]][, 1] == LOQ_key), 2])
    assign("LOQ_level", DA_LOQ, envir=.GlobalEnv)

    geo_col <- paste0("GeoLOQ", LOQ_level, "_", pool)
    curr_dfs[[3]][[geo_col]] <- as.numeric(curr_dfs[[3]][["LOQ"]])

    # Pull out single negative geometric mean from raw collapsed counts for DSP-DA v2.0
    neg_name <-
      curr_dfs[[5]][curr_dfs[[5]][["CodeClass"]] == "Negative" &
                      curr_dfs[[5]][["Pooling"]] == pool, "TargetName"]
    neg_geo <- curr_dfs[[1]][curr_dfs[[1]][["TargetName"]] == neg_name,
                               2:ncol(curr_dfs[[1]])]
    if(nrow(neg_geo) != 1) {
      stop("Negatives not mapping 1:1.
             Make sure negative target names unique for each pool and
             each negative exists in collapsed count data.")
    }
    neg_geo_col <- paste0("NegGeoMean_", pool)
    curr_dfs[[3]][[neg_geo_col]] <- unlist(neg_geo[, curr_dfs[[3]][["Sample_ID"]]])

    # Accomodates DSP-DA v2.0 with one norm factor for all pools
    # This feature only accurate for single panel results from DA v2.0
    if (norm_type == "Neg") {
      neg_col <- paste0("NormFactorNeg_", pool)
      curr_dfs[[3]][[neg_col]] <- as.numeric(curr_dfs[[3]][["NormalizationFactor"]])
    }
  }
  # Add column with correct column header for Q3 or HK normalization factors
  if (norm_type != "Neg") {
    norm_col <- paste0("NormFactor", norm_type)
    curr_dfs[[3]][[norm_col]] <- as.numeric(curr_dfs[[3]][["NormalizationFactor"]])
  }

  return(curr_dfs)
}

#' @title Reads in DSP-DA excel workbook for analysis and visualiation
#' @param xl_path string specifying path to dsp dataset excel workbooks
#' @param raw_xl string specifying name of raw collapsed count workbook
#' @param norm_xl string specifying name of normalized count workbook
#' @param norm_type string indicating type of normalization used in DSP-DA.
#'  Accepted values: "Q3", "Neg", "HK"
#' @return list of dataframes, one for each tab
#' @examples
#' read_dataset(xl_path="./testData/DA_data/Q3_normalization",
#'                raw_xl="ProbeQC - Default.xlsx",
#'                norm_xl="Q3 Norm.xlsx",
#'                norm_type="Q3")
read_dataset <- function(xl_path, raw_xl, norm_xl, norm_type){

  # Contains two tabs from initial dataset export, one containing raw counts
  raw_dfs <- parallel_read(paste(xl_path, raw_xl, sep="/"), get_cols=FALSE)
  # Contains five tabs from scaled/normalized data export
  norm_dfs <- parallel_read(paste(xl_path, norm_xl, sep="/"), get_cols=FALSE)

  # Associate dataframes with corresponding excel sheet data
  DA_dfs <- list()
  DA_dfs[[1]] <- data.frame(raw_dfs["TargetCountMatrix"]) # counts
  DA_dfs[[2]] <- data.frame(norm_dfs["TargetCountMatrix"]) # normalized counts
  DA_dfs[[3]] <- data.frame(norm_dfs["SegmentProperties"]) # ROI/AOI annotations
  DA_dfs[[4]] <- data.frame(norm_dfs["Dataset summary"]) # summary of experiment
  # Use probe QC target properties since negative dropped by scaling
  DA_dfs[[5]] <- data.frame(raw_dfs["TargetProperties"]) # target notes
  names(DA_dfs) <- c("RawCountMatrix", "TargetCountMatrix",
                        "SegmentProperties", "Dataset summary", "TargetProperties")

  # Reformat to match txt dataframes format
  DA_dfs <- unname(reformat_dataset(curr_dfs=DA_dfs, norm_type=norm_type))

  # return list of dataframes
  return(DA_dfs)

}


# remove gene outliers from dataset
remove_gene_outliers <- function(df, outlier_cutoff, thresh){
  #read in norm_counts
  norm_counts <- df[[2]]
  rownames(norm_counts) <- norm_counts[, "TargetName"]
  norm_counts[["TargetName"]] <- NULL

  #genes removed with total gene count larger than
  #given outlier_cutoff x average gene count
  w2kp <- which(rowSums(norm_counts) > (mean(rowSums(norm_counts) * outlier_cutoff)))
  if (length(w2kp) != 0){
    flog.info("Genes removed with high gene count across all ROIs", name="DSP_NGS_log")
    flog.info(str_wrap(paste(row.names(norm_counts[w2kp,]), collapse = ", "), width = 120),
              name="DSP_NGS_log")
    print(paste0("Removed ", length(w2kp), " out of ", nrow(norm_counts), " genes with high gene count"))
    norm_counts <- norm_counts[-w2kp,]
  }

  samp_notes <- df[[3]]

  geomean <- grep(names(samp_notes), pattern = "GeoMean")

  geomean <- apply(samp_notes[, geomean, drop=FALSE], 1, median)

  w2kp <- which(cor(y = geomean, x = t(norm_counts), method = "spearman") > thresh)

  if (length(w2kp) != 0){
    flog.info("Genes removed with high correlation to negative probes across all ROIs",
              name="DSP_NGS_log")
    flog.info(str_wrap(paste(row.names(norm_counts[w2kp,]), collapse = ", "), width = 120),
              name="DSP_NGS_log")
    print(paste0("Removed ", length(w2kp), " out of ", nrow(norm_counts), " genes with high correlation to negative probes"))
    norm_counts <- norm_counts[-w2kp,]

  }

  norm_counts <- cbind(row.names(norm_counts), norm_counts)
  colnames(norm_counts)[1] <- "Gene"

  return(norm_counts)
}

remove_LOQ <- function(dfs, LOQ_level = LOQ_level, LOQ_cutoff = LOQ_cutoff,
                       grp_var = grouping_var, ctrl_var = control_var){
  counts <- dfs[[1]][,-1]
  rownames(counts) <- dfs[[1]][,1]
  w2kp <- which(rownames(counts) %in% dfs[[2]]$Gene)
  counts <- counts[w2kp,]

  targets <- dfs[[5]]
  rownames(targets) <- dfs[[5]]$HUGOSymbol
  # Filter out targets with no target information including NegProbe
  counts <- counts[rownames(counts) %in% rownames(targets), ]
  AOIs <- dfs[[3]]
  rownames(AOIs) <- AOIs$Sample_ID

  if(any(rownames(AOIs) != colnames(counts))){
    counts <- counts[,match(rownames(AOIs),colnames(counts))]
    AOIs <- AOIs[match(colnames(counts), rownames(AOIs)),]
  }

  #same code as from run_seqQC
  LOQs <- lapply(rownames(counts), function(x) {
    # Get counts for target
    row <- counts[x,]

    # Get Pool Name and Column
    pool <- targets[x,'Pooling']
    poolcol <- paste0('GeoLOQ', LOQ_level, '_', pool)

    loqs <- AOIs[,poolcol]
    loqs[loqs < LOQ_floor] <- LOQ_floor
    loqtest <- row > t(loqs)

    return(loqtest)
  })
  LOQs <- data.frame(matrix(unlist(LOQs), nrow=length(LOQs), byrow=T,
                            dimnames = list(rownames(counts), colnames(counts))))

  percents <- apply(LOQs, 2, function(x){sum(x)})

  Percent <- as.data.frame((percents/nrow(LOQs))*100)
  names(Percent) <- "Percent"

  AOI <- data.frame(AOIs$Sample_ID, AOIs[[grp_var]], AOIs[[ctrl_var]])
  names(AOI) <- c("AOI", "grp", "ctrl")

  order <- match(AOI$AOI, rownames(Percent))
  percents <- percents[order]

  data <- data.frame(AOI, Percent[order,])
  colnames(data)[4] <- "Percent"

  if(nrow(data) > 75){
    labels <- theme(axis.text.y=element_blank())
  }else{
    labels <- NULL
  }

  #barplot of percent of genes above LOQ in every AOI
  #colored by AOI Type (grp, ctrl variables)
  #AOI labels are printed if there are less than 76
  gp <- ggplot(data, aes(y=Percent, x=reorder(AOI, Percent), fill=paste(grp, ctrl), color=paste(grp, ctrl))) +
    geom_bar(stat = "identity", position = position_stack(reverse=TRUE)) +
    scale_y_discrete(limits = seq(0,100,10)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    scale_fill_manual(values=colors$AOIs) +
    scale_color_manual(values=colors$AOIs) +
    ggtitle("Percent of Genes at or above LOQ in every AOI") +
    xlab("AOI") +
    coord_flip() +
    theme(axis.ticks.y=element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "darkgray"),
          legend.title = element_blank())+
    labels


  ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ.", fileType), device = fileType, height = 8, width = 12, scale = 1.5)


  #Boxplots of percents of genes above LOQ in different AOI Types (grp, ctrl variables)
  gp <- ggplot(data, aes(x= paste(grp, ctrl), y = Percent, fill = paste(grp, ctrl)))+
    geom_boxplot()+
    scale_fill_manual(name = 'AOI type', values = colors$AOIs)+
    labs(x = "", y = "Percent of Genes at or above LOQ", fill = "AOI Type", title = "Percent of Genes at or above LOQ in different AOI types")+
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

  ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ_percent_boxplot.", fileType),
         device = fileType, height = 8, width = 8, scale = 1.5)

  data$count <- percents

  #Boxplots of gene counts above LOQ in different AOI Types (grp, ctrl variables)
  gp <- ggplot(data, aes(x= paste(grp, ctrl), y = count, fill = paste(grp, ctrl)))+
    geom_boxplot()+
    scale_fill_manual(name = 'AOI type', values = colors$AOIs)+
    labs(x = "", y = "Number of Genes at or above LOQ", fill = "AOI Type", title = "Number of Genes at or above LOQ in different AOI types")+
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

  ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ_counts_boxplot.", fileType),
         device = fileType, height = 8, width = 12, scale = 1.5)


  #scatter plot of gene counts above LOQ vs sat_numcrit variable
  #colored by AOI Type (grp, ctrl variables)
  if(!is.null(sat_numcrit)){
    data[[sat_numcrit]] <- AOIs[[sat_numcrit]]
    gp <- ggplot(data, aes(x = data[[sat_numcrit]], y = count, col = paste(grp, ctrl)))+
      geom_point()+
      labs(x = sat_numcrit, y = "Gene count at or above LOQ", col = "AOI Type", title = paste0("Gene count at or above LOQ vs ", sat_numcrit))+
      scale_color_manual(name = 'AOI type', values = colors$AOIs)

    ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ_counts_vs_", sat_numcrit, ".", fileType),
           device = fileType, height = 8, width = 12, scale = 1.5)

  }

  #scatter plot of gene counts above LOQ vs cell_counts
  #colored by AOI Type (grp, ctrl variables)
  if(!is.null(cell_count)){
    if(cell_count != sat_numcrit){
      data[[cell_count]] <- AOIs[[cell_count]]
      gp <- ggplot(data, aes(x = data[[cell_count]], y = count, col = paste(grp, ctrl)))+
        geom_point()+
        labs(x = cell_count, y = "Gene count at or above LOQ", col = "AOI Type", title = paste0("Gene count at or above LOQ vs ", cell_count))+
        scale_color_manual(name = 'AOI type', values = colors$AOIs)

      ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ_counts_vs_", cell_count, ".", fileType),
             device = fileType, height = 8, width = 12, scale = 1.5)
    }
  }

  above_LOQ <- apply(LOQs, 1, function(x){sum(x)})

  above_LOQ <- above_LOQ/ncol(LOQs)

  below_LOQ <- 1-above_LOQ

  #which genes are below LOQ in a higher percentage of AOIs than LOQ_cutoff
  #genes to be removed
  w2kp <- which(below_LOQ >= LOQ_cutoff)

  if (length(w2kp) != 0){
    flog.info(paste0("Genes removed below LOQ in ", LOQ_cutoff * 100,"% of AOIs", "(",
                     length(w2kp), "/", nrow(counts), ")"), name="DSP_NGS_log")
    if (length(w2kp) < 300){
      flog.info(str_wrap(paste(row.names(counts[w2kp,]), collapse = ", "), width = 120),
                name="DSP_NGS_log")
    }
    print(paste0("Removed ", length(w2kp), " out of ", nrow(counts), " genes"))
    norm_counts <- dfs[[2]][-w2kp,]
  }else{
    flog.info(paste0("Genes removed below LOQ in ", LOQ_cutoff * 100,"% of AOIs", "(0/",
                     nrow(counts), ")"), name="DSP_NGS_log")
    norm_counts <- dfs[[2]]
  }


  plt_df <- as.data.frame(above_LOQ)

  # add plot colors for multiple genesets
  plt_df$color <- 'Not Specified'
  if(tolower(names(genes[1])) != "all probes") {
    for (i in 1:length(genes)){
      w2kp <- which(rownames(plt_df) %in% unlist(genes[i]))
      w2kp2 <- which(plt_df$color[w2kp] != "Not Specified")

      plt_df$color[w2kp] <- str_wrap(names(genes[i]), wrap_num)
      plt_df$color[w2kp2] <- "Multiple"
    }
  }

  w2kp <- which(rownames(plt_df) %in% fav_genes)
  plt_df$color[w2kp] <- "Fav.Genes"
  plt_df$color[tolower(plt_df$color) == "all probes"] <- "Not Specified"
  plt_df$color <- factor(plt_df$color)

  #plt_df <- plt_df[which(plt_df$color != "Not Specified"),]
  plt_df$gene <- rownames(plt_df)

  hist <- hist(plt_df$above_LOQ, breaks=seq(0, 1, by = 0.02), plot =FALSE, right = TRUE)

  plt_df$y = 0
  for(i in 1:nrow(plt_df)){
    bin <- max(which(hist$breaks <= plt_df$above_LOQ[i]))
    bin <- ifelse(bin==1, bin, bin-1)
    plt_df$y[i] <- runif(1, min = 1, max = hist$counts[bin])
  }


  #histogram of percent of AOIs genes are above LOQ in
  #overlayed by scatterplot of individual genes in given gene sets or fav genes
  #dropped genes are shaded on histogram and labeled in scatterplot
  gp <- ggplot(plt_df, aes(x = above_LOQ, label = gene))+
    geom_histogram(fill = "Green3", breaks = hist$breaks)+
    labs(x = "Percent of AOIs individual genes are at or above LOQ in", y = "Gene Count",
         title = "Distribution of genes above LOQ in AOIs")+
    geom_rect(aes(xmin = -Inf, ymin = -Inf, ymax = Inf, xmax = 1-LOQ_cutoff), alpha = 0.006, fill = "grey55")+
    annotate("text", x = (1-LOQ_cutoff)/2, y = max(hist$counts)/2,
             label = paste(length(which(above_LOQ < (1-LOQ_cutoff))), "genes\nremoved"))+
    geom_point(data = plt_df[which(plt_df$color != "Not Specified"),], aes(y = y, color = color, x = above_LOQ - 0.01))+
    scale_color_manual(name = 'GeneSet', values = colors$gene_groups)+
    geom_text_repel(data = plt_df[plt_df$above_LOQ < 1-LOQ_cutoff & plt_df$color != "Not Specified",],
                    aes(y = y), color = 'black', box.padding = .4, point.padding = .2,
                    fontface = 'bold', force = 4, nudge_y = 5, min.segment.length = 0.1)

  ggsave(plot = gp, filename = paste0(qc_dir, "/LOQ_genes_removed_histogram.", fileType),
         device = fileType, height = 8, width = 12, scale = 1.5)

  data <- data[order(-data$Percent),]
  high <- subset(data, Percent >= signif(max(Percent)-10, digits = 1))
  low <- subset(data, Percent <= signif(min(Percent)+10, digits = 1))

  flog.info(paste("ROIs with more than", signif(max(Percent)-10, digits = 1),
                  "% of probes at or above LOQ, highest to lowest: ", high$AOI, high$grp, high$ctrl,
                  sep = " "), name="DSP_NGS_log")
  flog.info("\n\n", name="DSP_NGS_log")
  flog.info(paste("ROIs with less than", signif(min(Percent)+10, digits = 1),
                  "% of probes at or above LOQ, highest to lowest: ", low$AOI, low$grp, low$ctrl,
                  sep = " "), name="DSP_NGS_log")

  # relevel norm_counts after dropping genes
  # Payman
  str(norm_counts)
  norm_counts$Gene <- as.factor(norm_counts$Gene)
  norm_counts$Gene <- droplevels(norm_counts$Gene)

  return(norm_counts)
}

remove_blacklist_genes <- function(dfs, rm_genes){
  w2rm <- which(dfs[[2]]$Gene %in% rm_genes)
  if (length(w2rm) > 0){
    dfs[[2]] <- dfs[[2]][-w2rm,]
  }

  return(dfs[[2]])
}

change_LOQ <- function(dfs, LOQ_level = 2.5){
  samp_notes <- dfs[[3]]

  pools <- unique(dfs[[5]]$Pooling)

  for (pool in pools){
    name <- paste0("GeoLOQ", LOQ_level, "_", pool)
    samp_notes[[name]] <- data.frame(samp_notes[[paste0("NegGeoMean_", pool)]] *
                                       (samp_notes[[paste0("NegGeoSD_", pool)]]^LOQ_level))
    colnames(samp_notes[[name]]) <- name
  }


  return(samp_notes)
}

drop_expected_negs <- function(df, neg_flag = 1){
  raw_counts <- df[[1]]
  norm_counts <- df[[2]]
  samp_notes <- df[[3]]

  negs <- samp_notes$Sample_ID[samp_notes$expected_neg == neg_flag]

  raw_drop <- which(colnames(raw_counts) %in% negs)
  norm_drop <- which(colnames(norm_counts) %in% negs)
  samp_drop <- which(samp_notes$Sample_ID %in% negs)

  df[[1]] <- raw_counts[,-raw_drop]
  df[[2]] <- norm_counts[,-norm_drop]
  df[[3]] <- samp_notes[-samp_drop,]

  return(df)
}

#' @title Function that drops samples that meet a given criteria
#' @param df list of data.frame objects
#' @param cols what column(s) can criteria be found in Segment Properties
#' @param conditions criteria to find in segment properties
#' @param keep should samples that match conditions be kept
#' @param and if multiple conditions, should samples match both conditions or just one
#' @param custom custom filtering equation
#' @return list of data.frame objects with samples dropped
#' @examples
#' drop_samples(cols = "tissue", conditions = "CRC7", keep = TRUE, and = TRUE, custom = NULL)

drop_samples <- function(df, cols, conditions, operators, keep = FALSE, and = TRUE, custom){
  if(is.null(custom)){
    drop <- NULL
    values <- 0
    for(value in 1:length(cols)){
      if(! cols[value] %in% colnames(df[[3]])){
        print(paste(cols[value], "is not a valid column; skipping matching for this pairing"))
      }else{
        if(operators[value] == "!="){
          drop <- c(drop, which(df[[3]][,cols[value]] != conditions[value]))
        }else if(class(df[[3]][[cols[value]]]) == "factor" | class(df[[3]][[cols[value]]]) == "character" |
                 operators[value] == "=="){
          drop <- c(drop, which(df[[3]][,cols[value]] == conditions[value]))
        }else if(operators[value] == ">"){
          drop <- c(drop, which(df[[3]][,cols[value]] > as.numeric(conditions[value])))
        }else if(operators[value] == ">="){
          drop <- c(drop, which(df[[3]][,cols[value]] >= as.numeric(conditions[value])))
        }else if(operators[value] == "<"){
          drop <- c(drop, which(df[[3]][,cols[value]] < as.numeric(conditions[value])))
        }else if(operators[value] == "<="){
          drop <- c(drop, which(df[[3]][,cols[value]] > as.numeric(conditions[value])))
        }
        values <- values + 1
      }
    }

    if(values > 1){
      if(and){
        drop <- as.numeric(names(table(drop)[table(drop) == values]))
      }else{
        drop <- unique(drop)
      }
    }
  }else{
    drop <- which(eval(parse(text = custom)))
  }


  if(keep){
    neg_flag <- 0
    norm <- 1
  }else{
    neg_flag <- 1
    norm <- 0
  }

  df[[3]]$expected_neg <- norm
  if(length(drop) > 0){df[[3]]$expected_neg[drop] <- neg_flag}

  num_drops <- length(which(df[[3]]$expected_neg == 1))

  df <- drop_expected_negs(df,neg_flag = 1)

  print(paste(num_drops, "samples were dropped;", nrow(df[[3]]), "samples kept"))
  flog.info(paste(num_drops, "samples were dropped by user;", nrow(df[[3]]), "samples kept"),
            name="DSP_NGS_log")

  #drop unused factor levels from segment notes
  df[[3]] <- droplevels(df[[3]])

  return(df)
}
