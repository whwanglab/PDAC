###############
### AZORIUS ###
###############

# Last Updated: January 7, 2021 by Nicole Ortogero

##############################
#### Variables to Assign #####
##############################

#
# BEFORE DOING ANYTHING: Set your working directory to where your annotation sheet,
# deduplicated count table, and summary table are. This will be where your tables
# and figures are output as well.
# https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces 
#

# Must have RTS_ID, Gene and Module columns
# Set to NULL if using PKC file(s)
lookup_table_for_rnaids <- '../azorius/cta2.0_probelookup.txt'

# What do you want your figures labeled with?
# Files will be saved with this prefix but all spaces replaced with underscores
# Figures will use this in their titles
experiment_name <- 'Broad_PDAC_WTA_AllSamples'
annotations_file <- '../All_annotationsMerged.xlsx' # Sample Annotation File
annotations_sheet_name <- 'Annotations' # Sample Annotation Sheet

# DCC Files input
# Set both to NULL if just using dedup & summary files
dcc_directory <- './All_DCCs/'
output_prefix <- 'Broad_PDAC_WTA_AllSamples' #prefix for output summary and dedup files
n_processors <- 1

# Summary and dedup files, optional
# Will get overwritten if using dcc folder
dedup_file <- '../MK017_dedup.txt' # Deduplicated Counts table
processing_summary_file <- '' # Summary Table

# Target Groups
# This is a new (optional) file as of 12/03/2019
# This file contains two columns. One with probe groups and one with gene names.
# You get this file from Erin or at P:\DSP-NGS\KiloPlex\Target_Groups_Tables
# This file is important if you would like to do DE in Hydra
# Set to NULL if using PKC file(s)
targetGroups_file <- 'CTA2.0_probe_groups_short.txt'

# PKC Files
# Overrides use of targetGroups_file and lookup_table_for_rnaids
# Leave vector empty if using targetGroups_file and lookup_table_for_rnaids
# List PKC file(s)' path(s) in relation to working directory
pkc_vec <- c('Dev_WTx_v1.6.pkc')

# Analysis Type
# String indicating type of Analysis
# Determines what type of outlier testing to perform
# Allowed options "CTA" or "WTA"
analysis_type = "WTA"

# Flag that identifies an expected negative sample in the
# expected_neg column of the annotations file
# Michelle usually uses 1, Zoey/Kristina use 0
expect_neg_flag = 1

# Sample Name Columns (Used to Identify PCR Duplicates)
# list them exactly as they appear in the sample annotations file
# The columns will be pasted together with '_' and used as the
# Sample_ID after the initial QC steps. 
samp_nm_cols <- c("Sample_ID")

# Desired columns from the annotations file
# These are the annotation columns you would like included in
# the 2 main output tables
# Set to NULL if you would like to keep all annotation columns
annot_cols <- c('Batch','Patient','PatientNumber','Slide_name','Scan_name',
                'ROI_number','X_coordinate','Y_coordinate','Segment','AOI_area',
                'ROI_mask_area','Not_malignant','Not_collected',
                'Nuclei','Histology','Infiltration','PrevDropped',
                'TreatmentClass','Treatment','Status')

# Figure Generation Options
skipgpdf <- TRUE # Set to true to skip making gene pdf

no.figures <- TRUE # Set to true to skip making all figures

figs.byset <- FALSE
setsize <- 48

# Filter & Stats Options
remove_dropouts = TRUE
percent_unique_filter = 70
read_minimum = 10000 # Samples with fewer than this many reads will be removed from the analysis
low_count_filter = 0.1 # Ratio of mean(probe count) / mean(gene count)
minimum_count = 10 # Minimum count used for outlier testing
local_alpha = 0.01 # alpha for local outlier testing
freq_outlier_thresh = 0.2 # Portion of AOIs a probe needs to be identified as an outlier in before being flagged as a global outlier
remove_locals = TRUE # whether or not local outlier probes (within a given gene in a single AOI) should be removed
local_both_ways = TRUE # boolean indicating to perform local outlier testing both probe- and sample-wise
minLOQ = 2 # set a minimum LOQ, any calculated LOQ below this value will be reset to the minimum

# Housekeeper list
# If using pkc files for probe information and you would like to 
# use the pkc file designated housekeepers, set this to an empty vector
housekeepers <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
                  'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')

# This tells it not to "flip underscores to dashes" in the sample annotation sheet. Switch to FALSE if using data from UW.
flip_annotation_underscores = T
#
# These are the original 12 housekeepers selected by Margaret + Michelle. 
# If you want to reset the housekeepers listed above to the original just
# copy this string up
# housekeepers <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
# 'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')
#
##############################################################################################################################
########################################################## STOP ##############################################################
##############################################################################################################################
#
##############################################################################################################################
###################################################### Don't Touch ###########################################################
##############################################################################################################################
#
#
#
# If you have a feature request or encounter a bug, let Michelle know.
#
#
#
#
#
#### Required packages ####
required_packages <- c('reshape2', 'ggplot2', 'scales',
                       'dplyr', 'stringr', 'readxl',
                       'EnvStats', 'outliers', 'stringi',
                       'purrr', 'RColorBrewer', 'xlsx',
                       'vioplot', 'rjson', 'plyr', 'parallel', 
                       'data.table')
for (pkg in required_packages) {
  if (!pkg %in% installed.packages()) install.packages(pkg) 
}

lapply(required_packages, require, character.only=TRUE)

#### Functions ####
# generate_pkc_lookup
# Function to generate lookup table from PKC file format
# INPUT
#   jsons_vec = vector of json target and probe properties with pools as names
# OUTPUT
#   data frame with same format as probe lookup table
generate_pkc_lookup <- function(jsons_vec) {
  lookup_df <- data.frame(RTS_ID=character(), 
                            Gene=character(), 
                            Module=character(), 
                            Codeclass=character(),
                            stringsAsFactors=FALSE)
  for (curr_idx in 1:length(jsons_vec)) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    for (targ in curr_json[["Targets"]]) {
      curr_gene <- targ[["DisplayName"]]
      curr_code_class <- targ[["CodeClass"]]
      for (prb in targ[["Probes"]]) {
        curr_RTS_ID <- prb$RTS_ID
        lookup_df[nrow(lookup_df) + 1, ] <- 
          list(curr_RTS_ID, curr_gene, curr_module, curr_code_class)
      }
    }
  }
  # Grab the suffix digits from codeclass
  lookup_df[["PoolNum"]] <- 
    substr(lookup_df[["Codeclass"]], 
             nchar(lookup_df[["Codeclass"]]) - 1, 
             nchar(lookup_df[["Codeclass"]]))
  # Add set suffix to pool name in case more than one set in pkc
  lookup_df[["Module"]] <- 
    paste(lookup_df[["Module"]], lookup_df[["PoolNum"]], sep="_")
  return(lookup_df)
}

# generate_pkc_targ_notes
# Function to generate target groups table from PKC file format
# INPUT
#   jsons_vec = vector of json target and probe properties with pools as names
# OUTPUT
#   data frame with same format as target groups file
generate_pkc_targ_notes <- function(jsons_vec, lookup_tab) {
  # Create non-duplicated map from gene to pool and codeclass
  sub_lookup <- unique(rnaid_lookup_df[, names(rnaid_lookup_df) != "RTS_ID"])
  rownames(sub_lookup) <- sub_lookup[["Gene"]]
  notes_df <- 
    data.frame(TargetName=rownames(sub_lookup),
                 HUGOSymbol=rownames(sub_lookup),
                 TargetGroup=rep("All Probes", length(rownames(sub_lookup))),
                 AnalyteType=rep("RNA", length(rownames(sub_lookup))),
                 Codeclass=sub_lookup[rownames(sub_lookup), "Codeclass"],
                 Pooling=sub_lookup[rownames(sub_lookup), "Module"],
                 stringsAsFactors=FALSE)
  for (curr_idx in 1:length(jsons_vec)) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    if(length(curr_json[["ProbeGroups"]]) > 0) {
      for (prb_group in curr_json[["ProbeGroups"]]) {
        curr_group <- prb_group[["Name"]]
        for (targ in prb_group[["Targets"]]) {
          notes_df[notes_df$TargetName == targ, "TargetGroup"] <-
            paste(notes_df[notes_df$TargetName == targ, "TargetGroup"], 
                    curr_group, sep=";")
        }
      }
    }
  }

  return(notes_df)
}

# pnum
# Function to make pretty numbers
# INPUT
#   n = number
# OUTPUT
#   number as a string with commas for hundred/thousands/etc
pnum <- function(n) {
  return(format(n, big.mark = ','))
}

# ngeoMean
# Robust geoMean Function
# INPUT
#   v = numeric vector
# OUTPUT
#   geometric mean with "NA"s removed
ngeoMean <- function(v) {
  v[v == 0] <- 1
  return(geoMean(v, na.rm = T))
}

# ngeoSD
# Robust geoSD Function
# INPUT
#   v = numeric vector
# OUTPUT
#   geometric standard deviation with "NA"s removed
ngeoSD <- function(v) {
  v[v == 0] <- 1
  return(geoSD(v, na.rm=T))
}

# rosner.flag 
# function to remove outliers from vector and return the names
# of removed outliers
# INPUT
#   x = named vector
#   p = percent of observations suspected of being outliers
#   alpha = indicator of the expected Type I erorr frequency
#   min_count = integer of minimum expected counts for testing
# OUTPUT
#   vector of TRUE / FALSE For flagged or not flagged
rosner.flag <- function(x, p=0.2, alpha_thresh=0.01, min_count=2) {
  returnvec <- rep("", length(x))
  # Skip analysis if vector is all the same value or below min count threshold
  if (all(x == x[1]) || max(x) < min_count) {
    return(returnvec)
  }
  
  flag <- rosnerTest(x, k = round(length(x)*p), alpha = alpha_thresh, warn = F)
  flag.idx <- subset(flag$all.stats, Outlier, select=Obs.Num, drop=T)
  if ( length(flag.idx) == 0) {
    return(returnvec)
  } else {
    # Outliers should lie to either side of the median with enough n and p < 1
    for (i in 1:length(flag.idx)) {
      if(x[flag.idx[i]] < median(x)) {
        returnvec[flag.idx[i]] <- "low"
      } else if (x[flag.idx[i]] > median(x)) {
        returnvec[flag.idx[i]] <- "high"
      }
    }
    return(returnvec) 
  }
}

# grubbs.flag
# helper function to remove outliers using Grubbs' test given the controlled type I error alpha
# modified from https://stackoverflow.com/questions/22837099/how-to-repeat-the-grubbs-test-and-flag-the-outliers
# INPUT
#   x = named vector
#   alpha = indicator of the expected Type I erorr frequency
#   logt = boolean to log 10 transform data prior to outlier test
#   min_count = integer of minimum expected counts for testing
# OUTPUT
#   vector of TRUE / FALSE For flagged or not flagged
grubbs.flag <- function(x, alpha_thresh=0.01, logt=TRUE, min_count=2) {
  outliers <- rep("", length(x))
  # Skip analysis if vector is all the same value or below min count threshold
  if (all(x == x[1]) || max(x) < min_count) {
    return(outliers)
  }
  if(logt == TRUE){
    initial <- log10(pmax(x, 1))
  } else {
    initial <- x
  }

  test <- sort(initial)
  grubbs.result <- outliers::grubbs.test(test, two.sided=TRUE)
  if(grubbs.result$p.value < alpha_thresh) {
    if(grepl("lowest", grubbs.result$alternative)){
      outliers[[1]] <- "low"
    } 
    if(grepl("highest", grubbs.result$alternative)){
      outliers[[length(test)]] <- "high"
    } 
    names(outliers) <- names(test)
    # reorder outliers to match initial
    outliers <- outliers[order(factor(names(outliers), levels=names(initial)))]
  }

  names(outliers) <- NULL
  return(outliers)
}

# transpose_local
# test that local outlier is also an outlier across all AOIs
# for a more conservative exclusion of local outliers
# INPUT
#   test_df = dataframe of counts from genes flagged as local outliers
#   locals = list of sample probes flagged as local outliers
# OUTPUT
#   list of local outliers that were flagged across sample IDs
transpose_local <- function(test_df, local_lows, local_highs, min_count=2, alpha_thresh=0.01) {
  # Do not perform test if sample size is too small, keep all local outliers
  # Also skip if there are no local outliers
  if (length(unique(test_df[["Sample_ID"]])) < 3 || length(c(local_lows, local_highs)) == 0) {
    return (c())
  }

  # Get names of genes with local outliers
  local_genes <-
    unique(c(as.character(test_df[test_df[["Outliers"]] == "Local Low Outlier", "Gene"]),
             as.character(test_df[test_df[["Outliers"]] == "Local High Outlier", "Gene"])))
  # Reduce size of test_df for gene subsetting
  test_df <- test_df[test_df[["Gene"]] %in% local_genes, ]
  # Filter out low count and global outlier probes
  test_df <- test_df[test_df[["LowFlag"]] != "Low Dropout", ,drop=FALSE]

  # Generate lists of local outliers found across sample IDs
  keep_lists <- 
    lapply(
      local_genes, 
      function(gene, full_df=test_df) {
        # Munge gene data for testing
        gene_df <- subset(full_df, Gene == gene)
        gene_df <- gene_df[,c("Sample_ID", "RNAID", "Count")]
        gene_df <- dcast(gene_df, Sample_ID ~ RNAID, value.var="Count")
        rownames(gene_df) <- gene_df[["Sample_ID"]]

        # Rosner Test with at least 10 sample IDs
        if (nrow(gene_df) >= 10) {
          outliercalls <- apply(gene_df[2:ncol(gene_df)], 2, rosner.flag,
                                  min_count=min_count, alpha_thresh=alpha_thresh)
          rownames(outliercalls) <- rownames(gene_df)
        # Grubbs Test with at least 3 sample IDs
        } else if (nrow(gene_df) > 2) {
          outliercalls <- apply(gene_df[2:ncol(gene_df)], 2, grubbs.flag,
                                  min_count=min_count, alpha_thresh=alpha_thresh)
          rownames(outliercalls) <- rownames(gene_df)
        }

        # Assign string value of local outlier to any TRUE cells
        outliercalls <- melt(outliercalls)
        lows <- 
          apply(outliercalls[which(outliercalls$value == "low"), c('Var1', 'Var2')],
                  1, paste0, collapse="")
        highs <- 
          apply(outliercalls[which(outliercalls$value == "high"), c('Var1', 'Var2')],
                  1, paste0, collapse="")
        # Get list of outliers that are local outliers both ways
        lows_to_keep <- local_lows[which(local_lows %in% lows)]
        highs_to_keep <- local_highs[which(local_highs %in% highs)]
        names(lows_to_keep) <- c()
        names(highs_to_keep) <- c()
        curr_keep <- c(lows_to_keep, highs_to_keep)
        return(curr_keep)
      })
  keep_local_tags <- unlist(keep_lists)
  return(keep_local_tags)
}

# INPUT: Dedup Counts data frame, lookupdf, ratio threshold
prbflag <- function(df, 
                      lookupdf, 
                      rt_thresh, 
                      aoi_thresh, 
                      min_count=2, 
                      alpha_thresh = 0.01,
                      expt_type="CTA") {
  # Initialize empty columns and list
  df$LowFlag <- ''
  df$Outliers <- ''
  globalOutliers <- c()
  localLows <- c()
  localHighs <- c()
  lowDrops <- c()

  # Get all target names
  gene_list <- unique(df$Gene)
  # Perform outlier test on negative probes only
  if (expt_type == "WTA") {
    gene_list <- gene_list[gene_list %in% df[df$Negative, "Gene"]]
  }

  for (gene in gene_list) {
    gene_df <- subset(df, Gene == gene)[,c('Sample_ID', 'RNAID', 'Count')]
   
    # Skip if less than 3 gene probes
    if (nrow(gene_df) < 3) { next }

    gene_df <- dcast(gene_df, RNAID ~ Sample_ID, value.var='Count')
    rownames(gene_df) <- gene_df$RNAID

    # Drop probes with low counts
    gene_mean <- ngeoMean(as.vector(as.matrix(gene_df[,2:ncol(gene_df)])))
    lowflag <- (apply(gene_df[,2:ncol(gene_df)], 1, ngeoMean) / 
                  gene_mean < rt_thresh)
    # Flag low count probes
    lowDrops <- append(lowDrops, names(which(lowflag)))
    # REMOVE LOW FLAG PROBES
    gene_df <- gene_df[!lowflag,]
    
    # Move on if less than 3 probes remaining
    if (nrow(gene_df) < 3) { next }

    # Rosner Test with at least 10 probes
    if (nrow(gene_df) >= 10) {
      outliercalls <- apply(gene_df[2:ncol(gene_df)], 2, rosner.flag, 
                              min_count=min_count, alpha_thresh=alpha_thresh)
      rownames(outliercalls) <- rownames(gene_df)
      outlierfreq <- apply(outliercalls, 1, function(prbflags) {
        return( length(which(prbflags != "")) / ncol(outliercalls))
      })
      globalOutliers <- append(globalOutliers, 
                                names(which(outlierfreq > aoi_thresh)))
    # Grubbs Test with at least 3 probes
    } else if (nrow(gene_df) > 2) {
      outliercalls <- apply(gene_df[2:ncol(gene_df)], 2, grubbs.flag, 
                              min_count=min_count, alpha_thresh=alpha_thresh)
      rownames(outliercalls) <- rownames(gene_df)
      outlierfreq <- apply(outliercalls, 1, function(prbflags) {
        return( length(which(prbflags != "")) / ncol(outliercalls))
      })
      globalOutliers <- append(globalOutliers, 
                                names(which(outlierfreq > aoi_thresh)))
    }

    # Assign string value of local outlier to any TRUE cells
    outliercalls <- melt(outliercalls)
    localLows <- 
      append(
        localLows, 
        apply(outliercalls[which(outliercalls$value == "low"), c('Var2', 'Var1')],
                1,
                paste0,
                collapse=''))
    localHighs <- 
      append(
        localHighs, 
        apply(outliercalls[which(outliercalls$value == "high"), c('Var2', 'Var1')],
                1,
                paste0,
                collapse=''))
  }
  names(localLows) <- c()
  names(localHighs) <- c()

  # Add flags to returning dataframe
  df$LowFlag[which(df$RNAID %in% lowDrops)] <- "Low Dropout"
  df$Outliers[which(paste0(df$Sample_ID, df$RNAID) %in% localLows)] <- 
    'Local Low Outlier'
  df$Outliers[which(paste0(df$Sample_ID, df$RNAID) %in% localHighs)] <- 
    'Local High Outlier'
  if (length(globalOutliers) != 0) {
    df$Outliers[which(df$RNAID %in% globalOutliers)] <- 'Global Outlier'
  }

  # Unflag local outliers not outliers across AOIs
  if(local_both_ways)
  {
    # Get list of outliers that are local outliers both along AOI and probe-wise
    keep_locals <- transpose_local(test_df=df, 
                                     local_lows=localLows,
                                     local_highs=localHighs,
                                     min_count=min_count,
                                     alpha_thresh=alpha_thresh)
    # Remove local outlier tag for outliers not found both ways
    if (length(keep_locals) != 0) {
      # Change flag to indicate this is two way outlier
      df$Outliers[
        intersect(which(paste0(df$Sample_ID, df$RNAID) %in% keep_locals),
                    which(df$Outliers == "Local Low Outlier"))] <- "Local Outlier"
      df$Outliers[
        intersect(which(paste0(df$Sample_ID, df$RNAID) %in% keep_locals),
                    which(df$Outliers == "Local High Outlier"))] <- "Local Outlier"
    }
    # Remove flag if local outlier only one way and still has old flag
    df[df[["Outliers"]] == "Local Low Outlier", "Outliers"] <- ""
    df[df[["Outliers"]] == "Local High Outlier", "Outliers"] <- ""
  } else {
    # Rename for downstream filtering
    df[df[["Outliers"]] == "Local Low Outlier", "Outliers"] <- "Local Outlier"
    df[df[["Outliers"]] == "Local High Outlier", "Outliers"] <- "Local Outlier"
  }

  return(df)
}

#### Check user inputs ####
if (!analysis_type %in% c("CTA", "WTA")) {
  stop(paste("Invalid analysis type designated.", 
               "Only \"CTA\" and \"WTA\" allowed as inputs."))
}

if(!is.null(targetGroups_file)) {
  if(!file.exists(targetGroups_file)) {
    stop("Target groups file designated, but does not exist. Check filepath. Aborting.")
  }
}

if(length(pkc_vec)>0) {
  if(length(pkc_vec) != length(Sys.glob(pkc_vec))) {
    stop("Some PKC file(s) do not exist. Check filepath. Aborting.")
  }
}

# Process DCC files if needed
if(!is.null(dcc_directory)){
  #copy (with trimming) of dnd_summary.R written by Tyler Hether 2020-10-12
  
  ### ##########################
  ### Check input path and files
  ### ##########################
  
  # Add a trailing / to directory name
  if(!endsWith(dcc_directory, "/")){
    dcc_directory <- paste0(dcc_directory, "/")
  }
  # Confirm that directory exists
  if(!dir.exists(dcc_directory)){
    stop(paste0("The input directory, ", dcc_directory, 
                " was not found. Please check path."))
  }
  # Get the dcc files in dcc_directory
  dcc_files <- Sys.glob(paste0(dcc_directory, "*dcc"))
  # Confirm that there is at least 1 dcc file in dcc_directory
  if(length(dcc_files)<1){
    stop(paste0("The input directory, ", dcc_directory, 
                " was found but no files ending with .dcc were located.\n", 
                "Are these files compressed?"))
  }
  
  ### #########
  ### Functions
  ### #########
  
  process_dcc <- function(the_dcc){
    # ARGS:
    # the_dcc is a character with a valid dcc file path.
    # 
    # returns a list of length two providing the parsed data
    
    # Here are the lines
    the_lines <- readLines(the_dcc)
    
    # Here is the Sample_ID
    scan_attributes_header_line <- which(the_lines=="<Scan_Attributes>")
    if(length(scan_attributes_header_line)!=1){
      stop(paste0("The Scan_Attributes header was not found in ", the_dcc))
    }
    if(grepl("^ID,", the_lines[scan_attributes_header_line+1])){
      Sample_ID <- strsplit(the_lines[scan_attributes_header_line+1], split="ID,")[[1]][2]
    } else {
      stop(paste0("The Sample_ID line was not found in ", the_dcc))
    }
    
    # Raw, Trimmed, Stitched, and Aligned are bounded withing
    # the NGS_Processing_Attributes header
    # Here are the bounding indeces
    ngs_prcoessing_header_line <- which(the_lines=="<NGS_Processing_Attributes>")
    if(length(ngs_prcoessing_header_line)!=1){
      stop(paste0("The NGS_Processing_Attributes header was not found in ", the_dcc))
    }
    ngs_prcoessing_footer_line <- which(the_lines=="</NGS_Processing_Attributes>")
    if(length(ngs_prcoessing_footer_line)!=1){
      stop(paste0("The NGS_Processing_Attributes footer was not found in ", the_dcc))
    }
    # And parse the Raw, Trimmed, Stitched, and Aligned
    rtsa_names <- c("Raw", "Trimmed", "Stitched", "Aligned")
    rtsa <- lapply(rtsa_names, function(x){
      x_line <- which(grepl(paste0("^",eval(x), ","), 
                            the_lines[ngs_prcoessing_header_line:ngs_prcoessing_footer_line]))
      if(length(x_line)!=1){
        stop(paste0("The ", x, " counts were not found in file ", the_dcc))
      } else {
        x_out <- as.numeric(
          strsplit(the_lines[ngs_prcoessing_header_line:ngs_prcoessing_footer_line][x_line], 
                   split=paste0(eval(x), ","))[[1]][2])
        return(x_out)
      }
    })
    rtsa_df <- as.data.frame(do.call(cbind, rtsa))
    colnames(rtsa_df) <- rtsa_names
    # Make summary_df by combining data (missing Uniques at this point)
    summary_df <- cbind(data.frame(Sample_ID=Sample_ID), rtsa_df)
    
    # Now pull out the RTS IDs and their associated counts.
    code_summary_header_line <- which(the_lines=="<Code_Summary>")
    if(length(code_summary_header_line)!=1){
      stop(paste0("The Code_Summary header was not found in ", the_dcc))
    }
    code_summary_footer_line <- which(the_lines=="</Code_Summary>")
    if(length(code_summary_footer_line)!=1){
      stop(paste0("The Code_Summary footer was not found in ", the_dcc))
    } 
    # If there are not entries, provide a warning
    # example:
    # <Code_Summary>
    # </Code_Summary>
    uniques <- 0
    if(code_summary_footer_line - code_summary_header_line < 2){
      warning(paste0("There are no entries in code summary for file ", the_dcc, 
                     " so all RTS values will be 0."))
      RTSs_df <- data.frame(Sample_ID=Sample_ID)
    } else {
      # Parse and transform
      RTSs <- the_lines[(code_summary_header_line+1):(code_summary_footer_line-1)]
      RTSs_long <- do.call(rbind, strsplit(RTSs, split=","))
      RTSs_wide <- matrix(as.numeric(RTSs_long[,2]), nrow=1, dimnames=list(NULL, RTSs_long[,1]))
      uniques <- apply(RTSs_wide, 1, sum)
      RTSs_df <- data.frame(Sample_ID=Sample_ID, as.data.frame(RTSs_wide))
    }
    
    # Add the uniqes to summary_df
    summary_df$Unique <- uniques
    
    # Return the two data.frame objects
    return(list(summary_df, RTSs_df))
    
  }
  
  ### ##########
  ### Processing
  ### ##########
  
  # Two steps: 1) Process each file 2) merge data
  
  ## Step 1 process each file
  # Set up cluster
  cl <- makeCluster(n_processors)
  clusterExport(cl=cl, varlist=c("dcc_files", "process_dcc"), envir=environment())
  # Execute the main function in parallel
  all_processed <- parLapply(cl, dcc_files, process_dcc)
  # Stop cluster
  stopCluster(cl)
  
  ## Step 2 merge data
  # The summary file.
  summary_out <- do.call(rbind, lapply(all_processed, "[[", 1L))
  
  # The dedup file. 
  # Note that the RTS IDs need to be in the same order so
  # we will use dplyr's bind_rows
  dedup_out <- bind_rows(lapply(all_processed, "[[", 2L))
  
  # Ensure there are no conversion issues with Sample_IDs.
  if(any(is.na(dedup_out$Sample_ID))){
    stop("Sample_ID was converted to NA, implying an error in processing.")
  }
  
  # Let the user know that how many NA counts are being converted to zero and convert.
  message(paste0("Converting ", length(which(is.na(dedup_out))), " NAs to zero."))
  dedup_out[is.na(dedup_out)] <- 0
  
  # reorder the column names
  cols_sorted <- c("Sample_ID", sort(colnames(dedup_out)[2:ncol(dedup_out)]))
  dedup_out <- dedup_out %>% select(eval(cols_sorted))
  
  # Why convert NA to zero? 
  # When columns do not exist, bind_rows places an NA in the position.
  # x <- data.frame("Sample_ID"="one", "RTS1"=1, "RTS2"=2)
  # y <- data.frame("Sample_ID"="two", "RTS3"=3, "RTS2"=2.1, "RTS4"=4)
  # z <- bind_rows(x,y)
  # z
  # z[is.na(z)] <- 0
  # z
  
  ### ##########
  ### Write data
  ### ##########
  
  processing_summary_file <- paste0(output_prefix, "_summary.txt")
  dedup_file <- paste0(output_prefix, "_dedup.txt")

  data.table::fwrite(summary_out, file = processing_summary_file, 
                     sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  data.table::fwrite(x=dedup_out, file = dedup_file, 
                     sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=0)
  
}

#### Read in Files ####
# Output Root File Name
fRoot = gsub(' |/', '_', experiment_name)

# Import Sequencing QC Info
processing_sum <- read.delim(processing_summary_file)
processing_sum$Percent_Unique <- 100*(processing_sum$Unique / 
                                        processing_sum$Aligned)

# Import Sample Annotations
annotations <- read_excel(annotations_file, sheet=annotations_sheet_name)
if (flip_annotation_underscores) {
  annotations$Sample_ID <- gsub('_', '-', annotations$Sample_ID)
}

if (any(!sapply(annot_cols, function(colnm) {colnm %in% colnames(annotations)}))) {
  
  stop(paste('Columns listed in "annot_cols" which are not found in the annotations file:',
       paste(annot_cols[which(!(annot_cols %in% colnames(annotations)))], collapse=', ')
       ))
}

#### Sequencing Depth QC ####
# Summary Statistics Blurb
sumstats <- paste0("Total Aligned Reads: ", round(sum(processing_sum$Aligned) /
                                                    1000000), 'M\n',
                   "Total Unique Reads: ", 
                   round(sum(processing_sum$Unique) / 1000000), 'M (',
                   round((sum(processing_sum$Unique) / 
                            sum(processing_sum$Aligned))*100 ), 
                   "%)\nAligned Median: ", 
                   pnum(median(processing_sum$Aligned)),
                   "\nAligned Range: ", pnum(min(processing_sum$Aligned)), 
                   ' - ', pnum(max(processing_sum$Aligned)))

# Identify outlier samples based on % dedup
dropouts <- as.character(processing_sum$Sample_ID[
  processing_sum$Percent_Unique > percent_unique_filter])

text_summary <- paste0("Summary for ", experiment_name, 
                       '\n', nrow(annotations), ' AOIs + NTCs Expected',
                       '\n', sumstats, "\nPercent Unique > ", 
                       percent_unique_filter, "% Dropouts: ", 
                       paste(dropouts, collapse=', '), '\n')

low_reads <- as.character(processing_sum$Sample_ID[
  processing_sum$Raw < read_minimum ])

text_summary <- paste0(text_summary, "Low Read Dropouts: ",
                       paste(low_reads, collapse=','),
                       '\n')

dropouts <- append(dropouts, low_reads)

# Also flag any expected negatives for removal
expected_neg = annotations$Sample_ID[annotations$expected_neg == 
                                       expect_neg_flag]
text_summary <- paste0(text_summary, "Expected negatives to be removed: ",
                       paste(expected_neg, collapse = ', '), '\n')
dropouts <- append(dropouts, expected_neg)

# Subset Annotations to just desired columns
if(!is.null(annot_cols)) {
  annotations <- annotations[,unique(c('Sample_ID', annot_cols))]
}
tmp <- processing_sum[match(annotations$Sample_ID, processing_sum$Sample_ID),]
tmp$Percent_Unique <- 100 - tmp$Percent_Unique
colnames(tmp) <- c('Sample_ID', 'RawReads', 'TrimmedReads', 'StitchedReads',
                   'AlignedReads', 'DeduplicatedReads', 'SequencingSaturation')

annotations <- left_join(annotations, tmp, by=c('Sample_ID'))
rm(tmp)

# Data frame for Plotting
processing_sum <- melt(processing_sum[,c('Sample_ID','Aligned','Unique')], 
                       id.vars = 'Sample_ID')

# Set up figure sets if desired
samplist <- unique(processing_sum$Sample_ID)
if (figs.byset) {
  samp_count <- length(samplist)
  figsets <- list()
  for (setnum in 1:ceiling(samp_count/setsize)) {
    figsets[[paste0('Set', formatC(setnum, digits=2, flag="0"))]] <- 
      samplist[((setnum-1)*setsize + 1):min((setnum*setsize),samp_count)]
  }
  
} else {
  figsets <- list('ALL'=samplist)
}

#### Median Aligned vs Unique Figure ####
# Plot Aligned and Unique for Each Sample
readcountplot <- function(sumdf, sampset, fname) {
  sumdf <- sumdf[sumdf$Sample_ID %in% sampset,]
  
  yaxis_max = max(sumdf$value) + max(sumdf$value)*0.05
  png(paste0(fname, '_read_depth.png'), height=6, width=max(nrow(sumdf)/6, 11),
      units='in', res=300)
  print(
    ggplot(data=subset(sumdf, variable=='Aligned'), 
           aes(x = Sample_ID, y = value, fill=variable)) + 
      geom_bar(stat = "identity") +
      geom_bar(data=subset(sumdf, variable=='Unique'), stat="identity") +
      scale_y_continuous(labels=comma, limits=c(0,yaxis_max), expand=c(0,0)) +
      labs(title=paste0("Read Depth of ", experiment_name), y = 'Reads',
           subtitle=sumstats, fill='Type') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust=0.5),
            legend.position = 'bottom',
            axis.text = element_text(color='black'))
  )
  dev.off()
}


if (!no.figures) {
  for (figset in names(figsets)) {
    readcountplot(processing_sum, figsets[[figset]], 
                  paste0(fRoot, '_', figset))
  }
}

# Recast summary data
processing_sum <- dcast(processing_sum, Sample_ID ~ variable)

# Import Deduplicated Counts
dedup_counts <- read.delim(dedup_file)
colnames(dedup_counts)[2:ncol(dedup_counts)] <- 
  gsub('_.*', '', colnames(dedup_counts))[2:ncol(dedup_counts)]

# Read in lookup table
if(length(pkc_vec) == 0) {
  rnaid_lookup_df <- read.delim(lookup_table_for_rnaids, stringsAsFactors=FALSE)
# Read in json file to fill lookup and target notes tables
} else {
  pkc_json_list <- lapply(pkc_vec, function(pkc_file) {fromJSON(file=pkc_file)})
  pkc_names <- 
    unlist(lapply( pkc_vec, 
             function(file) {
                               base_pkc_name <- gsub(".pkc", "", trimws(basename(file)))
                               return(base_pkc_name)
                             }))
  names(pkc_json_list) <- pkc_names
  rnaid_lookup_df <- generate_pkc_lookup(pkc_json_list)
  target_notes <- generate_pkc_targ_notes(pkc_json_list, rnaid_lookup_df)
  # Remove codeclass column, only needed for target notes generation
  rnaid_lookup_df <- 
    rnaid_lookup_df[, !colnames(rnaid_lookup_df) %in% c("Codeclass", "PoolNum")]
}

# Subset Dedup to just RTSIDs in lookuptable
dedup_counts <- dedup_counts[,c(1, 
                  which(colnames(dedup_counts) %in% rnaid_lookup_df$RTS_ID))]

# Determine if there are any probes with all zero counts (Probably not)
all_zero_tally <- length(which(sapply(dedup_counts[,2:ncol(dedup_counts)], 
                                      function(x) all(x == 0))))

probes_detected <- ncol(dedup_counts) - all_zero_tally - 1

text_summary <- paste0(text_summary, 'Probes Detected at least once: ', 
                       pnum(probes_detected), '\n')

text_summary <- paste0(text_summary, pnum(nrow(rnaid_lookup_df)), 
                       ' Probes listed in expected file.\n')

#### Report out NTC summary stats ####
if (length(expected_neg) > 0) {
  expect_neg_counts <- dedup_counts[dedup_counts$Sample_ID %in% expected_neg, ]
  rownames(expect_neg_counts) <- expect_neg_counts[["Sample_ID"]]
  expect_neg_counts <- expect_neg_counts[, colnames(expect_neg_counts) != "Sample_ID"]
  expect_neg_summ <-
    do.call(rbind, apply(expect_neg_counts, 1,
            function(x) {
              x <- as.numeric(x)
              x_summ <- list("Mean"=mean(x),
                                "Median"=median(x),
                                "Var"=var(x),
                                "Min"=min(x),
                                "Max"=max(x),
                                "Total"=sum(x))
            }))
  text_summary <- 
    paste0(text_summary, 
             "NTC Summary: \n\t",
             paste(colnames(expect_neg_summ), collapse="\t"), 
             "\n",
             paste0(
               apply(cbind(rownames(expect_neg_summ), expect_neg_summ),
                       1, paste0, collapse="\t"), 
               collapse="\n"),
             "\n")

  # Warn users if high counts in AOIs flagged as expected negative
  if (any(expect_neg_summ[,"Mean"] > 5)) {
    warning(
      paste0("High NTC counts in expected negatives. ",
               "The following NTCs exceed an average probe count of 5.\n  ",
               paste(
                 row.names(expect_neg_summ[expect_neg_summ[,"Mean"] > 5,]),
                 collapse="\n  ")))
  }
  if (any(expect_neg_counts > 100)) {
    warning(
      paste0("High NTC probe counts in expected negatives. ", 
               "The following NTCs have at least one probe with a count greater than 100.\n  ",
               paste(
                 names(apply(expect_neg_counts > 100, 1, sum) > 0),
                 collapse="\n  ")))
  }
}

#### Dedup Counts to Apply-able ####
# Remove %Unique & Expected Negative Dropouts
if (remove_dropouts) {
  text_summary <- paste0(text_summary, length(unique(dropouts)), 
                         ' AOIs / NTCs removed from analysis.\n')
  dedup_counts <- dedup_counts[!(dedup_counts$Sample_ID %in% dropouts),]
  processing_sum <- processing_sum[!(processing_sum$Sample_ID %in% dropouts),]
  annotations <- annotations[!(annotations$Sample_ID %in% dropouts),]
  
  text_summary <- paste0(text_summary,
    "Remaining Aligned Reads: ", round(sum(processing_sum$Aligned) /
                                         1000000), 'M',
    "\n\tRemaining Aligned Median: ", 
    pnum(median(processing_sum$Aligned)),
    "\n\tRemaining Aligned Range: ", pnum(min(processing_sum$Aligned)), 
    ' - ', pnum(max(processing_sum$Aligned)),
    "\nRemaining Unique Reads: ", 
    round(sum(processing_sum$Unique) / 1000000), 'M (',
    round((sum(processing_sum$Unique) / 
             sum(processing_sum$Aligned))*100 ), '%)\n')
}

if (nrow(dedup_counts) == 0) {
  stop('No samples passed initial filters. Execution halted.')
}

# Read in Annotations and add Gene information to dedup counts
dedup_counts <- melt(dedup_counts, id.vars = 'Sample_ID')
colnames(dedup_counts) <- c('Sample_ID', 'RNAID', 'Count')
dedup_counts$Gene <- rnaid_lookup_df$Gene[match(dedup_counts$RNAID, 
                                                rnaid_lookup_df$RTS_ID)]

if (any(is.na(dedup_counts$Gene))) {
  stop('Some RNAIDs were not found in the lookup table.')
}
#### STOP HERE
old_SID <- dedup_counts$Sample_ID
dedup_counts$Sample_ID <- gsub('_','-',dedup_counts$Sample_ID)

annotations$UID <- do.call(paste, c(annotations[,samp_nm_cols], sep='_'))
if (any(duplicated(annotations$UID))) {
  stop('Annotation columns selected to generate a sample ID do not generate',
       ' unique sample IDs.')
}

all(annotations$Sample_ID %in% dedup_counts$Sample_ID)
unique(dedup_counts$Sample_ID[!dedup_counts$Sample_ID %in% annotations$Sample_ID])
# samples not included in dedup file are NTCs and should be removed
# 
dedup_counts$Sample_ID <- annotations$UID[match(dedup_counts$Sample_ID, 
                                                    annotations$Sample_ID)]

if (any(is.na(dedup_counts$Sample_ID))) {
  warning(paste0("Not all Sample_IDs found in annotations sheet. ",
          "Removing unidentified data points."))
  dedup_counts <- dedup_counts[!is.na(dedup_counts$Sample_ID),]
}

## Saved WorkingAzorius - NTCs Removed but PP included.RData HERE

text_summary <- paste0(text_summary,
  pnum(length(unique(dedup_counts$Gene))),
  " targets found in dataset. ",
  pnum(length(unique(rnaid_lookup_df$Gene))),
  " targets listed in lookup table.\n"
)

# Add pooling and negative information
dedup_counts$Pool <- rnaid_lookup_df$Module[match(dedup_counts$RNAID, 
                                                  rnaid_lookup_df$RTS_ID)]
if(length(pkc_vec) == 0) {
  dedup_counts$Negative <- sapply(as.character(dedup_counts$Gene), startsWith, 
                                prefix='NegProbe')
} else {
  dedup_counts$Negative <- 
    dedup_counts$Gene %in% 
      target_notes[grep("Negative", target_notes$Codeclass), "TargetName"]
}

# Flag Low Count & Outlier Probes
dedup_counts <- prbflag(dedup_counts, 
                          rnaid_lookup_df, 
                          low_count_filter, 
                          freq_outlier_thresh, 
                          min_count=minimum_count, 
                          alpha_thresh=local_alpha,
                          expt_type=analysis_type)

# Separate out flagged counts
if (remove_locals == TRUE){
  flaggedRows <-c(which(dedup_counts$LowFlag != ''), 
                  which(dedup_counts$Outliers != ''))  
} else{
  # Always remove local outliers for negative probes
  flaggedRows <-c(which(dedup_counts$LowFlag != ''), 
                  which(dedup_counts$Outliers == 'Global Outlier'),
                  which(dedup_counts$Negative == TRUE &
                          dedup_counts$Outliers == 'Local Outlier'))
}
if (length(flaggedRows) != 0) {
  removed_probes <- dedup_counts[flaggedRows,]
  dedup_counts <- dedup_counts[-flaggedRows,]
} else {
  removed_probes <- dedup_counts[0, ]
}

text_summary <- paste0(text_summary, 
   pnum(length(unique(
     removed_probes$RNAID[which(removed_probes$LowFlag != '')]
     ))), ' low count probes identified.\n',
   pnum(length(unique(
     removed_probes$RNAID[which(removed_probes$Outliers == 'Global Outlier')]
     ))), ' probes identified as global grubbs outliers.\n'
   )

#### Raw Count Figure ####
# Plot all probe counts in black
# Plot geomean of negative probes / samp in red
# Plot geomean of experimental probes in green
# Subtitle as geomean summary stats
rawProbeplot <- function(dedupdf, sampset, fname) {
  dedupdf <- dedupdf[dedupdf$Sample_ID %in% sampset,]
  
  subtext <- paste0(
    "Range of Counts: ", pnum(min(dedupdf$Count)), " - ",
    pnum(max(dedupdf$Count)),
    "\nGeoMean Counts for Targets: ", 
    pnum(ngeoMean(dedupdf$Count[!dedupdf$Negative])),
    "\nGeoMean Counts for Negatives: ",
    pnum(ngeoMean(dedupdf$Count[dedupdf$Negative]))
  )
  
  # Create Geometric Mean Data Frame
  targetGeoMeans <- dcast(subset(dedupdf, !Negative, 
                                 select=c(Sample_ID, Pool, Count)),
                          Sample_ID ~ Pool, value.var='Count', 
                          fun.aggregate = ngeoMean)
  negativeGeoMeans <- dcast(subset(dedupdf, Negative, 
                                   select=c(Sample_ID, Pool, Count)),
                            Sample_ID ~ Pool, value.var='Count', 
                            fun.aggregate = ngeoMean)
  colnames(negativeGeoMeans)[2:ncol(negativeGeoMeans)] <- (
    paste0('Neg-', colnames(negativeGeoMeans)[2:ncol(negativeGeoMeans)])
  )
  
  geoMeanDf <- rbind(melt(targetGeoMeans), melt(negativeGeoMeans))
  col_num <- length(unique(geoMeanDf$variable))
  
  png(paste0(fname, '_counts_raw.png'), height=6, 
      width=max(length(sampset)/6, 11), units='in', res=300)
  
  rawplot <- ggplot(dedupdf, aes(x=Sample_ID, y=log2(Count + 1))) + 
    geom_jitter(width=0.1, height=0) +
    geom_point(data=geoMeanDf, aes(Sample_ID, log2(value + 1), col=variable)) +
    scale_color_manual(values=c(
      rev(brewer.pal(max(3, col_num/2), 'Greens'))[1:(col_num/2)],
      rev(brewer.pal(max(3, col_num/2), 'Reds'))[1:(col_num/2)]
    )) +
    labs(title=paste0("Probe Counts for ", experiment_name),
         y='log2(Count + 1)', color='GeoMetric Mean',
         subtitle=subtext) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5),
          legend.position = 'bottom',
          axis.text = element_text(color='black'))
  
  print(rawplot)
  
  dev.off()
  
  png(paste0(fname, '_boxplot_raw.png'), height=6, 
      width=max(length(sampset)/6, 11), units='in', res=300)
  
  boxdf <- bind_rows(lapply(unique(dedupdf$Sample_ID), function(sampid) {
    counts <- log2(dedupdf$Count[which(dedup_counts$Sample_ID == sampid)] + 1)
    return(data.frame(Sample_ID=as.character(sampid), 
             y0=min(counts), 
             y25=quantile(counts, 0.25),
             y50=median(counts),
             y75=quantile(counts, 0.75),
             y100=max(counts),
             stringsAsFactors = FALSE
             ))
  }))
  
  rawboxes <- 
  ggplot(boxdf, aes(x=Sample_ID)) + 
    geom_boxplot(width=0.8, aes(ymin=y0, lower=y25, middle=y50, upper=y75, ymax=y100),
                 stat='identity', size=0.8) +
    geom_point(data=geoMeanDf, aes(Sample_ID, log2(value + 1), col=variable),
               size=2) +
    scale_color_manual(values=c(
      rev(brewer.pal(max(3, col_num/2), 'Greens'))[1:(col_num/2)],
      rev(brewer.pal(max(3, col_num/2), 'Reds'))[1:(col_num/2)]
    )) +
    labs(title=paste0("Probe Counts for ", experiment_name),
         y='log2(Count + 1)', color='GeoMetric Mean',
         subtitle=subtext) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5),
          legend.position = 'bottom',
          axis.text = element_text(color='black'))
  
  print(rawboxes)
  
  dev.off()
  
  png(paste0(fname, '_distributions_raw.png'), height=6, 
      width=max(length(sampset)/6, 11), units='in', res=300)
  
  rawdistributions<- ggplot(subset(dedupdf, !Negative), 
         aes(x=log2(Count+1), color=Sample_ID, fill=Sample_ID)) +
    geom_density(alpha=0.5) +
    geom_density(data=subset(dedupdf, Negative), aes(y=-..density..),
                 alpha=0.5) +
    theme_bw()
  
  print(rawdistributions)
  
  dev.off()
  
}



if (!no.figures) {
  for (figset in names(figsets)) {
    rawProbeplot(dedup_counts, figsets[[figset]], 
                 paste0(fRoot, '_', figset))
  }
}

# Negative Data frame w/ GeoMean + LOQ by Sample + Pool
genNegFrame <- function(df) {
  negFrame <- dcast(subset(df, Negative, 
                           select = c('Sample_ID', 'Pool', 'Count')),
                    Sample_ID ~ Pool, fun.aggregate = ngeoMean, 
                    value.var = 'Count')
  negFrame <- melt(negFrame, id.vars = c('Sample_ID'))
  colnames(negFrame)[2:3] = c('Pool', 'NegGeoMean')
  tmp <- dcast(subset(df, Negative, select = c('Sample_ID', 'Pool', 'Count')),
               Sample_ID ~ Pool, fun.aggregate = ngeoSD, value.var = 'Count')
  tmp <- melt(tmp, id.vars = c('Sample_ID'))
  colnames(tmp)[2:3] = c('Pool', 'NegGeoSD')
  
  negFrame <- left_join(negFrame, tmp, by=c('Sample_ID', 'Pool'))
  negFrame$GeoLOQ2.5 <- negFrame$NegGeoMean * (negFrame$NegGeoSD^2.5)
  negFrame$GeoLOQ2.5 <- ifelse(negFrame$GeoLOQ2.5 < minLOQ, minLOQ, negFrame$GeoLOQ2.5)
  return(negFrame)
}

negFrame <- genNegFrame(dedup_counts)

collapsedCounts <- dcast(dedup_counts, Gene + Pool ~ Sample_ID, 
                         value.var = 'Count', fun.aggregate = ngeoMean)

for (pool in unique(negFrame$Pool)) {
  text_summary <- paste0(text_summary,
     'Counts for pool: ', pool,
     '\n\tMean of Collapsed Negatives: ',
     pnum(mean(negFrame$NegGeoMean[which(negFrame$Pool == pool)])),
     '\n\tMean of Collapsed Targets: ',
     pnum(mean(as.matrix(collapsedCounts[
       which(collapsedCounts$Pool == pool),
       3:ncol(collapsedCounts)
     ]))), '\n'
  )
}


if ( length(unique(collapsedCounts$Gene)) != nrow(collapsedCounts) ) {
  stop('Some genes are listed in multiple pools.')
}

# Write the Raw Collapsed Count Matrix
colnames(collapsedCounts)[1] <- 'TargetName'
write.table(collapsedCounts[,c(1,3:ncol(collapsedCounts))],
            paste0(fRoot, '_TargetCountMatrix.txt'),
            sep='\t', quote=F, row.names = F)
colnames(collapsedCounts)[1] <- 'Gene'

# Gene by Gene Plot
gene_plot <- function(count_df, fail_df, collapsedCounts, neg_df, gene) {
  gene_df <- rbind(
    subset(count_df, Gene == gene),
    subset(fail_df, Gene == gene)
  )
  gene_df$Flag <- paste0(gene_df$LowFlag, gene_df$Outliers)
  gene_df$Flag[which(gene_df$Flag == '')] <- 'Used'
  
  geneGeoMeans <- melt(subset(collapsedCounts, Gene == gene), 
                       id.vars=c('Gene', 'Pool'))
  
  pool = geneGeoMeans$Pool[1]
  
  rawplot <- ggplot(gene_df, aes(x=Sample_ID, y=log2(Count)) ) +
    geom_point(aes(fill=RNAID, shape=Flag), size=1.25) +
    geom_path(data=geneGeoMeans, aes(x=variable, y=log2(value)), 
              group='Exp_Mean', col='black', size=0.75) +
    geom_path(data=subset(neg_df, Pool == pool), 
              aes(x=Sample_ID, y=log2(GeoLOQ2.5)), group='Neg_LOQ', 
              col='red', size=0.75) +
    labs(title=gene) +
    scale_shape_manual(values=c('Used'=21, 'Global Outlier'=4, 
                                'Local Outlier'=8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5),
          legend.position = 'bottom',
          legend.title = element_blank(),
          axis.text = element_text(color='black'))
  
  return(rawplot)
}

gene_splitplots <- function(count_df, fail_df, collapsedCounts, neg_df, 
                            sampset, fname) {
  
  pdf(paste0(fname, '_counts_bytile.pdf'), width=max(length(sampset)/9, 11), 
      height=6)

  for (genename in unique(count_df$Gene)){
    print(gene_plot(count_df, fail_df, collapsedCounts, neg_df, genename))
  }
  
  dev.off()
}


if (!skipgpdf && !no.figures) {
  for (figset in names(figsets)) {
    gene_splitplots(dedup_counts, removed_probes, collapsedCounts, negFrame, 
                    figsets[[figset]], paste0(fRoot, '_', figset))
  }
}

# Generate Target properties data frame from lookup table
# Correlation to negatives from collapsed counts
# Global Outliers from removed_probes
if (!is.null(targetGroups_file)) {

  targetGroups <- read.delim(targetGroups_file, header = F,
                             stringsAsFactors = F)
  target_notes <- data.frame(TargetName=unique(rnaid_lookup_df$Gene),
                             HUGOSymbol=unique(rnaid_lookup_df$Gene))
  target_notes$TargetGroup <- sapply(target_notes$TargetName, 
      function(gene) {return(paste(unique(c('All Probes', 
           targetGroups$V2[which(targetGroups$V1 == gene)])),
           sep=';', collapse=';'))})
  target_notes$AnalyteType <- 'RNA'
  
  target_notes$Codeclass <- 'Endogenous'
  
  target_notes$Codeclass[
    which(target_notes$TargetName %in% housekeepers)] <- 'Control'
  target_notes$Codeclass[
    grep('NegProbe', target_notes$TargetName)] <- 'Negative'
  #rm(targetGroups)
} else if (length(pkc_vec) > 0) {

  target_notes$Codeclass[
    grep("Endogenous", target_notes$Codeclass)] <- "Endogeous"
  # Replace controls with user-designated housekeeper list
  if (length(housekeepers) > 0) {
    target_notes$Codeclass[
      grep("Control", target_notes$Codeclass)] <- "Endogeous"
    target_notes$Codeclass[which(target_notes$TargetName %in% housekeepers)] <- "Control"
  # Or use the pkc file designated housekeepers
  } else {
    target_notes$Codeclass[
      grep("Control", target_notes$Codeclass)] <- "Control"
  }
  target_notes$Codeclass[
      grep("Negative", target_notes$Codeclass)] <- "Negative"
}

# Reassign PoolNames and Numbers
poolNames <- unique(rnaid_lookup_df$Module)
poolNum <- formatC(seq(1, length(poolNames)), width=2, flag="0")
text_summary <- paste0(text_summary,
   'Pool mappings\n',
   paste0('\t', poolNames, ' = ', poolNum, collapse='\n'),
   '\n'
   )

if (!is.null(targetGroups_file) || length(pkc_vec) > 0) {
  target_notes$Pooling <- poolNum[match(
    rnaid_lookup_df$Module[
      match(target_notes$TargetName, rnaid_lookup_df$Gene)], 
    poolNames)]

  # Add names to the poolNum list
  short_pool <- poolNum
  names(short_pool) <- poolNames

  # Separate negative probe counts by pool
  negative_cors <- lapply(poolNames, 
                            function(curr_pool) { 
                              # Get negative counts for current pool
                              neg_count <- 
                                unlist(collapsedCounts[
                                         collapsedCounts$Pool == curr_pool & 
                                         collapsedCounts$Gene %in% 
                                         dedup_counts[dedup_counts$Negative, "Gene"],
                                         3:ncol(collapsedCounts)])
                              # Get list of all targets in pool
                              pool_targs <- target_notes[target_notes$Pooling == short_pool[curr_pool], "TargetName"]
                              targ_order <- which(collapsedCounts$Gene %in% pool_targs)
                              targ_names <- as.character(collapsedCounts$Gene[targ_order])
                              # Create data frame of current pool negative correlations
                              neg_cor <- cor(x=t(collapsedCounts[targ_order, 3:ncol(collapsedCounts)]), y=neg_count)
                              neg_cor_df <- data.frame(NegCor=neg_cor, row.names=targ_names)
                              return(neg_cor_df)
                            })
  pooled_cors <- do.call(rbind, negative_cors)
  # Append negative correlation values to target notes
  target_notes[["CorrelationToNegatives"]] <-
    pooled_cors[as.character(target_notes$TargetName), "NegCor"]

  # Append number of global outliers to target notes
  target_notes$GlobalOutliers <- sapply(target_notes$TargetName,
    function(gene) {
      
      globeCount <- length(unique(
        removed_probes$RNAID[intersect(which(removed_probes$Gene == gene),
            c(which(removed_probes$LowFlag != ''),
              which(removed_probes$Outliers == 'Global Outlier')))]
        ))
      
      totalCount <- length(unique(c(
        removed_probes$RNAID[which(removed_probes$Gene == gene)],
        dedup_counts$RNAID[which(dedup_counts$Gene == gene)]
        )))
      
      return(globeCount / totalCount)
    })
  target_notes <- target_notes[order(target_notes$TargetName), ]

  write.table(target_notes, file=paste0(fRoot, '_TargetProperties.txt'),
              sep='\t', quote=F, row.names=F)
  #rm(target_notes)
}

# Add the necessary columns to probe matrix
probe_matrix <- dcast(rbind(dedup_counts, removed_probes),
                      RNAID + Gene ~ Sample_ID, value.var = 'Count')
probe_matrix <- cbind(data.frame(ProbeName='', 
                           ProbeDisplayName=probe_matrix$RNAID,
                           TargetName=probe_matrix$Gene,
                           HUGOSymbol=probe_matrix$Gene, 
                           Accessions='',
                           GenomeBuild='', GenomicPosition='',
                           GlobalOutlier=FALSE,
                           GlobalOutlierReason='',
                           OutlierFrequency='',
                           OutlierAOIs='',
                           stringsAsFactors = FALSE),
                      probe_matrix[,2:ncol(probe_matrix)])

lowCountProbes <- as.character(unique(removed_probes$RNAID[
  which(removed_probes$LowFlag != '')]))
outlierProbes <- as.character(unique(removed_probes$RNAID[ 
    which(removed_probes$Outliers == 'Global Outlier')]))

probe_matrix$GlobalOutlier[
  which(probe_matrix$ProbeDisplayName %in% 
          c(lowCountProbes, outlierProbes))] <- TRUE

probe_matrix$GlobalOutlierReason[
  which(probe_matrix$ProbeDisplayName %in% lowCountProbes)] <- 'Low Count'
probe_matrix$GlobalOutlierReason[
  which(probe_matrix$ProbeDisplayName %in% outlierProbes)] <- 'Grubbs Outlier'
probe_matrix$OutlierAOIs <- sapply(probe_matrix$ProbeDisplayName, 
    function(probe) {
      localrows <- intersect(
        which(removed_probes$RNAID == probe),
        which(removed_probes$Outliers == 'Local Outlier')
        )
      return(paste(
        removed_probes$Sample_ID[localrows], collapse=';'
      ))
    })

# Create column with list of AOIs where Probe was a local outlier
write.table(probe_matrix, paste0(fRoot, '_BioProbeCountMatrix.txt'),
            sep='\t', quote=F, row.names=F)
#rm(probe_matrix, removed_probes, dedup_counts)

# Make the Overall count stripplot with target geomeans instead of probe counts
collapsedCounts <- melt(collapsedCounts, id.vars = c('Gene', 'Pool'))
colnames(collapsedCounts) <- c('Gene', 'Pool', 'Sample_ID', 'Count')

rawTargetplot <- function(countdf, deduped_df, sampset, fname) {
  countdf <- countdf[countdf$Sample_ID %in% sampset,]
  countdf$Negative <- countdf$Gene %in% deduped_df[deduped_df$Negative, "Gene"]
  subtext <- paste0(
    "Range of Collapsed Counts: ", pnum(min(countdf$Count)), " - ",
    pnum(max(countdf$Count)),
    "\nGeoMean Collapsed Counts for Targets: ", 
    pnum(ngeoMean(countdf$Count[!countdf$Negative])),
    "\nGeoMean Collapsed Counts for Negatives: ",
    pnum(ngeoMean(countdf$Count[countdf$Negative]))
  )
  
  # Create Geometric Mean Data Frame
  targetGeoMeans <- dcast(subset(countdf, !Negative, 
                                 select=c(Sample_ID, Pool, Count)),
                          Sample_ID ~ Pool, value.var='Count', 
                          fun.aggregate = ngeoMean)
  negativeGeoMeans <- dcast(subset(countdf, Negative, 
                                   select=c(Sample_ID, Pool, Count)),
                            Sample_ID ~ Pool, value.var='Count', 
                            fun.aggregate = ngeoMean)
  colnames(negativeGeoMeans)[2:ncol(negativeGeoMeans)] <- (
    paste0('Neg-', colnames(negativeGeoMeans)[2:ncol(negativeGeoMeans)])
  )
  
  geoMeanDf <- rbind(melt(targetGeoMeans), melt(negativeGeoMeans))
  
  png(paste0(fname, '_collapsedcounts_raw.png'), height=6, 
      width=max(length(sampset)/6, 11), units='in', res=300)
  
  col_num <- length(unique(geoMeanDf$variable))
  
  rawplot <- ggplot(countdf, aes(x=Sample_ID, y=log2(Count + 1))) + 
    geom_jitter(width=0.1, height=0) +
    geom_point(data=geoMeanDf, aes(Sample_ID, log2(value + 1), col=variable)) +
    scale_color_manual(values=c(
      rev(brewer.pal(max(3, col_num/2), 'Greens'))[1:(col_num/2)],
      rev(brewer.pal(max(3, col_num/2), 'Reds'))[1:(col_num/2)]
    )) +
    labs(title=paste0("Collapsed Counts for ", experiment_name),
         y='log2(Count + 1)', color='GeoMetric Mean',
         subtitle=subtext) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5),
          legend.position = 'bottom',
          axis.text = element_text(color='black'))
  
  print(rawplot)
  
  dev.off()
  
}

if (!no.figures) {
  for (figset in names(figsets)) {
    rawTargetplot(collapsedCounts, dedup_counts, figsets[[figset]], 
                 paste0(fRoot, '_', figset))
  }
}

# Report Gene above LOQ (2.5) summary stats
collapsedCounts$GeoLOQ2.5 <- negFrame$GeoLOQ2.5[match(
  paste0(collapsedCounts$Sample_ID, collapsedCounts$Pool),
  paste0(negFrame$Sample_ID, negFrame$Pool))]
collapsedCounts$LOQCheck <- collapsedCounts$GeoLOQ2.5 < collapsedCounts$Count
processing_sum$OverLOQ <- sapply(unique(as.character(
  processing_sum$Sample_ID)), function(samp) {
    return(length(intersect(
      which(collapsedCounts$Sample_ID == samp),
      which(collapsedCounts$LOQCheck)
      )))
  })
text_summary <- paste0(text_summary, 
   length(which(processing_sum$OverLOQ > 0)), 
   ' AOIs out of ', nrow(processing_sum),
   ' had at least one target with expression greater than',
   ' an LOQ of GeoMean(Negatives) * GeoSD(Negatives)^2.5\n',
   'Min Significant: ', min(processing_sum$OverLOQ), '\n',
   'Median Significant: ', median(processing_sum$OverLOQ), '\n',
   'Max Significant: ', max(processing_sum$OverLOQ), '\n'
   )


# Calculate Normalization Factors
nanoNormFactors <- function(v) {
  return(v / ngeoMean(v))
}

negFrame$NormFactorNeg <- 0
for (pool in unique(negFrame$Pool)) {
  
  poolrows <- which(negFrame$Pool == pool)
  negFrame$NormFactorNeg[poolrows] <- 
    nanoNormFactors(negFrame$NegGeoMean[poolrows])
  
}

for (pool in unique(negFrame$Pool)) {
  # Pull Pool rows
  subframe <- negFrame[which(negFrame$Pool == pool),]
  # Put Rows in same order as annotations frame and grab relevant columns
  subframe <- subframe[match(annotations$Sample_ID, subframe$Sample_ID),
                       c('NegGeoMean', 'NegGeoSD', 'GeoLOQ2.5', 'NormFactorNeg')]
  # Rename columnes
  newpool <- poolNum[which(poolNames == as.character(pool))]
  colnames(subframe) <- paste0(colnames(subframe), '_', newpool)
  # Add to annotations
  annotations <- cbind(annotations, subframe)
  rm(subframe)
}

collapsedCounts <- dcast(collapsedCounts, Gene + Pool ~ Sample_ID, 
                         value.var = 'Count')

collapsedCounts$Pool <- poolNum[match(as.character(collapsedCounts$Pool), poolNames)]

q3s <- nanoNormFactors(apply(collapsedCounts[3:ncol(collapsedCounts)], 2, 
                             quantile, probs = 0.75,  na.rm=T))

annotations$NormFactorQ3 <- q3s[match(annotations$Sample_ID, names(q3s))]

hks <- nanoNormFactors(
  apply(collapsedCounts[which(collapsedCounts$Gene %in% housekeepers),
                        3:ncol(collapsedCounts)], 2, ngeoMean)
  )

annotations$NormFactorHK <- hks[match(annotations$Sample_ID, names(hks))]

# Write Negative Normalized Sheet
negNorm <- function(count_df, segmentProperties) {
  countcols <- 3:ncol(count_df)
  return(cbind(count_df[,1:2], bind_rows(lapply(1:nrow(count_df), function(rownum) {
    pool <- count_df$Pool[rownum]
    return( 
      count_df[rownum, countcols] /
      segmentProperties[match(colnames(count_df[countcols]), 
                              segmentProperties$Sample_ID), 
                        paste0('NormFactorNeg_', pool)]
    )
  }))))
}

negNormTable <- negNorm(collapsedCounts, annotations)

negNormTable <- negNormTable[,c(1, 3:ncol(negNormTable))]

colnames(negNormTable)[1] <- 'TargetName'

write.table(negNormTable, 
            paste0(fRoot, '_NegNorm_TargetCountMatrix.txt'),
            sep='\t', quote=FALSE, row.names=FALSE)

# Write SegmentProperties
write.table(annotations, 
            paste0(fRoot, '_SegmentProperties.txt'),
            sep='\t', quote=FALSE, row.names=FALSE)

cat(text_summary, file=paste0(fRoot, '_DatasetHistory.txt'))
