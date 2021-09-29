#' @title DSP differential expression analysis (parallel)
#'
#' @description Wrapper function to run differential expression in parallel. This function
#' is used prior to \code{\link{dsp_de_analysis}}.
#'
#' @param \code{de_results} object that can be one of the following:
#' \enumerate{
#'    \item{NULL. If not provided, the differential expression analysis will run}
#'    \item{object of type \code{character} specifying the location of the previous saved DE results}
#'    \item{object of type \code{data.frame} providing the loaded DE results}
#' }
#'
#' @param \code{norm_counts} is \code{data.frame} containing the normalized count data
#' to be analyzed. Each row corresponds to a gene of interest. The first column is assumed
#' to be labeled 'Gene'. Each additional column corresponds to a Sample_ID. Sample_ID values
#' have a matching row in the annotations objects (see \code{samp_notes} argument below).
#'
#' @param \code{samp_notes} is a \code{data.frame} that provides annotations. If DE is to be
#' run (i.e., \code{!is.null(de_results)}), \code{samp_notes} must have
#' the grouping_var column, with containing base_level factor, and it must have a control_var
#' column.
#'
#' @param \code{grouping_var} is the column in the annotations data.frame (i.e., \code{samp_notes})
#' to use for grouping. This is the level 2 factor in the Mixed Effects model for which estimates
#' are made. This can be a global object generated from \code{read_config()} if running
#' the Hydra pipeline.
#'
#' @param \code{base_level} is a level in the \code{grouping_var} column used as base. This can be
#' a global object generated from \code{read_config()} if running the Hydra pipeline.
#'
#' @param \code{control_var} is the column name in the annotations data.frame (i.e., \code{samp_notes})
#' used at the level 1 random effect. For DSP, this is usually an individual ID for which several
#' AOIs are nest within. If it is invariable (e.g., a vector of 1s), the default behavor is to run
#' a linear model using R's base::lm function. Otherwise, a Mixed-effects model will be used via
#' lme4::lmer.
#'
#' @param \code{de_dir} the parent directory to store the de_results
#'
#' @param \code{n_processors}. The number of processors to use. Default is 4.
#'
#' @param \code{the_formula}. Either the classic formula (default) or a custom formula. Do to customization possibilities,
#' there is no guarantee that this will run, converge, or make sense. Example of custom formula =
#' 'a_gene ~ Disease * Tissue + (1 + Disease | DSP_scan)'. Note that 'a_gene' should be left as is.
#'
#' @param \code{show_pb}. Logical to show the progress bar (TRUE) or not (FALSE). Default is TRUE.
#'
#' @param \code{min_sample_size}. An integer of 1 or more than specifies the minimum number of samples
#' needed in a group for that group to be used in DE. This is similar to the IO360 report parameter
#' \code{minSampleSize}. Default is 3.
#'
#' @details The function strips down the main \code{\link{dsp_de_analysis}} function and performs the
#' differential expression analysis over a specified number of processors. It uses parallel::parLapply to process the
#' normalized expression for each gene across different groups of interest (i.e., as specified by the \code{grouping_var}
#' object). By default, this function runs the same as \code{\link{dsp_de_analysis}} in terms of the formulas used.
#' If the \code{control_var} column is invariant. A simple linaer model will be used. Otherwise, a Mixed effects model
#' will be used of the form 'a_gene ~ grouping_var + (1 + grouping_var | control_var)'. In both of these scenarios,
#' a data.frame will be saved to disk and returned. If, however,  \code{the_formula} differs from
#' its default value of "classic", a custom formula will be used. In such a case, a list of list where the final list
#' provides the 1) gene, the 2) test (e.g., 'normal vs disease'), and the model. There is no error checking for convergence
#' etc. if the_formual != "classic". In all Mixed effect models, REML=TRUE.
#'
#' @author Tyler Hether
#'
#' @reference None.
#'
#' @return Objects and files returned depend on the input:
#' \enumerate{
#'    \item{If \code{the_formula=="classic"}, an object of class \code{data.frame} is returned. The
#'    dataframe is also sent to disk.}
#'    \item{Otherwise, a list of list is return (see details)}
#' }
#'
#' @seealso \code{\link{dsp_de_analysis}} for plotting.
#'
#' @export dsp_de_analysis_parallel
#'
#' @examples
#'
#'# Simulate 2 genes in log2 space. Look at the true fold
#'# change (for gene_a) and compare to the estimated fold change.
#'n_samples <- 50
#'n_aois <- 10
#'constant <- 4.5
#'# Simulate normalized, log2
#'make_aois <- function(mu, sd, n_aois, group){
#'  out <- data.frame(Group=rep(group, n_aois),
#'      gene_a=rnorm(n=n_aois, mean=mu, sd=sd))
#'      return(out)
#'}
#'# Different effect size of group A from mean
#'diff_a <- rnorm(n=n_samples, mean=2, sd=1)
#'# Make samples.
#'A <- do.call(rbind, lapply(1:n_samples, function(i){
#'  samp <- make_aois(constant+diff_a[i], 1, n_aois, "A")
#'  samp$ind <- i
#'  return(samp)
#'}))
#'diff_b <- rnorm(n=n_samples, mean=-2, sd=1)
#'B <- do.call(rbind, lapply(1:n_samples, function(i){
#'  samp <- make_aois(constant+diff_b[i], 1, n_aois, "B")
#'  samp$ind <- i+max(A$ind)
#'  return(samp)
#'}))
#'df <- rbind(A, B)
#'require(plyr)
#'require(dplyr)
#'df <- ddply(df, .(Group, ind), function(x){
#'  x$aoi <- 1:nrow(x)
#'  return(x)
#'})
#'df <- df %>% mutate(Sample_ID=paste(Group, ind, aoi, sep="_"), gene_b=gene_a+rnorm(1))
#'# Make sure there's no negative values.
#'df$gene_a[which(df$gene_a<=0)] <- 1e-03
#'df$gene_b[which(df$gene_b<=0)] <- 1e-03
#'library(tidyr)
#'norms <- rbind(spread(df %>% dplyr::select(-aoi, -gene_b, -Group, -ind), Sample_ID, gene_a),
#'  spread(df %>% dplyr::select(-aoi, -gene_a, -Group, -ind), Sample_ID, gene_b))
#'# "Convert" to linear scale since the function will log2 transform internally
#'norms <- 2^norms
#'# Combine
#'norms <- cbind(data.frame(Gene=c("gene_a", "gene_b")), norms)
#'row.names(norms) <- norms$Gene
#'# Convert the df to respective dataframes
#'annotations <- df %>% select(Sample_ID, Group, ind)
#'library(ggplot2)
#'ggplot(df, aes(x=factor(ind), y=gene_a)) +
#'geom_jitter(height=0) + ylab("log2 expression for simulated gene_a")
#'# Means of the groups in log2 space
#'means <- ddply(df, .(Group), summarize, mean(gene_a))[,2]
#'# Fold change where vec is the means in log2 space
#'# returns FC in linear space
#'calc_fc <- function(vec){
#'  vec2 <- 2^vec
#'  return(vec2[2] / vec2[1])
#'}
#'# This is the simulated fold change in log2 space
#'print(log2(calc_fc(means)))
#'#Compare that fold change to the estimate for gene_a
#'library(dspNgs)
#'dsp_de_analysis_parallel(de_results=NULL, norm_counts=norms,
#'samp_notes=annotations, grouping_var="Group", base_level="A",
#'  control_var="ind", de_dir="./", n_processors=4, the_formula="classic")
#'# End run.

dsp_de_analysis_parallel <- function(de_results=NULL, norm_counts, samp_notes, grouping_var, base_level,
                                control_var, de_dir, n_processors=4, the_formula="classic", show_pb=TRUE, min_sample_size=3){
  ### ###############
  ### For development
  ### ###############
  # Toggle these ON/OFF for development.
  # library(dspNgs)
  # message("Waring! Dev testing is *ON*.")
  # load("./inst/testData/cta_test/cornell_collab.rds")
  # de_results <- NULL
  # norm_counts <- dfs[[2]]
  # samp_notes <- dfs[[3]]
  # grouping_var <- "Disease"
  # base_level <- "Covid-19"
  # control_var <- "DSP_scan"
  # de_dir <- "./"
  # n_processors <- 6
  # the_formula="classic" # A custom example: the_formula='a_gene ~ Disease + (1 + Disease | dsp_collection_plate/DSP_scan)'
  # show_pb = TRUE
  # min_sample_size=3
  #~~~~~~~~~~~~~~~~~~

  # Logging the execution of this run with provided parameters
  send_to_log(grouping_var, base_level, control_var)

  # Check the input.
  check_input(de_results, norm_counts, samp_notes, grouping_var, base_level, control_var, show_pb, min_sample_size)

  # Run the DE if needed.
  de_res <- de_logic(de_results, norm_counts, samp_notes, grouping_var, base_level, control_var,
              n_processors, the_formula, show_pb, min_sample_size)

  # Write results to disk
  if(class(de_res)=="data.frame"){
    # Write results to disk for classic formula cases.
    # Note, some WTA gene names have a comma. To generalize, change "," to "." for all
    # cases so that parsing/reading the ensuing csv is easier. This will not change
    # the original results.
    de_res_to_write <- de_res
    de_res_to_write$gene <- gsub(pattern=",", replacement=".", x=de_res_to_write$gene)
    write.csv(de_res_to_write, file=paste(de_dir, "de_results.csv", sep="/"))
  } # else {
    # # Note. For larger datasets, this can take prohibitively long to save so
    # # this feature is not turned on.
    # # return the list of lists to disk. (Tried two methods on a CTA dataset with ~300 AOIs)
    # save(x=de_res, file=paste(de_dir, "de_results_list.rds", sep="/"))
    # rlist::list.save(x=de_res, file=paste(de_dir, "de_results_list.rds", sep="/"))
    # }

  # Return the de results.
  return(de_res)

}


#' @title Internal lower-level function that actually computes the DE.
#' returns the data.frame of results.
#' This function does the following:
#' 1. Transposes and transforms the entire norm_counts dataset.
#' 2. Create all pairwise combinations that have at least min_sample_size observations per group.
#' 3. lapply over each element of #2.
#' 4. Process over each gene in #3.
#' 5. Collate data.frames (or lists) and return.
de_internal <- function(norm_counts, samp_notes, grouping_var, base_level, control_var, n_processors, the_formula, show_pb, min_sample_size){
  # See main wrapper function for object attributes.

  # 1. log2 normalize counts and transpose dataframe
  rownames(norm_counts) <- norm_counts$Gene
  gene_log2 <- as.data.frame(log2(t(norm_counts %>% select(-Gene))))
  the_genes <- colnames(gene_log2) # Used later
  gene_log2$Sample_ID <- rownames(gene_log2)

  # Merge log2 data with annotations. Note that this explicitly keeps all columns of samp_notes in case a custom
  # formula is desired at some ponit.
  if(!("Sample_ID" %in% colnames(samp_notes))){
    # This ensures that Sample_ID is present as a column in dfs[[3]] prior to merging.
    stop("A column named Sample_ID must be included in the sample annotations (i.e., dfs[[3]]).\n\tPlease add Sample_ID and re-run.")
  }
  df_merged <- merge(gene_log2, samp_notes, by="Sample_ID")

  # 2. Create a list of combinations to run.
  df_combos <- create_combos(vec=as.character(df_merged[[grouping_var]]), base=base_level, the_min=min_sample_size)

  # 3. Loop through each element of df_list
  res <- lapply(df_combos, function(a_combo){
    # For testing.
    # a_combo <- c("Covid-19", "Flu")

    # Get index for message.
    index <- which((as.data.frame(do.call(rbind,df_combos)) %>% mutate(combo = paste0(V1, "_vs_", V2)))$combo == paste0(a_combo[1], "_vs_", a_combo[2]))

    # Subset the data to only contain the combination of interest
    df_x <- df_merged[df_merged[[grouping_var]] %in% a_combo, ]

    # Make the level2 factor: 0 == first (base_line) and 1 is the other group
    # Note: this won't be used if the_formula != "classic"
    df_x$level2 <- factor(ifelse(df_x[[grouping_var]]==a_combo[1], 0, 1))

    # STDOUT msg
    msg <- paste0("Working on ", a_combo[1], " vs ", a_combo[2], " (comparison ", index, " of ", length(df_combos), ")")
    message(msg)
    flog.info(msg, name="DSP_NGS_log")

    # 4. Process each gene in parallel
    cl <- makeCluster(n_processors)
    clusterExport(cl=cl, varlist=c( "df_x", "the_genes", "a_combo", "the_formula", "control_var", "LMER"), envir=environment())
    if(show_pb){
      require(pbapply)
      inner_res <- pbapply::pblapply(the_genes, LMER, cl=cl)
    } else{
      inner_res <- parallel::parLapply(cl, the_genes, LMER)
    }

    stopCluster(cl)

    # If the_formula is 'classic', we know the structure of the data and we can rbind them together
    # and add additional metrics.
    if(the_formula=="classic"){
      inner_res <- do.call(rbind, inner_res) %>% arrange(gene)
      inner_res$test <- paste0(a_combo[1], " vs ", a_combo[2])
      # adjust p-values
      inner_res$fdr <- p.adjust(inner_res$Pval, 'fdr')
      inner_res$fwer <- p.adjust(inner_res$Pval, 'bonferroni')

    }

    # Return the DE results.
    return(inner_res)

  })

  # 5. Depending on the output structure, collate into a single data frame or a list (non-classic formula)
  if(the_formula=="classic"){
    to_return <- do.call(rbind, res)
  } else {
    # Not classic so res is a list of list of models
    to_return <- res
  }

  return(to_return)
}

# Inner function to create combinations:
# For n groups, there will be up to n choose 2 combinations.
# A combination will not appear if one or more of its groups has fewer than the_min samples.
# the focal baseline group gets priority and others are filled in based on the default order.
# Exammple: 3 groups (A, B, C) with B as baseline. Each group as >=the_min observations:
#           The result will be a list with elements c(B, A), c(B, C), and c(A, C)
create_combos <- function(vec, base, the_min){
  # Testing.
  # vec=as.character(df_merged[[grouping_var]]) # a character vector with groupings
  # base=base_level # The focal level (base_level) to use. This will be prioritize as the first position
  # the_min = 3 # This is the minmum number of samples that each group in the given pairwise comparison needs to be considered

  # There must be at least the_min observations for a group to be considered
  tab <- as.data.frame(table(vec))
  levs <- as.character(tab$vec) # This is a character vector that might get trimmed downstream
  tab_too_small <- dplyr::filter(tab, Freq<the_min)
  rm_base <- FALSE
  if(nrow(tab_too_small)>0){
    # There's at least one group that needs to be culled.

    if(nrow(tab) <= (nrow(tab_too_small)+1)){
      # This condition is TRUE if and only if there
      # are 1 or fewer groups with at least the_min observations.
      # In that case, DE should not run. Issue a stop.
      msg <- paste0("There were ", nrow(tab), " groups detected but ", nrow(tab_too_small),
      " of these groups had fewer than ", the_min, " observations. Check to ensure the grouping variable is formatted properly with multiple observations per group and at least 2 groups with >=", the_min, " observations. Aborting.")
      flog.info(msg, name="DSP_NGS_log")
      stop(msg)
    }

    msg <- paste0("The following groups had fewer than ", the_min, " observations and will be culled: ", paste0(tab_too_small$vec, collapse=", "))
    flog.info(msg, name="DSP_NGS_log")
    warning(msg)

    if(base %in% tab_too_small$vec){
      # The base_level itself needs to be culled due to too few observations.
      # Specifically call this situation out since the user is expecting this as baseline.
      rm_base <- TRUE
      msg <- paste0("The base level ", base, " has fewer than ", the_min, " observations and will not be used in DE but other combinations will be used. Please check to ensure the grouping variable is formatted properly.")
      flog.info(msg, name="DSP_NGS_log")
      warning(msg)
    }

    # Actually setting the levs
    levs <- setdiff(tab$vec, tab_too_small$vec)
  } else {
    msg <- paste0("All groups had a sufficient number of observations (>=", the_min, ") for DE analysis")
    flog.info(msg, name="DSP_NGS_log")
  }

  # Put base in front (if it wasn't removed)
  if(!rm_base){
    levs <- c(base, setdiff(levs, base))
  }

  # Return combinations
  to_return <- combn(levs, 2, simplify=FALSE)

  return(to_return)

}

# Inner portion of de_internal to call, agnositic with respect to show_pb
LMER <- function(a_gene){
  require(lme4)
  require(lmerTest)
  require(dplyr)
  df_x$a_gene <- df_x[,which(colnames(df_x)==a_gene)] # renaming avoids gene names with special characters (e.g., "-") in the formula.
  if(the_formula=="classic"){
    # the 'classic', non-custom formula is desired:
    if(length(unique(df_x[[control_var]]))==1){
      # All have the same control so this is just a linear model.
      form <- paste0('a_gene ~ level2')
      mod <- lm(as.formula(form), data=df_x)
      cf <- data.frame(coefficients(summary(mod)))[2,]
      colnames(cf)[1] <- "FC"
      colnames(cf)[4] <- "Pval"
      cf$df <- summary(mod)$df[2]
      cf <- cf[,c(1,2,5,3,4)]
      cf$gene <- a_gene
      cf$Significance <- -log10(cf$Pval)
    } else {
      # This has a multi-level structure:
      form <- paste0('a_gene ~ level2 + (1+level2|', control_var, ")")
      mod <- lmerTest::lmer(as.formula(form), data=df_x, REML=TRUE,
                            # hard prevent singular boundary warnings:
                            control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
      # For the classic format, here is a summary of the results:
      cf <- data.frame(coefficients(summary(mod)))[2,]
      cf$gene <- a_gene
      colnames(cf)[1] <- "FC"
      colnames(cf)[5] <- "Pval"
      cf$Significance <- -log10(cf$Pval)
    }
    # reuturn the cf
    row.names(cf) <- NULL
    return(cf)
  } else {
    # A custom formula is desired.
    # the formla = "a_gene ~ Disease * Tissue + (1 + Disease | dsp_collection_plate/DSP_scan)"
    # the_formula = "a_gene ~ level2 + (1 | dsp_collection_plate/DSP_scan)" # Used for testing.
    mod <- lmerTest::lmer(as.formula(the_formula), data=df_x, REML=TRUE,
                          # hard prevent singular boundary warnings:
                          control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    # Since this could multiple things. Return the *model* itself and let the user decide downstream analysis.
    the_test <- paste(a_combo, collapse = " vs ")

    return(list(test=the_test, gene=a_gene, model=mod))
  }
}

#' @title Internal lower-level function that checks de_ressults attributes and runs de_internal
#' if needed. Returns the de_results (whether it was computed, referenced by file, or referenced
#' by object).
de_logic <- function(de_results, norm_counts, samp_notes, grouping_var, base_level, control_var,
                         n_processors, the_formula, show_pb, min_sample_size){
  # See parent function for argument details.

  # There's a cascade of logic depending on what de_results is ccomprised of.
  if(is.null(de_results)) {
    # No de_results were provided so run the DE analysis via the de_internal function
    de_res <- de_internal(norm_counts, samp_notes, grouping_var, base_level, control_var, n_processors, the_formula, show_pb, min_sample_size)
  } else if(is.character(de_results)) {
    # Make sure this file is present and a valid filepath.
    if(length(Sys.glob(de_results))==1){
      # Read in the file an warn user if only 1 column was found since it might affect downstream analysis.
      message("Reading in csv file")
      de_res <- read.delim(de_results, sep = ',', header = TRUE, as.is = TRUE)
      if(ncol(de_res)==1){
        warning(paste0("de_results file only had 1 column. Was it in csv format?"))
      }
    } else {
      msg <- paste0("The de_results file, ", de_results, ", was provide but not found. Check path. Aborting.")
      flog.info(msg, name="DSP_NGS_log")
      stop(msg)
    }
  } else if(is.data.frame(de_results) & all(c('FC','Pval','gene') %in% colnames(de_results))) {
    de_res <- de_results
  } else {
    # If the logic made it to this point, de_results was an object but neither a character nor the expected data.frame. Abort.
    msg <- paste0("The de_results file, ", de_results, ", was a non-character object but not a data.frame with FC, Pval, and gene. Check object. Aborting.")
    flog.info(msg, name="DSP_NGS_log")
    stop(msg)
  }

  # return the de_res
  return(de_res)
}

#' @title Internal helper function for writing to log.
# The arguments here feed directly from the main function expect for
# the following, which are pulled from Global: gene_group
send_to_log <- function(grouping_var, base_level, control_var){
  flog.info("\n \n ############   Experimental Parameters   ################## \n", name="DSP_NGS_log")
  flog.info(paste0("Gene grouping variable: ", grouping_var), name="DSP_NGS_log")
  flog.info(paste0("Base Level (focal): ", base_level), name="DSP_NGS_log")
  flog.info(paste0("Control: ", control_var), name="DSP_NGS_log")
}

#' @title Internal helper function for checking input paratmeters
check_input <- function(de_results, norm_counts, samp_notes, grouping_var, base_level, control_var, show_pb, min_sample_size){
  # Check to ensure that grouping_var, base_level, and control_var are elements
  # samp_notes if DE needs to be performed.
  if(is.null(de_results)){
    # Checking grouping_var and base_level within
    if(!(grouping_var %in% colnames(samp_notes))){
      msg <- paste0("The grouping_var, ", grouping_var, ", is not found in samp_notes. Aborting.")
      flog.info(msg, name="DSP_NGS_log")
      stop(msg)
    } else if(!(base_level %in% samp_notes[[grouping_var]])){
      msg <- paste0("The grouping_var, ", grouping_var, ", was found but ", base_level, " was not found within. Aborting.")
      flog.info(msg, name="DSP_NGS_log")
      stop(msg)
    }
    # Checking that control_var is present as a column in samp_notes
    if(!(control_var %in% colnames(samp_notes))){
      msg <- paste0("The control_var, ", control_var, ", is not present in samp_notes. Aborting.")
      flog.info(msg, name="DSP_NGS_log")
      stop(msg)
    }
  }
  # Ensure that "Gene" is the first column in norm_counts since it will be
  # referenced by index downstream.
  if(colnames(norm_counts)[1] != "Gene"){
    msg <- "Column Gene is expected to be the first column in norm_counts but it is not. Check format. Aborting."
    flog.info(msg, name="DSP_NGS_log")
    stop(msg)
  }
  # Ensure that show_pb is logical (TRUE/FALSE) or equal to 0 or 1
  if(!is.logical(show_pb)){
    if(is.numeric(show_pb)){
      if(!(show_pb %in% c(0,1))){
        stop(paste0("parameter show_pb, ", show_pb, " ,is numeric but not 0/1. Please check format. Aborting."))
      }
    } else {
      stop(paste0("parameter show_pb, ", show_pb, " ,is neither logical (TRUE/FALSE) or 0/1. Please check format. Aborting."))
    }
  }
  # Ensure that min_sample_size is an integer of at least 1.
  if(!is.numeric(min_sample_size) & !is.integer(min_sample_size)){
    stop(paste0("parameter min_sample_size, ", min_sample_size, " ,is not an integer >=1. Aborting."))
  } else {
    # Make sure it's >=1
    if(is.numeric(min_sample_size) & (min_sample_size != round(min_sample_size))){
      stop(paste0("parameter min_sample_size, ", min_sample_size, " ,is numeric but not integer >=1. Aborting."))
    }
    if(min_sample_size<1){
      stop(paste0("parameter min_sample_size, ", min_sample_size, " ,is numeric/integer but < 1. Aborting."))
    }
  }
  flog.info("Finished checking input", name="DSP_NGS_log")
}

