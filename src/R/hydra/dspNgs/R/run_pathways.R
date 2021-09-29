#' @title run_pathways
#'
#' Pathway Analysis giving either coverage or GSEA
#'
#' @param df list of data frames containing input data
#' @param de_results the differential expression data
#' @param genListrank the way to rank the gene list either by foldchange, pval or FDR
#' @param path_fc the fold change cutoff to use for gene list
#' @param path_pval the pvalue to threshold the enrichment results by
#' @param enrichment the type of enrichment to be returned
#' @param outdir the output directory for plots
#' @param custom_pathways_path the directory where custom gmt files are located
#' @param exclude_reactome to exclude reactome pathways from custom pathway analysis
#'
#' @return enrich_res The results from your chosen enrichment
#'
#' @examples
#'  run_pathways(df = dfs, de_data = de_results, geneListrank=geneListrank, path_fc=path_fc, path_pval=path_pval, enrichment=enrichment, outdir = path_dir, custom_pathways_path = NULL, exclude_reactome = FALSE)
#'
#' @export custom_error
#' @export run_pathways
#' @export load_gmts
#' @export symbol_to_entrez
#' @export filter_custom_pathways

run_pathways <- function(df,
                         de_results,
                         geneListrank,
                         path_fc,
                         path_pval,
                         enrichment,
                         outdir,
                         custom_pathways_path = NULL,
                         exclude_reactome = FALSE,
                         ...) {

  # Convert to entrez id and set as rownames
  norm_data <- df[[2]]
  norm_data <- symbol_to_entrez(path_df=norm_data, gene_col="Gene", species = species)
  rownames(norm_data) <- norm_data$entrez
  norm_data <- as.matrix(norm_data[, -c(1, dim(norm_data)[2])])

  # Convert de results to entrez id
  de <- de_results
  de <- symbol_to_entrez(path_df=de, gene_col="gene", species = species)
  gene_ids <- as.character(de[,"entrez"])

  # Make lists for your results
  ora_list <- list()
  fgsea_list <- list()
  # Saved for future feature: DE on ssGSEA
  #ssGSEA_list <- list()

  # Pull reactome pathway for default analysis
  while(!is.null(custom_pathways_path)) {
    # Check if path exists otherwise run default reactome
    if(!file.exists(custom_pathways_path)) {
      custom_pathways_path <- custom_error(custom_path=custom_pathways_path,
                                              msg="does not exist")
      break
    }

    # Generate list of all custom gmt files
    custom_pathways_files <- list.files(path=custom_pathways_path,
                                          pattern="\\.gmt$",
                                          full.names=TRUE)

    # Check if list is not empty otherwise run default reactome
    if(length(custom_pathways_files) < 1) {
      custom_pathways_path <- custom_error(custom_path=custom_pathways_path,
                                              msg="has no gmt files")
      break
    }

    # Remove empty gmt files
    gmt_info <- file.info(custom_pathways_files)
    custom_pathways_files <- rownames(gmt_info[!(gmt_info$size == 0), ])

    # Check if custom files list is empty after empty file removal
    if(length(custom_pathways_files) < 1) {
      custom_pathways_path <- custom_error(custom_path=custom_pathways_path,
                                              msg="only has empty gmt files")
      break
    }

    # Load and append custom pathway lists based on genes analyzed
    custom_pathways_list <- lapply(custom_pathways_files,
                                      load_gmts,
                                      list_of_genes=gene_ids)
    pathways <- do.call(c, custom_pathways_list)

    # Append reactome pathways if not excluding
    if(!exclude_reactome) {
      # Set up your Reactome Pathways based on genes analyzed
      react_pathways <- reactomePathways(gene_ids)
      pathways <- c(react_pathways, pathways)
    }

    # To break out of loop
    break
  }

  # Use default reactome pathways analysis
  if(is.null(custom_pathways_path)){
    # Set up your Reactome Pathways, pulling from your ranked list from above
    pathways <- reactomePathways(gene_ids)
  }

  wrap_pathways <- pathways

  names(wrap_pathways) <- str_wrap(names(wrap_pathways), 25)

  # For loop to handle multiple tests
  for(test in unique(de$test)){

    test_dir <- paste(outdir, test, sep = "/")

    dir.create(test_dir)

    # Subset for each test
    sub_de <- de[de$test == test, ]

    # Create a ranked list according to fold change for the DE results
    geneList <- sub_de[,geneListrank]
    names(geneList) <- as.character(sub_de[,"entrez"])
    geneList <- sort(geneList, decreasing = FALSE)

    # Pull the genes that are below pval cutoff
    if(!any(abs(geneList) > path_fc)) {
      stop('No genes significant at the set cutoff, please adjust the cutoff and re-run read_config()\n')
    }

    # Run ORA
    print(paste0("Performing Over Representation Analysis (ORA) for test: ", test))
    ora_res <- run_ora(geneList = geneList, 
                       pathways = wrap_pathways, 
                       test = test, 
                       path_fc = path_fc, 
                       path_pval = path_pval, 
                       outdir = test_dir)
    ora_list <- append(ora_list, ora_res)

    # Run GSEA
    print(paste0("Performing GSEA for test: ", test))
    fgsea_res <- suppressMessages(run_GSEA(geneList = geneList, pathways = pathways, test = test, outdir = test_dir))
    fgsea_list <- append(fgsea_list, list(fgsea_res))
  }

  # Run ssGSEA
  print(paste0("Performing ssGSEA"))

  ssGSEA_res <- run_ssGSEA(norm_data = norm_data, annots=df[[3]], pathways = pathways,
                           test = 'All Samples', outdir = outdir)

  # Add pathway analysis results to global environment
  names(ora_list) <- unique(de$test)
  assign("ora_list", ora_list, envir=.GlobalEnv)

  #names(ssGSEA_list) <- unique(de$test)
  #assign("ssGSEA_list", ssGSEA_list, envir=.GlobalEnv)
  names(fgsea_list) <- unique(de$test)
  assign("fgsea_list", fgsea_list, envir=.GlobalEnv)

  # Return the selected enrichment results
  #if (enrichment == "ORA") {
  #  enrich_res <- ora_list
  #}else if(enrichment == "ssGSEA"){
  #  enrich_res <- ssGSEA_list
  #} else {
  # Only support for gsea in downstream volcano plot
  enrich_res <- fgsea_list
  #}

  return(enrich_res)

}

# Print error messages regarding custom pathways files
custom_error <- function(custom_path, msg) {
    cat(paste0("Custom pathways directory, \"", custom_path, "\", ", msg, ".\n",
                 "Running default pathway analysis.\n"))
  return(NULL)
}

# Load custom .pathway files with gmt extension
load_gmts <- function(custom_gmt, list_of_genes) {
  # Load pathways from current gmt
  curr_path <- gmtPathways(custom_gmt)

  # Filter pathways by gene list
  filtered_pathways <- filter_custom_pathways(curr_pathway=curr_path,
                                                genes_list=list_of_genes)

  return(filtered_pathways)
}

# Convert symbols to entrez ids
symbol_to_entrez <- function(path_df, gene_col, species = "human") {

  if(tolower(species) == "mouse"){
    db <- org.Mm.eg.db
  }else{
    db <- org.Hs.eg.db
  }
  # Add entrez ids to custom pathway data frame
  path_df$entrez <- mapIds(db, as.character(path_df[[gene_col]]),
                             'ENTREZID', 'SYMBOL')

  # Find which entrez ids are NAs and drop them
  na_idx <- which(is.na(path_df$entrez))
  if(length(na_idx) != 0) {
    path_df <- path_df[-na_idx, ]
  }
  return(path_df)
}

# Look for pathways that include genes of interest
filter_custom_pathways <- function(curr_pathway, genes_list) {
  # Melt pathway lists into single gene data frame
  melted_path <- data.frame(melt(curr_pathway))
  names(melted_path) <- c("genes", "path")

  # If strings not numeric then need to convert symbols to entrez ids
  if (is.na(as.numeric(curr_pathway[[1]][[1]]))) {
    # Avoid mix up with expected entrez ID column genes
    names(melted_path)[names(melted_path) == "genes"] <- "symbol"
    melted_path <- symbol_to_entrez(path_df=melted_path, gene_col="symbol")
    # Make entrez ID the gene name for analysis
    names(melted_path)[names(melted_path) == "entrez"] <- "genes"
  }

  pathways_keep <- list()
  # Find pathways that include gene of interest
  for (gene in genes_list) {
    curr_keep <- do.call(c, list(melted_path[melted_path[,"genes"]==gene, "path"]))
    pathways_keep <- c(unlist(pathways_keep), unlist(curr_keep))
  }
  # Remove duplicate pathways
  pathways_keep <- unique(pathways_keep)
  # Filter pathway list by only those to keep
  found_paths <- curr_pathway[pathways_keep]

  # If strings not numeric then need entrez ids in filtered list
  if (is.na(as.numeric(found_paths[[1]][[1]]))) {
    for (kept_path in pathways_keep) {
      found_paths[[kept_path]] <-
        melted_path[melted_path[["path"]] == kept_path, "genes"]
    }
  }

  return(found_paths)
}
