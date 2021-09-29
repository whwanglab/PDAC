#' @title run_ora
#'
#' Pathway Analysis helper function to run Over Representation Analysis
#'
#' @param geneList a list of ranked genes, the way that they are ranked is determined in run_pathways.R
#' @param pathways are the genesets that we are working with, default is Reactome
#' @param test is the test or contrast that is being compared for pathway analysis, from DE results
#' @param path_fc the fold change cutoff to use for gene list
#' @param path_pval the pvalue to threshold the enrichment results by
#' @param outdir the output directory for plots
#'
#' @return ora The results from over representation analysis
#'
#' @examples
#'  run_ora(geneList = geneList, pathways = pathways, test = test, path_fc = path_fc, path_pval = path_pval, outdir = outdir)
#'
#' @export run_ora

run_ora <- function(geneList = geneList, pathways = pathways, test = test, path_fc = path_fc, path_pval = path_pval, outdir = outdir, sp = species){

  # Grab the top genes
  top_genes <- names(geneList)[abs(geneList) > path_fc]

  # Threshold
  if(length(top_genes) < 50){
    print("Top genes defined by Fold Change is not more than 50, will not do ORA pathway Analysis")
  }else {

    # Reformat to accommodate consistent pathways used across the different analyses
    melted_df <- data.frame(melt(pathways))
    names(melted_df) <- c("entrez", "path")

    # Create two column data frame for builing annotation environment
    path_to_entrez <- melted_df[,c("path", "entrez")]

    # Build annotations environment for dose package
    anno_env <- DOSE:::build_Anno(path2gene=path_to_entrez)

    if(tolower(sp) == "mouse"){
      db <- org.Mm.eg.db
    }else{
      db <- org.Hs.eg.db
    }

    # Run ORA
    ora <- DOSE:::enricher_internal(gene=top_genes, pvalueCutoff=path_pval, USER_DATA=anno_env)
    if(as.integer(dim(ora)[1]) < 1){
      print("No significantly over-represented pathways found in ORA.")
      return(ora)
    }
    ora <- DOSE:::setReadable(ora, db, keyType="ENTREZID")

    # Make the enrichment plots
    p1 <- barplot(ora, x = "GeneRatio", showCategory = 8, title = "Overall Coverage")
    p2 <- dotplot(ora, x = "GeneRatio", showCategory = 8, title = "Overall Coverage with Gene Ratios", orderBy="GeneRatio")
    p3 <- emapplot(ora)
    p4 <- cnetplot(ora, categorySize = "pvalue", foldChange = geneList)

    # Make Plots and Save Files
    width <- 700
    height <- 500

    if (fileType == "pdf" | fileType == "svg"){
      width <- width/85
      height <- height/85
    }

    # Save to file
    match.fun(fileType)(paste0(outdir, "/", test, "_Overall_Coverage_ORA.", fileType), width = width, height = height)
    print(p1)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/", test, "_Overall_Coverage_ORA_with_Ratios.", fileType), width = width, height = height)
    print(p2)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/", test, "_Overall_Coverage_ORA_Map.", fileType), width = width, height = height)
    print(p3)
    dev.off()

    match.fun(fileType)(paste0(outdir, "/", test, "_Mulitple_Annot_ORA.", fileType), width = width, height = height)
    print(p4)
    dev.off()

    ora <- as.data.frame(ora)
    ora$ID <- gsub(pattern = "\n", replacement = " ", ora$ID)

    # Write to file
    write.csv(ora, file=paste0(outdir, "/", test, "_ORA_results.csv"), row.names = F)

    return(ora)
  }
}
