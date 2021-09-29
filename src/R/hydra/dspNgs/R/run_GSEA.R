#' @title run_GSEA
#'
#' Pathway Analysis helper function to run GSEA
#'
#' @param geneList a list of ranked genes, the way that they are ranked is determined in run_pathways.R
#' @param pathways are the genesets that we are working with, default is Reactome
#' @param test is the test or contrast that is being compared for pathway analysis, from DE results
#' @param outdir the output directory for plots
#'
#' @return fgseaRes The results from GSEA analysis
#'
#' @examples
#'  run_GSEA(geneList = geneList, pathways = pathways, test = test, outdir = outdir)
#'
#' @export run_GSEA
#' @export entrez_to_symbol

run_GSEA <- function(geneList = geneList, pathways = pathways, test = test, outdir = outdir){

  # Run GSEA
  fgseaRes <- fgsea(pathways, geneList, nperm=1000, minSize = 15, maxSize=500)
  fgseaRes <- as.data.frame(fgseaRes)

  # Grab top pathways up and down
  fgseaRes <- fgseaRes[which(!is.na(fgseaRes$pval)),] # removes NAs

  topPathwaysUp <- fgseaRes[fgseaRes$ES > 0,]
  topPathwaysUp <- topPathwaysUp[order(topPathwaysUp$padj, -abs(topPathwaysUp$NES)), ]
  topPathwaysUp <- topPathwaysUp$pathway[1:10]

  topPathwaysDown <- fgseaRes[fgseaRes$ES < 0,]
  topPathwaysDown <- topPathwaysDown[order(topPathwaysDown$padj, abs(topPathwaysDown$NES)), ]
  topPathwaysDown <- topPathwaysDown$pathway[1:10]

  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  # Make Classic GSEA enrichment Plots
  p5 <- plotEnrichment(pathways[[topPathwaysUp[1]]], geneList) + labs(title = topPathwaysUp[1])
  p6 <- plotEnrichment(pathways[[topPathwaysDown[1]]], geneList) + labs(title = topPathwaysDown[1])

  # Make a plot to look at only independent pathways
  collapsedPathways <- fgseaRes[order(fgseaRes$padj), ]
  collapsedPathways <- collapsePathways(collapsedPathways[collapsedPathways$padj < 0.01, ], pathways, geneList)
  mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways$mainPathways,]
  mainPathways <- mainPathways[order(-(mainPathways$NES)), ]$pathway

  # Make Plots and Save Files
  width <- 900
  height <- 500

  if (fileType == "pdf" | fileType == "svg"){
    width <- width/85
    height <- height/85
  }

  # Save to file
  match.fun(fileType)(paste0(outdir, "/", test, "_Top_Pathway_Up.", fileType), width = width, height = height)
  print(p5)
  dev.off()

  match.fun(fileType)(paste0(outdir, "/", test, "_Top_Pathway_Down.", fileType), width = width, height = height)
  print(p6)
  dev.off()

  if(length(topPathways) > 0){
    match.fun(fileType)(paste0(outdir, "/", test, "_Top_Pathways_Table.", fileType), width = width, height = height)
    plotGseaTable(pathways[topPathways], geneList, fgseaRes, gseaParam = 0.5,colwidths=c(10, 3, 0.8, 0, 1))
    dev.off()
  }

  if(length(mainPathways) > 0){
    match.fun(fileType)(paste0(outdir, "/", test, "_Main_Independent_Pathways.", fileType), width = width, height = height)
    plotGseaTable(pathways[mainPathways], geneList, fgseaRes, gseaParam = 0.5,colwidths=c(10, 3, 0.8, 0, 1))
    dev.off()

    write.table(as.data.frame(mainPathways), file=paste0(outdir, "/", test, "_MainPathways.tsv"), row.names = FALSE, sep = "\t")
  }

  # Write to file
  write.table(as.matrix(fgseaRes), file=paste0(outdir,"/", test, "_GSEA_results.tsv"), row.names = FALSE, sep = "\t")

  fgseaRes2 <- entrez_to_symbol(fgseaRes)

  write.table(as.matrix(fgseaRes2)[,c(1:7, 9)], file=paste0(outdir,"/", test, "_GSEA_results_names.tsv"), sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)

  return(fgseaRes)
}


#' @title entrez_to_symbol
#'
#' Function to change entrez ID to gene symbol
#'
#' @param X GSEA dataframe with leading edge data
#' @param sp name of the species for which analysis applies, supports "human" or "mouse", defaults to species in config
#'
#' @return gene symbols
#'
#' @examples
#'
#' entrez_to_symbol(X = fgseaRes, sp = species)
#'
entrez_to_symbol <- function(X, sp = species){
  if(tolower(sp) == "mouse"){
    db <- org.Mm.eg.db
  }else{
    db <- org.Hs.eg.db
  }

  X$leadingEdgeName <- sapply(X$leadingEdge, function(x){
    return(paste0(as.character(mapIds(db, x, 'SYMBOL','ENTREZID')), collapse=";"))
  })
  return(X)
}
