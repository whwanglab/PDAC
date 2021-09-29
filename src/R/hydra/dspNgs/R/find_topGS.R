#' @title find_topGS
#'
#' Generic call to heatmap functionality
#'
#' @param targetNotes target annotations - list of lists of genes output from parse_GeneSets
#' @param de_results table of DE results
#' @param width file width
#' @param height file height
#' @param name annotations to include on heatmap
#' @param outdir where to print file
#' @param fileType function to call for file printing
#'
#' @return
#'
#' @examples
#'
#'
#' @export find_topGS





# top GeneSets
find_topGS <- function(targetNotes = NULL,
                       de_results = NULL,
                       name = NULL,
                       outdir = NULL,
                       fileType = NULL,
                       width = 1500,
                       height = 850) {
  GS_hits <- data.frame(row.names = names(targetNotes))
  GS_hits$FC = GS_hits$pval = GS_hits$both <- NA
  for(gs in names(targetNotes)) {
    gs_DE <- subset(de_results, gene %in% targetNotes[[gs]])
    GS_hits[gs, 'FC'] <- sum(abs(gs_DE$FC) > fc_cutoff & gs_DE$Pval > pval_cutoff)
    GS_hits[gs, 'pval'] <- sum(abs(gs_DE$FC) < fc_cutoff & gs_DE$Pval < pval_cutoff)
    GS_hits[gs, 'both'] <- sum(abs(gs_DE$FC) > fc_cutoff & gs_DE$Pval < pval_cutoff)
  }
  GS_hits$all <- rowSums(GS_hits)
  GS_hits$nGenes <- as.numeric(lapply(targetNotes, length))
  GS_hits_frac <- GS_hits / GS_hits$nGenes
  GS_hits_frac <- GS_hits_frac[order(GS_hits_frac$all, decreasing = TRUE), ]
  GS_hits_frac$GeneSet <- rownames(GS_hits_frac)
  GS_hits_m <- melt(GS_hits_frac, id.vars = 'GeneSet', measure.vars = c('both','pval','FC'), variable.name = 'ModeHit', value.name = 'Percentage')

  # make plot
  bp <- ggplot(data = subset(GS_hits_m, GeneSet %in% rownames(GS_hits_frac)[1:10]),
               aes(y = Percentage*100,
                   x = reorder(GeneSet, Percentage),
                   fill = ModeHit)) +
    geom_bar(position="stack", stat="identity") +
    coord_flip() + scale_fill_brewer(palette = 'Dark2', name = 'Identified by') +
    labs(y = 'Genes Detected, %', x = '',
         title = paste0('Top DE Annotation Sets: ', name)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))
  if (is.null(preset_theme)) {
    bp <- bp + theme_light(base_size = 16)
  }

  # print plot
  match.fun(fileType)(paste0(outdir, "/DEGeneSetSummary.", name, '.',
                             fileType),
                      width = width,
                      height = height, res = 130)
  print(bp)
  dev.off()

  # return data table of GS_hits_frac
  return(GS_hits_frac)
}
