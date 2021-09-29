#' Functions for pathway overlap and visualization
#'
#' @export pathway_overlap
#' @export pathway_intersect
#' @export pathway_chord

#' @title pathway_overlap
#'
#' Function used to observe ledge gene overlap for pathways from fgsea results
#'
#' @param pathways a subset of distinct fsgea pathway results
#' @param all boolean of whether or not to calculate pairwise overlaps for all ordered pathway pairs
#' @param species_data name of the species for which analysis applies, supports "human" or "mouse", defaults to "human"
#'
#' @return df containing the intersection of leading edge genes for selected pathways
#'
#' @examples
#'
#' pathways <- distinct(temp_res, pathway, .keep_all = TRUE)
#' pathways <- head(pathways[order(pathways$NES),], n = 10)
#' overlapOutput <- pathway_overlap(pathways, all = TRUE, species_data = "human")

# Pathways overlap function
# Built to handle fgsea enrichment results (preferably filtered)
# Add column of all genes in pathway to fgsea results if desired

pathway_overlap <- function(pathways,
                            all = TRUE,
                            species_data = "human",
                            for_plot = FALSE) {

  # If all genes columns is not present use a dummy column
  if (is.null(pathways$genes)) {
    pathways$genes <- c("")
    drop_genes <- T
    }

  # Select species for converting geneid to gene name
  if (species_data == "human"){
    data <- 'org.Hs.eg'
  }
  else if (species_data == "mouse"){
    data <- 'org.Mm.eg'
  }

  # Create list of pairs to run intersection on, either all or just unique pairs
  if (all) {overlaps <- expand.grid(pathways$pathway, pathways$pathway)}
  else {overlaps <- t(combn(pathways$pathway, 2))}

  # Apply intersect function
  check <- apply(overlaps, 1,
                 function(x){c(x[1], x[2],
                               list(intersect(unlist(pathways[pathways$pathway == x[1],]$genes),
                                              unlist(pathways[pathways$pathway == x[2],]$genes))),
                               list(intersect(unlist(pathways[pathways$pathway == x[1],]$leadingEdge),
                                              unlist(pathways[pathways$pathway == x[2],]$leadingEdge))),
                               pathways[pathways$pathway == x[1],]$genes,
                               pathways[pathways$pathway == x[2],]$genes,
                               pathways[pathways$pathway == x[1],]$leadingEdge,
                               pathways[pathways$pathway == x[2],]$leadingEdge)
                   })

  # Prep vectors for handling further parsing of the intersection output
  source <- c()
  target <- c()
  overlapGenes <- c()
  ledgeGenes <- c()
  numOverlap <- c()
  numLedgeOverlap <- c()

  percnt_source_genes <- c()
  percnt_target_genes <- c()

  percnt_source_ledge <- c()
  percnt_target_ledge <- c()

  # Parse results for each intersection
  for (row in check){

    # Calculate percentage overlap of intersection with either the source or target pathway
    if (drop_genes) {
      percnt_source_genes <- append(percnt_source_genes, " ")
      percnt_target_genes <- append(percnt_target_genes, " ")
    }
    else {
      percnt_source_genes <- append(percnt_source_genes, sprintf("%1.2f%%", 100*(length(row[[3]])/length(row[[5]]))))
      percnt_target_genes <- append(percnt_target_genes, sprintf("%1.2f%%", 100*(length(row[[3]])/length(row[[6]]))))
    }

    percnt_source_ledge <- append(percnt_source_ledge, sprintf("%1.2f%%", 100*(length(row[[4]])/length(row[[7]]))))
    percnt_target_ledge <- append(percnt_target_ledge, sprintf("%1.2f%%", 100*(length(row[[4]])/length(row[[8]]))))

    # Convert geneids to gene names, prepare output for inclusion in hover over output using breaks <br>
    if (length(row[[3]]) == 0) {genenames <- c()}
    else {
      genenames <- unname(getSYMBOL(row[[3]], data=data))
      #genenames <- genenames[!is.na(genenames)]
      genenames[is.na(genenames)] <- "NA"
    }

    if (all) {
      genenames <- paste(genenames, collapse = " ")
      if (for_plot){
        genenames <- gsub("(.{21,}?)\\s", "\\1<br>", genenames)
      }
    }

    if (length(row[[4]]) == 0) {ledgenames <- c()}
    else {
      ledgenames <- unname(getSYMBOL(row[[4]], data=data))
      #ledgenames <- ledgenames[!is.na(ledgenames)]
      ledgenames[is.na(ledgenames)] <- "NA"
    }

    if (all) {
      ledgenames <- paste(ledgenames, collapse = " ")
      if (for_plot){
        ledgenames <- gsub("(.{21,}?)\\s", "\\1<br>", ledgenames)
      }

    }

    # Append final results to various vectors
    source <- append(source, row[[1]])
    target <- append(target, row[[2]])
    overlapGenes <- append(overlapGenes, list(genenames))
    ledgeGenes <- append(ledgeGenes, list(ledgenames))
    numOverlap <- append(numOverlap, length(row[[3]]))
    numLedgeOverlap <- append(numLedgeOverlap, length(row[[4]]))
  }

  # Build df from output vectors depending on the exisitence of an initial genes column
  if (drop_genes) {
    output <- data.frame("source" = source,
                         "target" = target,
                         "ledgeGenes" = matrix(ledgeGenes, byrow = T),
                         "numLedgeOverlap" = numLedgeOverlap,
                         "percntSourceOverlapLedge" = percnt_source_ledge,
                         "percntTargetOverlapLedge" = percnt_target_ledge)
  }

  else {
    output <- data.frame("source" = source,
                       "target" = target,
                       "overlapGenes" = matrix(overlapGenes, byrow = T),
                       "ledgeGenes" = matrix(ledgeGenes, byrow = T),
                       "numOverlap" = numOverlap,
                       "numLedgeOverlap" = numLedgeOverlap,
                       "percntSourceOverlapGenes" = percnt_source_genes,
                       "percntTargetOverlapGenes" = percnt_target_genes,
                       "percntSourceOverlapLedge" = percnt_source_ledge,
                       "percntTargetOverlapLedge" = percnt_target_ledge)
  }

  return(output)
  }


#' @title pathway_intersect
#'
#' Function used to plot the ledge gene overlap for pathways from fgsea results
#'
#' @param pathways a subset of distinct fsgea pathway results
#' @param species_data name of the species for which analysis applies, supports "human" or "mouse", defaults to "human"
#' @param outfile name of file for output
#' @param save_plot boolean of whether to save the plot as a file or just output to the console, defaults to TRUE and uses config information for file format
#'
#' @return image file of pathway intersection plot
#' @return plotly interactive pathway intersection plot
#'
#' @examples
#'
#' pathway_intersect(pathways, species_data = "human", outfile = paste0('Pathways_Intersect.', fileType), save_plot = TRUE)
#'

# Pathways intersection plot function
# Currently hardcoded for leading edge genes only but functionality can easily be added for all genes
# Uses all = TRUE for pathway_overlap function output

pathway_intersect <- function(pathways,
                              species_data = "human",
                              outfile = 'Pathways_Intersect',
                              save_plot = TRUE) {

  output <- pathway_overlap(pathways, all = TRUE, species_data = species_data, for_plot = TRUE)

  # Use the order of enrichment for plotting
  order <- unique(output$source)

  # Build symetric matrix for overlap
  output_wide <- spread(output[c("source", "target", "numLedgeOverlap")], target, numLedgeOverlap)
  output_wide <- column_to_rownames(output_wide, var = "source")
  output_wide <- output_wide[order, order]

  # Function to "unspread" symetric matrix
  gather.matrix <- function(mat) {
    if (is.null(dimnames(mat))) {
      grid <- expand.grid(seq.int(nrow(mat)), seq.int(ncol(mat)))
    } else {
      grid <- expand.grid(dimnames(mat))
    }
    cbind(grid, value = as.vector(mat))
  }

  # Build df for plotting numeric overlap with hover over information
  output_low <- output_wide
  output_low[upper.tri(output_low, diag = TRUE)] <- NA
  output_low <- gather.matrix(data.matrix(output_low))
  names(output_low) <- c('source', 'target', 'numLedgeOverlap')
  output_low$numLedgeOverlapAll <- output$numLedgeOverlap
  output_low$ledgeGenes <- output$ledgeGenes
  output_low$percntSourceOverlapLedge <- output$percntSourceOverlapLedge
  output_low$percntTargetOverlapLedge <- output$percntTargetOverlapLedge

  # Build df for plotting overlap circles with hover over information
  output_high <- output_wide
  output_high[lower.tri(output_high, diag = FALSE)] <- NA
  output_high <- gather.matrix(data.matrix(output_high))
  names(output_high) <- c('source', 'target', 'numLedgeOverlap')
  output_high$numLedgeOverlapAll <- output$numLedgeOverlap
  output_high$ledgeGenes <- output$ledgeGenes
  output_high$percntSourceOverlapLedge <- output$percntSourceOverlapLedge
  output_high$percntTargetOverlapLedge <- output$percntTargetOverlapLedge

  # Build plot using ggplot2, warning suppresed since text hover over information not used in ggplot
  p <- suppressWarnings(ggplot(output_high, aes(y = source, x = target)) +
    geom_point(aes(colour = numLedgeOverlap,
                   size = numLedgeOverlap,
                   text = sprintf("Number Overlap: %s\nPercent Source: %s\nPercent Target: %s\nGenes: %s",
                                  numLedgeOverlapAll, percntSourceOverlapLedge, percntTargetOverlapLedge, ledgeGenes))) +
    scale_color_gradientn(colours = viridis(10, option = "D")) +
    labs(colour = "# Overlap") +
    theme_bw() +
    geom_text(data = output_low, aes(y = source,
                                     x = target,
                                     label = numLedgeOverlap,
                                     text = sprintf("Number Overlap: %s\nPercent Source: %s\nPercent Target: %s\nGenes: %s",
                                                    numLedgeOverlapAll, percntSourceOverlapLedge, percntTargetOverlapLedge, ledgeGenes))) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)))
  # Save ggplot to file

  if(save_plot) {
    if(fileType %in% c('tiff','png','jpeg','bmp')) {
      match.fun(fileType)(outfile, width = 1400, height = 1000, res = 150)
    } else {
      match.fun(fileType)(outfile)
    }
    suppressWarnings(print(p))
    dev.off()
  }

  # Convert ggplot results to interactive plotly with hover over information
  g <- ggplotly(p, tooltip="text")
  print(g)

  return(g)
}

#' @title pathway_chord
#'
#' Function used to plot a chord diagram of the ledge gene overlaps for pathways from fgsea results
#'
#' @param pathways a subset of distinct fsgea pathway results
#' @param species_data name of the species for which analysis applies, supports "human" or "mouse", defaults to "human"
#' @param outfile name of file for output
#' @param save_plot boolean of whether to save the plot as a file or just output to the console, defaults to TRUE and uses config information for file format
#'
#' @return image file of chord diagram pathway plot
#'
#' @examples
#'
#' pathway_chord(pathways, species_data = "human", outfile = paste0('GSEA_Pathways_Chord.', fileType), save_plot = TRUE)
#'

# Pathways chord diagram function
# Currently hardcoded for leading edge genes only but functionality can easily be added for all genes
# Uses all = TRUE for pathway_overlap function output

pathway_chord <- function(pathways,
                          species_data = "human",
                          outfile = 'Pathways_Chord',
                          save_plot = TRUE) {

  output <- pathway_overlap(pathways, all = TRUE, species_data = species_data, for_plot = TRUE)

  output_chord <- output

  # Keep enrichment order in tact for plotting
  order <- unique(output$source)

  # Convert output to symetric matrix
  output_chord <- spread(output_chord[c("source", "target", "numLedgeOverlap")], target, numLedgeOverlap)
  output_chord <- column_to_rownames(output_chord, var = "source")

  # Reorder based on enrichment order
  output_chord <- output_chord[order, order]

  # Add enrichment order to pathway names
  col_names <-c()

  i <- 1
  for (name in colnames(output_chord)) {
    num <- paste("#", i, sep = "")
    col_names <- c(col_names, paste(num, name, sep = " "))
    i <- i + 1
  }

  colnames(output_chord) <- col_names
  rownames(output_chord) <- col_names

  # Build copy dataframe with labels to be used in the chord diagram
  nums <- paste("#", 1:length(colnames(output_chord)), sep = "")

  output_chord_diag <- output_chord
  colnames(output_chord_diag) <- nums
  rownames(output_chord_diag) <- nums

  # Set color scheme
  mycolor <- viridis(length(colnames(output_chord)), alpha = 1, begin = 0, end = 1, option = "D")
  names(mycolor) <- colnames(output_chord_diag)
  # mycolor <- mycolor[c(c(2:length(mycolor)),1)]

  # Uses R base plot, viewports used to adjust where the chord diagram is plotted leaving room for legend

  plot_the_chords <- function() {

    plot.new()
    circle_size <- unit(0.7, 'npc')
    pushViewport(viewport(x = unit(0, "npc"), y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)

    chordDiagram(data.matrix(output_chord_diag),
                 order = union(rownames(output_chord_diag), colnames(output_chord_diag)),
                 grid.col = mycolor,
                 symmetric = TRUE,
                 link.sort = TRUE,
                 transparency = 0.6,
                 annotationTrack = "grid",
                 preAllocateTracks = 1)

    # Plot chord diagram labels in preAllocatedTrack
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.2, CELL_META$sector.index,
                  facing = "outside", niceFacing = TRUE, cex = 1.2, adj = c(0, 0.5))
    }, bg.border = NA)

    # Define legend
    l <- Legend(at = colnames(output_chord),
           legend_gp = gpar(fill = mycolor[match(nums, get.all.sector.index())]), title_position = "topleft",
           title = "Top Pathways")

    # Move out of chord viewport for plotting legend
    upViewport()
    pushViewport(viewport(x = unit(0.3, "npc"), y = 0.5, width = circle_size, height = circle_size, just = c("left", "center"), clip = "off"))
    grid.draw(l)
    upViewport()

  }

  if(save_plot) {
    if(fileType %in% c('tiff','png','jpeg','bmp')) {
      match.fun(fileType)(outfile, width = 2250, height = 1000, res = 150)
    } else {
      match.fun(fileType)(outfile)
    }
    plot_the_chords()
    dev.off()
  }

  plot_the_chords()

}

