#' @title convert_target_notes
#'
#' Converts target notes from Erin to DSP-DA format
#'
#' @param targetNotes target notes data frame
#'
#' @return targetNotes data frame formatted like DSP-DA
#'
#' @examples
#' targetNotes <- convert_target_notes(targetNotes)
#'
#' @export convert_target_notes

#library(reshape2) <- included in DSP_hydra.R

convert_target_notes <- function(targetNotes) {
  # Check if conversion is necessary
  if (ncol(targetNotes) == 2 && !('Gene' %in% colnames(targetNotes))) {
    # when set to header the probe group has spaces replaced by '.'
    # so fix that
    row_1 <- gsub('\\.', ' ', colnames(targetNotes))
    colnames(targetNotes) <- c('Group', 'Target')

    # Flip to strings incase either value in the header was unique
    targetNotes <- data.frame(apply(targetNotes, 2, as.character),
                           stringsAsFactors = FALSE)
    targetNotes <- rbind(row_1, targetNotes)

    # Replace all the weird symbols allowed in probe group names
    targetNotes$Group <- gsub('-|[/ ,]+', '.', targetNotes$Group)

    # Reshape
    targetNotes <- dcast(targetNotes, Target ~ Group, value.var = 'Target')

    # Store Gene Info
    Gene <- targetNotes$Target

    # Replace NAs and not NAs with appropriate symbol
    targetNotes <- targetNotes[,2:ncol(targetNotes)]
    targetNotes[!is.na(targetNotes)] <- '+'
    targetNotes[is.na(targetNotes)] <- '-'

    # Get Cell.Type Info
    # The below list comes from Erin and is said to be static for our purposes
    celltypingDf <- data.frame(
      'Gene'= c("BLK", "CD19", "MS4A1", "TNFRSF17",
      "FCRL2", "KIAA0125", "FAM30A", "PNOC", "SPIB", "TCL1A", "PTPRC", "CD8A",
      "CD8B", "CTSW", "GNLY", "GZMA", "GZMB", "GZMH", "KLRB1", "KLRD1",
      "KLRK1", "PRF1", "NKG7", "CCL13", "CD209", "HSD11B1", "CD244", "EOMES",
      "LAG3", "PTGER4", "CD163", "CD68", "CD84", "MS4A4A", "MS4A2", "TPSAB1",
      "CPA3", "HDC", "TPSB2", "CSF3R", "S100A12", "CEACAM3", "FCAR", "FCGR3A",
      "FCGR3B", "FPR1", "SIGLEC5", "IL21R", "KIR2DL3", "KIR3DL1", "KIR3DL2",
      "NCR1", "XCL2", "XCL1", "CD3D", "CD3E", "CD3G", "CD6", "SH2D1A",
      "TRAT1", "TBX21", "FOXP3", "Blk", "Cd19", "Ms4a1", "Tnfrsf17", "Fcrlb",
      "Pnoc", "Spib", "Tcl1", "Ptprc", "Cd8a", "Cd8b1", "Ctsw", "Gzma", "Gzmb",
      "Klrb1", "Klrd1", "Klrk1", "Prf1", "Nkg7", "Ccl2", "Cd209e", "Hsd11b1",
      "Cd244", "Cd244a", "Eomes", "Lag3", "Ptger4", "Cd163", "Cd68", "Cd84",
      "Ms4a4a", "Ms4a2", "Tpsb2", "Tpsab1", "Cpa3", "Hdc", "Csf3r", "Ceacam3",
      "Fcgr4", "Fpr1", "Il21r", "Kir3dl1", "Kir3dl2", "Ncr1", "Xcl1", "Cd3d",
      "Cd3e", "Cd3g", "Cd6", "Sh2d1a", "Trat1", "Tbx21", "Foxp3"),

      'Cell.Type' = c("B-cells", "B-cells", "B-cells", "B-cells", "B-cells",
      "B-cells", "B-cells", "B-cells", "B-cells", "B-cells", "CD45",
      "CD8 T cells", "CD8 T cells", "Cytotoxic cells", "Cytotoxic cells",
      "Cytotoxic cells", "Cytotoxic cells", "Cytotoxic cells",
      "Cytotoxic cells", "Cytotoxic cells", "Cytotoxic cells",
      "Cytotoxic cells", "Cytotoxic cells", "DC", "DC", "DC", "Exhausted CD8",
      "Exhausted CD8", "Exhausted CD8", "Exhausted CD8", "Macrophages",
      "Macrophages", "Macrophages", "Macrophages", "Mast cells", "Mast cells",
      "Mast cells", "Mast cells", "Mast cells", "Neutrophils", "Neutrophils",
      "Neutrophils", "Neutrophils", "Neutrophils", "Neutrophils",
      "Neutrophils", "Neutrophils", "NK CD56dim cells", "NK CD56dim cells",
      "NK CD56dim cells", "NK CD56dim cells", "NK cells", "NK cells",
      "NK cells", "T-cells", "T-cells", "T-cells", "T-cells", "T-cells",
      "T-cells", "Th1 cells", "Treg", "B-cells", "B-cells", "B-cells",
      "B-cells", "B-cells", "B-cells", "B-cells", "B-cells", "CD45",
      "CD8 T cells", "CD8 T cells", "Cytotoxic cells", "Cytotoxic cells",
      "Cytotoxic cells", "Cytotoxic cells", "Cytotoxic cells",
      "Cytotoxic cells", "Cytotoxic cells", "Cytotoxic cells", "DC", "DC",
      "DC", "Exhausted CD8", "Exhausted CD8", "Exhausted CD8",
      "Exhausted CD8", "Exhausted CD8", "Macrophages", "Macrophages",
      "Macrophages", "Macrophages", "Mast cells", "Mast cells",
      "Mast cells", "Mast cells", "Mast cells", "Neutrophils", "Neutrophils",
      "Neutrophils", "Neutrophils", "NK CD56dim cells", "NK CD56dim cells",
      "NK CD56dim cells", "NK cells", "NK cells", "T-cells", "T-cells",
      "T-cells", "T-cells", "T-cells", "T-cells", "Th1 cells", "Treg"),
      stringsAsFactors = FALSE
    )

    Cell.Type <- celltypingDf$Cell.Type[match(Gene, celltypingDf$Gene)]
    Cell.Type[is.na(Cell.Type)] <- ''

    # Add gene info back
    targetNotes <- cbind(Gene, Cell.Type, targetNotes)

  }

  return(targetNotes)
}

