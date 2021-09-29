#' DSP-NGS Sequencing QC
#'
#' Generates sequecing QC metrics for provided dataset
#'
#' @param dataframelist List of all 5 dataframes
#' @param loq Number of standard deviations desired for loq. Column must exist as GeoLOQ|SD|_|Pool| for each pool in the dataset, ex. GeoLOQ2.5_01. Default = 2.5
#' @param qc_dir Directory to output plot + QC table. Default = './'
#' @param facet_column Column to facet the line plot on if desired. Default = No facet.
#' @param fileType Type of file for ggsave to generate for the QC Line plot.
#' @param return_plots Boolean. If TRUE, a list of QC plots will be returned.
#'
#' @return Write QC metrics to file in QC directory
#' @return Write QC plots to file QC directory (QC of reads, QC of sequencing saturation)
#' @return QC plots as list if return_plots = TRUE.
#'       list('reads' = QC of reads, 'saturation' = QC of sequencing saturation)
#'

#'
#' @examples
#'    run_seqQC(dataframelist=dfs, loq=2.5, fileType='pdf')
#' @export run_seqQC
#'

run_seqQC <- function(dataframelist, loq=2.5, qc_dir='./', facet_column='',
                      fileType='pdf', return_plots = FALSE) {

  annot <- dataframelist[[3]]

  # Names for legibility
  names(dataframelist) <- c('TargetCountMatrix', 'NegNorm_TargetCountMatrix',
                            'SegmentProperties', 'DatasetHistory',
                            'TargetProperties')

  # Check the pooling information was read in as character
  if (class(dataframelist[['TargetProperties']]$Pooling) != "character") {
    stop(paste('Pooling information in the target properties sheet was not',
               'parsed correctly. Pooling should be a character column, but',
               'was read in as a',
               class(dataframelist[['TargetProperties']]$Pooling)), '.')
  }

  # The idColumn in the segment properties sheet might change in the DSPDA export
  # pre-building this in to easily adapt to a change
  idColumn <- 'Sample_ID'
  segment_properties <- dataframelist[['SegmentProperties']]
  colnames(segment_properties)[which(colnames(segment_properties) == idColumn)] <- 'Sample_ID'

  # Make things easy to work with
  rownames(dataframelist[['TargetCountMatrix']]) <- dataframelist[['TargetCountMatrix']][,1]
  dataframelist[['TargetCountMatrix']] <- dataframelist[['TargetCountMatrix']][,-1]
  rownames(dataframelist[['TargetProperties']]) <- dataframelist[['TargetProperties']][,'TargetName']

  # Get matching orders by sample ID
  order <- match(dataframelist[['SegmentProperties']]$Sample_ID,
                 colnames(dataframelist[['TargetCountMatrix']]))
  dataframelist[['TargetCountMatrix']] <- dataframelist[['TargetCountMatrix']][,order]
  # Remove any probes not in target properties (e.g. NegProbe in DSP-DA data)
  dataframelist[['TargetCountMatrix']] <- dataframelist[['TargetCountMatrix']][rownames(dataframelist[['TargetProperties']]), ]

  # LOQ Calculations
  loqCalcs <- lapply(rownames(dataframelist[['TargetCountMatrix']]), function(x) {
    # Get counts for target
    row <- dataframelist[['TargetCountMatrix']][x,]

    # Get Pool Name and Column
    pool <- dataframelist[['TargetProperties']][x,'Pooling']
    poolcol <- paste0('GeoLOQ', loq, '_', pool)

    loqs <- dataframelist[['SegmentProperties']][,poolcol]
    loqs[loqs < LOQ_floor] <- LOQ_floor
    loqtest <- row > t(loqs)

    return(loqtest)
  })
  loqCalcs <- data.frame(matrix(unlist(loqCalcs), nrow=length(loqCalcs), byrow=T))

  # Summarize LOQ information
  genesOverLOQ <- length(which(apply(loqCalcs, 1, any)))
  genesPerAOI <- apply(loqCalcs, 2, function(x) { return(length(which(x))) })
  names(genesPerAOI) <- segment_properties$Sample_ID
  medianSigGenes <- median(apply(loqCalcs, 2, function(x) { return(length(which(x))) }))

  # Calculate generic sequencing QC Metrics
  qc_vals <- list(
    'Segments' = format(nrow(segment_properties), big.mark = ','),
    'Median Raw Reads per Segment' = format(median(segment_properties$RawReads), big.mark = ','),
    'Median UMIs per Segment' = format(median(segment_properties$DeduplicatedReads),
                           big.mark = ','),
    'Median Genes > LOQ per Segment' = format(medianSigGenes, big.mark = ','),
    'Genes > LOQ in at least one Segment' = format(genesOverLOQ, big.mark = ','),
    'Total Raw' = format(sum(segment_properties$RawReads), big.mark = ','),
    'Total Mapped to RTS' = paste0(
      round(sum(segment_properties$AlignedReads) /
              sum(segment_properties$RawReads) * 100), '%'),
    'Percent of Mapped that are Unique' = paste0(
      round(sum(segment_properties$DeduplicatedReads) /
              sum(segment_properties$AlignedReads) * 100), '%'),
    'Sequencing Saturation' = paste0(
      round( (1 - sum(segment_properties$DeduplicatedReads) /
              sum(segment_properties$AlignedReads)) * 100), '%')
  )

  cat(paste0(names(qc_vals), ': ', qc_vals, collapse = '\n'),
      file=paste(qc_dir, 'seqQC_Values.txt', sep='/'))

  cat(paste0(names(genesPerAOI), ',', genesPerAOI, collapse = '\n'),
      file=paste(qc_dir, 'genesPerAOI.txt', sep='/'))

  # List of columns for the line plot
  measures <- c('RawReads', 'TrimmedReads', 'StitchedReads', 'AlignedReads',
                'DeduplicatedReads')

  # Melt data frame and create a facet column if desired
  if (facet_column != '') {
    sumdf <- melt(segment_properties, id.vars=c('Sample_ID', facet_column),
                  measure.vars = measures)
    colnames(sumdf)[2] <- 'Facet'
  } else {
    sumdf <- melt(segment_properties, id.vars=c('Sample_ID'),
                  measure.vars = measures)
  }

  # Reformat to make things look better
  sumdf$variable <- gsub('Reads', ' Reads', sumdf$variable)
  sumdf$variable <- factor(sumdf$variable,
                           levels = gsub('Reads', ' Reads', measures))

  # Flow style log 10 y axis breaks and Read Count labels
  yBreaks <- unique(c(1:10 %o% 10^(0:7)))
  yLabs <- formatC(yBreaks, format = 'd', big.mark = ',')
  yLabs[-grep(1, yLabs)] <- ''
  yLabs <- gsub(',000,000', 'M', yLabs)
  yLabs <- gsub(',000', 'K', yLabs)

  # QC Line Plot
  qcLine <- ggplot(sumdf, aes(x=Sample_ID, y=as.numeric(value), color=variable, group=variable, linetype=variable)) +
    geom_line(size=1.2, alpha=0.8) +
    scale_y_log10(breaks=yBreaks, labels=yLabs) +
    scale_color_manual(values=brewer.pal(5, 'Set1')) +
    labs(y='Reads')
    if (is.null(preset_theme)) {
      qcLine <- qcLine + theme_bw()
    }
    qcLine <- qcLine +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, size=6),
            strip.background = element_rect(fill='white'),
            axis.title.x = element_blank(),
            legend.title=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size=6))

  # Facet if necessary and determine figure dimensions
  if (facet_column != '') {
    qcLine <- qcLine + facet_wrap(~Facet, scales='free', ncol=1)
    maxFacetWidth <- max(table(segment_properties[facet_column])) / 11
    numFacets <- length(unique(sumdf$Facet))
  } else {
    maxFacetWidth <- nrow(segment_properties) / 11
    numFacets <- 1
  }

  # Save the plot
  ggsave(paste(qc_dir, paste0('seqQC_linePlot.', fileType), sep='/'),
         qcLine, device=fileType, dpi=320,
         width=max(maxFacetWidth, 5),
         height=3*numFacets, units="in", limitsize = FALSE, scale = 1)

  # qcSeqSat plot
  annot[["SequencingSaturation"]] <- as.numeric(annot[["SequencingSaturation"]])
  annot <- annot[order(annot$SequencingSaturation),]

  qcSeqSat <- ggplot(annot,
                     aes(x = reorder(Sample_ID, SequencingSaturation),
                         y = SequencingSaturation,
                         group = NA)) +
    geom_line(size = 1, color = 'black') +
    labs(x = 'Sample', y = 'Sequencing saturation')
    if (is.null(preset_theme)) {
      qcSeqSat <- qcSeqSat + theme_minimal(base_size = 12)
    }
    qcSeqSat <- qcSeqSat + theme(panel.grid = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.line.x = element_line(),
                                   axis.line.y = element_line())

  ggsave(paste(qc_dir, paste0('seqQC_seqSat_line.', fileType), sep='/'),
           qcSeqSat, device=fileType, dpi=300, width=4.5, height=3, limitsize = FALSE, scale = 1)

  #Output a list of AOIs with low sequencing saturation, based on user-set cutoff, as a warning, samples not removed from analysis
  if(length(annot$Sample_ID[annot$SequencingSaturation < sat_cutoff])>0) {
    loSeqSat_warn_df <-
      data.frame(LoSeqSaturation_Warning = annot$Sample_ID[annot$SequencingSaturation < sat_cutoff],
                   SequencingSaturation = annot$SequencingSaturation[annot$SequencingSaturation < sat_cutoff])
    write.table(loSeqSat_warn_df, paste0(qc_dir, '/LowSeqSaturation_warning.txt'), sep="\t", row.names=FALSE)
    print(loSeqSat_warn_df)
  }

  # produce scatter plot to investigate seq saturation as a function of AOI size or cell count, etc.
  # color by defined grouping
  annot[[sat_numcrit]] <- as.numeric(annot[[sat_numcrit]])
  qcSeqSat_Scatter <- ggplot(annot,
                             aes(x = .data[[sat_numcrit]],
                                 y = SequencingSaturation,
                                 col = as.factor(.data[[sat_factcrit]]))) +
    geom_point(size = 3, alpha=0.4) +
    geom_hline(yintercept = sat_cutoff, size = 1, color='hotpink', lty="21", alpha=0.5) +
    xlim(min(annot[[sat_numcrit]])*0.8, max(annot[[sat_numcrit]])*1.2) +
    ylim(0,100) +
    ggtitle(paste0('Sequencing Saturation: Dimensionally')) +
    labs(x = sat_numcrit, y = 'Sequencing saturation', color = sat_factcrit)
  if (is.null(preset_theme)) {
      qcSeqSat_Scatter <- qcSeqSat_Scatter + theme_minimal(base_size = 12)
  }
  qcSeqSat_Scatter <- qcSeqSat_Scatter + theme(panel.grid = element_blank(),
                                             axis.line.x = element_line(),
                                             axis.line.y = element_line())
  qcSeqSat_Scatter

  ggsave(paste(qc_dir, paste0('seqQC_seqSat_scatter.', fileType), sep='/'),
         qcSeqSat_Scatter, device=fileType, dpi=300, width=6, height=3, limitsize = FALSE, scale = 1.5)

  # produce graph to investigate distribution of seq saturation by group in a violin plot by group
  qcSeqSat_Viol <- ggplot(annot,
                          aes(x = as.factor(.data[[sat_factcrit]]),
                              y = SequencingSaturation,
                              fill = as.factor(.data[[sat_factcrit]]))) +
    geom_violin(show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.colour=NA, show.legend = FALSE) +
    geom_hline(yintercept = sat_cutoff, size = 1, color='hotpink', lty="21", alpha=0.5) +
    ylim(0,100) +
    ggtitle(paste0('Sequencing Saturation by ', sat_factcrit)) +
    labs(x = sat_factcrit, y = 'Sequencing saturation', fill = '')
  if (is.null(preset_theme)) {
      qcSeqSat_Viol <- qcSeqSat_Viol + theme_minimal(base_size = 12)
  }
  qcSeqSat_Viol <- qcSeqSat_Viol + theme(panel.grid = element_blank(),
                                           axis.line.x = element_line(),
                                           axis.line.y = element_line())
  qcSeqSat_Viol

  ggsave(paste(qc_dir, paste0('seqQC_seqSat_violin.', fileType), sep='/'),
         qcSeqSat_Viol, device=fileType, dpi=300,
         width=max(length(unique(annot[[sat_factcrit]])), 5),
         height=3, units="in", limitsize = FALSE, scale = 2.5)

  # If return_plots = TRUE, return list of plots
  if (return_plots == TRUE) {
    seqqc = list('reads' = qcLine, 'saturation' = qcSeqSat)
    return(seqqc)
    }
}
