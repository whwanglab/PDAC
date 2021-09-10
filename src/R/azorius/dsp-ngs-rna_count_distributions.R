# 

##############################
#### Variables to Assign #####
##############################

#
# BEFORE DOING ANYTHING: Set your working directory to where your azorius output files are
# https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces 
#


# Summarized Count File Output by Azorius
sumcountfile <- 'Broad_CRC_summarized_counts.txt'

# Desired Count Column
# One of "Count_Val", "Neg_Norm", "HK_Norm"
countcol <- "HK_Norm"

# X Axis Log transformed? TRUE or FALSE
logX <- T

# Label for X Axis
xaxislab <- 'HK Normalized Counts'

# Column with values you would like to compare detection frequency across
comparison_column <- 'timepoint'

# Order you would like the x-axis to appear in the charts
# Leave as empty vector if you don't know / care <- c()
comparison_sort_order <- c('B', 'R', 'D')

# Plot labels
fill_label <- 'Time Point'
expname <- "Stanford RNA Batch #1 Log10"

fillcolors <- c('#A6CE38', '#194983', '#56565B')

#
#
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
# If you have a feature request or encounter a bug, let Zach know.
#
#
#
#
#
#### Required packages ####
required_packages <- c('ggplot2', 'scales')

for (pkg in required_packages) {
  if (!pkg %in% installed.packages()) install.packages(pkg) 
}

library(ggplot2)
library(scales)

#### Read in and Format Data ####
# Summary Count file
sumdf <- read.delim(sumcountfile)
colnames(sumdf)[which(colnames(sumdf) == comparison_column)] <- 'CompCol'
if (length(comparison_sort_order) != 0) {
  sumdf$CompCol <- factor(sumdf$CompCol, levels = comparison_sort_order)
}
colnames(sumdf)[which(colnames(sumdf) == countcol)] <- 'CountCol'

xBreaks <- c(1:10 %o% 10^(0:7))
xLabs <- formatC(xBreaks, format = 'd', big.mark = ',')
xLabs[-grep(1, xLabs)] <- ''

pdf(paste0(gsub(' ', '_', expname), '_count_distributions.pdf'),
    height=6, width = 7)
for (gene in unique(sumdf$Gene)) {
  genedf <- sumdf[sumdf$Gene == gene,]
  genedf$Count_Val[genedf$Count_Val == 0] <- 1
  
  geneplot <- ggplot(genedf, aes(x=CountCol, fill=CompCol)) +
    geom_density(alpha=0.6, size=1) +
    scale_fill_manual(values = fillcolors) +
    labs(x=xaxislab, y = 'Density',
         title = expname,
         subtitle = gene) +
    theme_classic() +
    theme(axis.text = element_text(color='black', size=16),
          axis.title = element_text(size=20),
          plot.title = element_text(size=22),
          plot.subtitle = element_text(size=16),
          legend.title = element_blank(),
          legend.position = 'bottom')
  
  if (logX) {
    geneplot <- geneplot + scale_x_log10(breaks=xBreaks, labels=xLabs)
  }
  
  print(geneplot)
}
dev.off()
