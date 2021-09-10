# dsp-ngs-rna_detection_comparisons.R

##############################
#### Variables to Assign #####
##############################

#
# BEFORE DOING ANYTHING: Set your working directory to where your azorius output files are
# https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces 
#

# Significant Tallies File Output by Azorius
sigcountfile <- 'Pancreatic_HiSeq2_significant_tallies.txt'

# Summarized Count File Output by Azorius
sumcountfile <- 'Pancreatic_HiSeq2_summarized_counts.txt'

# Number of genes in panel used (for calculating %detected)
target_count <- 1412

# Column with values you would like to compare detection frequency across
comparison_column <- 'tissue'

# Order you would like the x-axis to appear in the charts
# Leave as empty vector if you don't know / care <- c()
comparison_sort_order <- c()

# Negative Probe 'Gene' Names
negnames <- c('NegProbe', '50merNeg')

# Plot labels
xaxisname <- 'Tissue'
expname <- "Broad Pancreatic"

# USE LOD For significant 
lod <- T

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
required_packages <- c('ggplot2', 'Rmisc', 'scales')

for (pkg in required_packages) {
  if (!pkg %in% installed.packages()) install.packages(pkg) 
}

library(ggplot2)
library(Rmisc)
library(scales)

#### Read in and Format Data ####
# Significant File
df <- read.delim(sigcountfile)
colnames(df)[which(colnames(df) == comparison_column)] <- 'CompCol'

# Summary Count file
sumdf <- read.delim(sumcountfile)
colnames(sumdf)[which(colnames(sumdf) == comparison_column)] <- 'CompCol'
sumdf$NegOrNot <- sumdf$Gene %in% negnames
if (length(comparison_sort_order) != 0) {
  sumdf$CompCol <- factor(sumdf$CompCol, levels = comparison_sort_order)
}

# Formatting significant stuff
if (lod) {
  lod_thresh <- which(sumdf$Gene == 'NegProbe')
}


df$Percent_Detected <- df$Significant_Count / target_count
perplotdf <- summarySE(df, 'Percent_Detected', groupvars = c('CompCol'), conf.interval = 0.95)
if (length(comparison_sort_order) != 0) {
  perplotdf$CompCol <- factor(perplotdf$CompCol, levels = comparison_sort_order)
}


#### Bar Plot ####
ymax <- ceiling(max(perplotdf$Percent_Detected + perplotdf$se) * 10) / 10

pdf(paste0(gsub(' ', '_', expname), '_sig_over_back.pdf'),
    height=5, width = max(8, nrow(perplotdf)/ 0.4))
print(
ggplot(perplotdf, aes(x = CompCol, y=Percent_Detected)) +
  geom_bar(stat='identity', fill='#A6CE38', color='black',
           size=1, width=0.6) +
  geom_errorbar(aes(ymin=Percent_Detected - se,
                    ymax=Percent_Detected + se),
                width=0.4, position = position_dodge(0.9),
                size=1) +
  scale_y_continuous(labels=percent, limits = c(0,ymax),
                     expand = c(0,0), breaks = seq(0,ymax,0.1)) +
  labs(x=xaxisname, y= "Percent Above Background",
       title=expname,
       subtitle='Percent of Significantly Expressed Genes') +
  theme_classic() +
  theme(axis.text = element_text(color='black', size=16),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=16))
)
dev.off()

#### Box Plot ####
ymax <- 10^(nchar(as.character(ceiling(max(sumdf$Count_Val)))))
yBreaks <- c(1:10 %o% 10^(0:7))
yLabs <- formatC(yBreaks, format = 'd', big.mark = ',')
yLabs[-grep(1, yLabs)] <- ''

pdf(paste0(gsub(' ', '_', expname), '_exp_vs_neg_boxplot.pdf'),
    height=6, width = 20 / nrow(perplotdf))
print(
ggplot(sumdf, aes(x = CompCol, y = Count_Val)) +
  geom_boxplot(aes(fill=NegOrNot)) +
  scale_y_log10(breaks=yBreaks, labels=yLabs) +
  scale_fill_manual(values=c('#A6CE38', '#194983'),
                    labels=c('Targets', 'Negatives')) +
  labs(x=xaxisname, y="log10(Reported Count)",
       title=expname,
       subtitle='Gene Counts Compared to Negative Counts') +
  theme_classic() +
  theme(axis.text = element_text(color='black', size=16),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=16),
        legend.title = element_blank(),
        legend.position = 'bottom')
)
dev.off()