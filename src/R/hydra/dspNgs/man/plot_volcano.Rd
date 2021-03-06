% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_volcano.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Volcano Plot}
\usage{
plot_volcano(
  de_results = NULL,
  base_level = NULL,
  target_group = "All Probes",
  fav_targets = NULL,
  targets = NULL,
  outfile = "VolcanoePlot-AllGenes",
  point_color = "RdPu",
  top_method = "Significance",
  n_targets = 15,
  plt_title = "Volcano Plot for All Genes",
  target_ID = "gene",
  FC_ID = "FC",
  Pval_ID = "Pval",
  color_list = NULL,
  save_plot = TRUE
)
}
\arguments{
\item{de_results}{output of dsp_de_analysis function with list of up and downregulated targets, e.g. genes or gene sets}

\item{base_level}{baseline level for DE analysis, used for plotting x-axis label}

\item{target_group}{name of target groups, or , used for adding color contrast to graph}

\item{targets}{a list of targets, e.g. genes, within the target_group, used for subsetting graph}

\item{outfile}{name of file for output}

\item{point_color}{a specific color to be used for plotting points, should be updated to a gradient method}

\item{top_method}{column within the de_results that should be used to select the top genes, selecting highest values from the column. Use Significance rather than Pval for this purpose}

\item{n_targets}{number of genes to show in the cloud}

\item{target_ID}{column within the de_results data frame that contains the target names to be used for labeling graphs}

\item{FC_ID}{column within the de_results data frame that contains the FC x-axis values to be used for plotting}

\item{Pval_ID}{column within the de_results data frame that contains the P-value y-axis values to be used for plotting}

\item{color_list}{the targetset color values if multiple targetsets are being used, can be NULL if only one targetset passed}

\item{save_plot}{boolean of whether to save the plot as a file or just output to the console, defaults to TRUE and uses config information for file format}
}
\value{
image file of volcano plot

plotly interactive volcano plot
}
\description{
Takes results of differential expression analysis and generates a volcano plot
}
\examples{
 plot_volcano(de_results, base_level)

}
