% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_PCA.R
\name{run_PCA}
\alias{run_PCA}
\title{run_PCA

Runs principle component analysis (PCA) on dsp data}
\usage{
run_PCA(
  df,
  outdir,
  color = "cluster",
  symbol = "dsp_slide",
  size = TRUE,
  colors = colors
)
}
\arguments{
\item{df}{read in dsp excel workbook dataset}

\item{outdir}{output directory for figures}

\item{color}{category to color the graph by}

\item{symbol}{category to change symbols in the graph}

\item{size}{should the size of the points change by gene_count}

\item{outlier_cutoff}{multiplier for outlier cutoff}
}
\value{
clustering information
}
\description{
run_PCA

Runs principle component analysis (PCA) on dsp data
}
\examples{
pca <- run_PCA(df, "file/path", outlier_cutoff)

}
