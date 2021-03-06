% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_ora.R
\name{run_ora}
\alias{run_ora}
\title{run_ora

Pathway Analysis helper function to run Over Representation Analysis}
\usage{
run_ora(
  geneList = geneList,
  pathways = pathways,
  test = test,
  path_fc = path_fc,
  path_pval = path_pval,
  outdir = outdir,
  sp = species
)
}
\arguments{
\item{geneList}{a list of ranked genes, the way that they are ranked is determined in run_pathways.R}

\item{pathways}{are the genesets that we are working with, default is Reactome}

\item{test}{is the test or contrast that is being compared for pathway analysis, from DE results}

\item{path_fc}{the fold change cutoff to use for gene list}

\item{path_pval}{the pvalue to threshold the enrichment results by}

\item{outdir}{the output directory for plots}
}
\value{
ora The results from over representation analysis
}
\description{
run_ora

Pathway Analysis helper function to run Over Representation Analysis
}
\examples{
 run_ora(geneList = geneList, pathways = pathways, test = test, path_fc = path_fc, path_pval = path_pval, outdir = outdir)

}
