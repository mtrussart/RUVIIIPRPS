% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeGenesVariableCorrelation.R
\name{computeGenesVariableCorrelation}
\alias{computeGenesVariableCorrelation}
\title{Computes the correlation between individual gene expression and a continuous variable.}
\usage{
computeGenesVariableCorrelation(
  se.obj,
  assay.names = "all",
  variable,
  method = "spearman",
  a = 0.05,
  rho = 0,
  apply.log = TRUE,
  pseudo.count = 1,
  plot.top.genes = FALSE,
  nb.top.genes = 3,
  apply.round = TRUE,
  assess.se.obj = TRUE,
  remove.na = "both",
  save.se.obj = TRUE,
  override.check = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object.}

\item{assay.names}{Character. A character string or vector of character strings specifying the name(s) of the assay(s)
in the SummarizedExperiment object for computing the correlation. By default is set to 'all', which indicates all
assays in the SummarizedExperiment object will be selected.}

\item{variable}{Character. A character string indicating a column name in the SummarizedExperiment object that contains
a continuous variable, such as library size, tumor purity, etc.}

\item{method}{Character. Specifies which correlation method should be used. The options are 'pearson', 'kendall', or
"spearman". The default is set to 'spearman'.}

\item{a}{Numeric. The significance level used for the confidence intervals in the correlation. The default is 0.05.
Refer to the 'correls' function from the Rfast R package for more details.}

\item{rho}{Numeric. The hypothesized correlation value to be used in hypothesis testing. The default is 0.
Refer to the 'correls' function from the Rfast R package for more details.}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data before computing the correlation.
By default, the log transformation is applied.}

\item{pseudo.count}{Numeric. A numeric value to be added as a pseudo count to all measurements before log transformation.
The default is set to 1.}

\item{plot.top.genes}{Logical. Indicates whether to plot the gene expression for the top or bottom genes based on
correlation values. The default is set to 'FALSE'.}

\item{nb.top.genes}{Numeric. Defines the number of genes with the highest or lowest correlation to the variable to plot.
The default is 3.}

\item{apply.round}{Logical. Indicates whether to round the correlation coefficients. The default is TRUE.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment object.
Refer to the 'checkSeObj' function for more details. The default is TRUE.}

\item{remove.na}{Character. Specifies whether to remove NA or missing values from the assays and variable.}

\item{save.se.obj}{Logical. Indicates whether to save the result in the metadata of the current SummarizedExperiment
object or to output the result. By default, it is set to TRUE.}

\item{override.check}{Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
determine if the correlation has already been computed with the current parameters. If so, the metric will not be recalculated.
The default is FALSE.}

\item{verbose}{Logical. If 'TRUE', displays the messages for different steps of the function.}
}
\value{
Either a SummarizedExperiment object or a list containing the correlation coefficients for compuated for individual
genes.
}
\description{
This function computes Spearman or Pearson correlations between individual gene-level expression of each assay and
a continuous variable in SummarizedExperiment object.
}
\references{
Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
Nature Biotechnology, 2023
}
\author{
Ramyar Molania
}
