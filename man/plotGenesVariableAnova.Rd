% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGeneVariableAnova.R
\name{plotGenesVariableAnova}
\alias{plotGenesVariableAnova}
\title{Plots F-statistics obtained from ANOVA.}
\usage{
plotGenesVariableAnova(
  se.obj,
  assay.names = "all",
  variable,
  anova.method = "aov",
  plot.ncol = 4,
  plot.nrow = 1,
  plot.output = TRUE,
  save.se.obj = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object.}

\item{assay.names}{Character or vector. A character string or a vector of strings for selecting the name(s) of the
assay(s) in the SummarizedExperiment object. By default, all assays in the SummarizedExperiment object will be selected,
If the specifed ANOVA has already been computed for them using 'computeGenesVariableANOVA' function.
The ANOVA results for each specified assay will be plotted.}

\item{variable}{Character. Indicates the column name in the SummarizedExperiment object that contains a categorical
variable such as sample types or batches. The association between individual gene expression and the specified levels
must have already been computed by the 'computeGenesVariableANOVA' function.}

\item{anova.method}{Logical. Specifies which ANOVA method results should be used for plotting.}

\item{plot.ncol}{Numeric. A numeric value specifying the number of columns in the plot grid. The default is set to 4.}

\item{plot.nrow}{Numeric. A numeric value specifying the number of rows in the plot grid. The default is set to 1.}

\item{plot.output}{Logical. Indicates whether to display the boxplot(s) of the ANOVA results. By default, this
is set to FALSE.}

\item{save.se.obj}{Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object
or to output the result as a list. By default, this is set to TRUE.}

\item{verbose}{Logical. If TRUE, displays process messages during execution.}
}
\value{
Either a SummarizedExperiment object or a list containing the plots of ANOVA results for the specified
categorical variable for each assay.
}
\description{
This functions plots the results of the ANOVA analysis computed by the 'computeGenesVariableANOVA' function.
}
\author{
Ramyar Molania
}
