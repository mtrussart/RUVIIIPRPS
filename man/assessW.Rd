% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assessW.R
\name{assessW}
\alias{assessW}
\title{Assess the W matrix of RUV-III}
\usage{
assessW(
  se.obj,
  variables,
  bio.variables = NULL,
  uv.variables = NULL,
  compare.w = FALSE,
  plot.output = TRUE,
  save.se.obj = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A summarized experiment object.}

\item{variables}{Symbol. A symbol or symbols representing the label of variable(s), such as cancer subtypes, tumour
purity, librar size ... within the SummarizedExperiment object. This can comprise a vector containing either categorical,
continuous, or a combination of both variables.}

\item{bio.variables}{Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
subtypes, tumour purity, ... within the SummarizedExperiment object. This can comprise a vector containing either
categorical, continuous, or a combination of both variables.}

\item{uv.variables}{Symbol. A symbol or symbols representing the label of unwanted variable(s), such as cancer
batch effects, library size, ... within the SummarizedExperiment object. This can comprise a vector containing either
categorical, continuous, or a combination of both variables.}

\item{compare.w}{Logical.  Specifies whether to compare different W matrices in their ability to capture unwanted
variation and lack correlation with biological variables. See the function's details for more information.}

\item{plot.output}{Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.}

\item{save.se.obj}{Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.}

\item{verbose}{Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.}
}
\description{
This functions assesses the correlation between the W of RUV-III normalized data with known variables. Tt calculates
the association of the columns of the W matrix with both specified biological and unwanted variables. The correlation
values obtained for the biological variables will be subtracted from 1. Finally, the average of all these correlations
is computed to determine the final correlation values. A higher correlation indicates better performance.
}
\author{
Ramyar Molania
}
