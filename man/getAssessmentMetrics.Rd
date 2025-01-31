% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAssessmentMetrics.R
\name{getAssessmentMetrics}
\alias{getAssessmentMetrics}
\title{Create all possible assessment plots and metrics for the variables.}
\usage{
getAssessmentMetrics(se.obj, variables, plot.output = TRUE, save.se.obj = TRUE)
}
\arguments{
\item{se.obj}{A summarized experiment object.}

\item{variables}{Symbol. A symbol and a vector of symbols indicating the columns names of variables in the samples
annotation in the SummarizedExperiment object. The 'variables' can be categorical and continuous variables.}

\item{plot.output}{Logical. Whether to print the plot of all possible assessment plots or not. The default is set to 'TRUE'.}

\item{save.se.obj}{Logical. Whether to save the results into the SummarizedExperiment object. The default is TRUE.}
}
\value{
A list of all possible assessment metrics for the variables.
}
\description{
This functions provides the names of all possible assessment plots and metrics for the given variable(s). The list
will be used in the 'assessVariation' and 'assessNormalization' functions.
}
\author{
Ramyar Molania
}
