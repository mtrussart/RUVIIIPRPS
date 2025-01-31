% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createHomogeneousUVGroups.R
\name{createHomogeneousUVGroups}
\alias{createHomogeneousUVGroups}
\title{Create all possible homogeneous groups with respect to unwanted variables.--Finalized}
\usage{
createHomogeneousUVGroups(
  se.obj,
  uv.variables,
  clustering.method = "kmeans",
  nb.clusters = 3,
  assess.se.obj = TRUE,
  remove.na = "sample.annotation",
  save.se.obj = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A 'SummarizedExperiment' object.}

\item{uv.variables}{Character. A character string or a vector of character strings specifying the column names of unwanted
variables in the sample annotation of the 'SummarizedExperiment  object. These 'uv.variables' can be either categorical
or continuous variables.}

\item{clustering.method}{Character. A character string specifying the clustering method to be applied for grouping each
continuous biological variable. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans' clustering.}

\item{nb.clusters}{Numeric. A value indicating the number of groups for continuous sources of unwanted variation.
The default is 3. This implies that each continuous variable will be split into 3 groups using the specified
'clustering.method'.}

\item{assess.se.obj}{Logical. Whether to assess the 'SummarizedExperiment' object or not. If 'TRUE', the function
'checkSeObj' will be applied. The default is set to ' TRUE' .}

\item{remove.na}{Character. Indicates whether to remove missing values from the specified variables. The options are
'sample.annotation' or 'none'. The default is set to 'sample.annotation', meaning that missing values in the variables
will be removed.}

\item{save.se.obj}{Logical. Indicates whether to save the results to the metadata of the ' SummarizedExperiment'  object
or not. If ' TRUE' , all the possible homogeneous groups will be saved into 'se.obj->metadata->HomogeneousGroups->
BiologicalVariables ; otherwise, the results will be returned as a vector. The default is set to ' TRUE' .}

\item{verbose}{Logical. If ' TRUE' , displays messages for different steps of the function.}
}
\value{
Either a ' SummarizedExperiment' object containing all possible homogeneous groups in 'metadata->HomogeneousGroups
->UnwantedVariables' or a vector of all possible homogeneous sample groups.
}
\description{
This function generates all possible homogeneous sample groups based on the specified unwanted variation variables.
}
\details{
The function generates all possible homogeneous sample groups based on the specified biological variables. If continuous
variables are provided, the function splits each into a number of clusters determined by ’nb.clusters’, using the
clustering method specified in ’clustering.method’. Ultimately, all combinations of all clusters are created and
each such combination is regarded as a homogeneous sample group concerning biological variables.
}
\author{
Ramyar Molania
}
