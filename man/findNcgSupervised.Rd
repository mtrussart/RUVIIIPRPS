% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findNcgSupervised.R
\name{findNcgSupervised}
\alias{findNcgSupervised}
\title{Finds a set of negative control genes using supervised approaches.}
\usage{
findNcgSupervised(
  se.obj,
  assay.name,
  bio.variables,
  uv.variables,
  approach = "TwoWayAnova",
  ncg.selection.method = "non.overlap",
  nb.ncg = 0.1,
  top.rank.bio.genes = 0.5,
  top.rank.uv.genes = 0.5,
  bio.percentile = 0.2,
  uv.percentile = 0.2,
  grid.group = "uv",
  grid.direction = "increase",
  grid.nb = 20,
  bio.groups = NULL,
  bio.clustering.method = "kmeans",
  nb.bio.clusters = 3,
  uv.groups = NULL,
  uv.clustering.method = "kmeans",
  nb.uv.clusters = 3,
  normalization = "CPM",
  regress.out.uv.variables = NULL,
  regress.out.bio.variables = NULL,
  apply.log = TRUE,
  pseudo.count = 1,
  anova.method = "aov",
  min.sample.for.aov = 3,
  corr.method = "spearman",
  a = 0.05,
  rho = 0,
  min.sample.for.correlation = 10,
  assess.ncg = TRUE,
  variables.to.assess.ncg = NULL,
  nb.pcs = 5,
  center = TRUE,
  scale = FALSE,
  assess.se.obj = TRUE,
  remove.na = "none",
  save.se.obj = TRUE,
  output.name = NULL,
  ncg.group = NULL,
  plot.output = TRUE,
  use.imf = FALSE,
  save.imf = FALSE,
  imf.name = NULL,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object.}

\item{assay.name}{Character. Indicates the name of the assay in the SummarizedExperiment object. This assay should
be the one that will be used for RUV-III-PRPS normalization. We recommend using raw data.}

\item{bio.variables}{Character. Indicates the column names that contain known biological variable(s) in the
SummarizedExperiment object. These biological variables can be categorical or continuous. Continuous variables will be
divided into 'nb.bio.clusters' groups based on a clustering method selected in the 'bio.clustering.method' argument.
This argument cannot be empty.}

\item{uv.variables}{Character. Indicates the column names that contain unwanted variable(s) in the SummarizedExperiment
object. These unwanted variables can be categorical or continuous. Continuous variables will be divided into
'nb.uv.clusters' groups based on a clustering method selected in the 'uv.clustering.method' argument. This argument
cannot be empty.}

\item{approach}{Character. Indicates which NCGs selection method should be used. The options are 'AnovaCorr.PerBatchPerBiology',
'AnovaCorr.AcrossAllSamples' and 'TwoWayAnova'. The default is set to 'TwoWayAnova'. Refer to details for more information.}

\item{ncg.selection.method}{Character. Indicates how to summarize F-statistics obtained from two-way ANOVA and select a
set of genes as negative control genes. The options are 'prod', 'average', 'sum', 'non.overlap', 'auto', and 'quantile'.
The default is set to 'non.overlap'. For more information, refer to the details of the function.}

\item{nb.ncg}{Numeric. Indicates how many genes should be selected as NCG. The value represents the proportion of the
total genes in the SummarizedExperiment object. The default is set to 0.1.}

\item{top.rank.bio.genes}{Numeric. Indicates the top-ranked genes that are highly affected by the biological variation.
This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'. The default
is set to 0.5.}

\item{top.rank.uv.genes}{Numeric. Indicates the top-ranked genes that are highly affected by the unwanted variation.
This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'. The default
is set to 0.5.}

\item{bio.percentile}{Numeric. The percentile cut-off of F-statistics to select genes that are highly affected by
the biological variation. The default is set to 0.8.}

\item{uv.percentile}{Numeric. The percentile cut-off of F-statistics to select genes that are highly affected by
the unwanted variation. The default is set to 0.8.}

\item{grid.group}{Character. Indicates whether the grid search should be performed on biological ('top.rank.bio.genes'),
unwanted ('top.rank.uv.genes') or both factors, when the 'ncg.selection.method' is set to 'auto'. The options are
'bio', 'uv', or 'both'. If set to 'both', the grid search will be performed on both biological and unwanted factors.
If set to 'bio' or 'uv', the grid search will be performed only on biological or unwanted factors. The default is
set to 'uv'.}

\item{grid.direction}{Character. Indicates whether the grid search should be performed in decreasing or increasing order. The
options are 'increase' or 'decrease'. The default is set to 'decrease'.}

\item{grid.nb}{Numeric. Indicates the number of genes for grid search when the 'ncg.selection.method' is set to 'auto'.
In the 'auto' approach, the grid search increases or decreases the initial values of 'top.rank.bio.genes' or
'top.rank.uv.genes' or 'both' to find ~'nb.ncg' of genes as NCGs. The default is set to 20.}

\item{bio.groups}{Character. A character string or a vector of strings indicating the columns names that contain
biological variables in the SummarizedExperiment object. If specified, the 'bio.groups' will be used for grouping
samples into different homogeneous biological groups. If set to 'NULL', the 'bio.variables' will be used for grouping
samples into different homogeneous biological groups.}

\item{bio.clustering.method}{Character. Indicates which clustering methods should be used to group continuous sources
of biological variation. Refer to the 'createHomogeneousBioGroups' function for more details. The default is set to
'kmeans' clustering.}

\item{nb.bio.clusters}{Numeric. Indicates the number of clusters for each continuous source of biological variation.
By default, it is set to 2. This means individual continuous sources of biological variation will be divided into two
groups.}

\item{uv.groups}{Character. A character string or a vector of strings indicating the columns names that contain biological variables
in the SummarizedExperiment object. If specified, the 'uv.groups' will be used for grouping samples into possible
homogeneous sample groups with respect to unwanted variables. If set to 'NULL', the 'uv.variables' will be used for grouping
samples.}

\item{uv.clustering.method}{Character. Indicates which clustering methods should be used to group continuous sources
of unwanted variation. Refer to the 'createHomogeneousUvGroups' function for more details. The default is set to
'kmeans' clustering.}

\item{nb.uv.clusters}{Numeric. Indicates the number of clusters for each continuous source of unwanted variation.
By default, it is set to 2. This means individual continuous sources of unwanted variation will be divided into two
groups.}

\item{normalization}{Character. Indicates which normalization method should be applied to the data before finding genes
that are affected by biological variation. The default is set to 'CPM'. If set to 'NULL', no normalization will be applied.
Refer to the 'applyOtherNormalizations' function for more details.}

\item{regress.out.uv.variables}{Character. Indicates the column names that contain unwanted variation variables in the
SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
affected by biological variation. The default is 'NULL', indicating regression will not be applied.}

\item{regress.out.bio.variables}{Character. Indicates the column names that contain biological variables in the
SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
affected by unwanted variation. The default is 'NULL', indicating regression will not be applied.}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
analysis. The default is set to 'TRUE'.}

\item{pseudo.count}{Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
default is 1.}

\item{anova.method}{Character. Indicates which ANOVA method should be used to compute association between gene-level
expression and a continuous variable.}

\item{min.sample.for.aov}{Numeric. Indicates the minimum number of samples to be present in each group before applying
the ANOVA. The default is 3.}

\item{corr.method}{Character. Indicates which correlation method should be used to compute association between gene-level
expression and a continuous variable. The default is 'spearman'.}

\item{a}{Numeric. The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.}

\item{rho}{Numeric. The value of the hypothesized correlation to be used in the hypothesis testing, by default it is set to 0.}

\item{min.sample.for.correlation}{Numeric. Indicates the minimum number of samples to be considered in each group before
applying the correlation analysis. The default is 10.}

\item{assess.ncg}{Logical. Indicates whether to assess the performance of selected NCG or not. This analysis involves
principal component analysis on only the selected NCG and then exploring the R² or vector correlation between the 'nb.pcs'
first principal components and the specified variables.}

\item{variables.to.assess.ncg}{Character. Indicates the column names of the SummarizedExperiment object that contain
variables whose association with the selected genes as NCG needs to be evaluated. The default is 'NULL'. This means all
the variables specified in the 'bio.variables' and 'uv.variables' will be assessed.}

\item{nb.pcs}{Numeric. Indicates the number of the first principal components of selected NCG to be used to assess
the performance of NCGs. The default is 5.}

\item{center}{Logical. Indicates whether to center the data before applying principal component analysis or not.
Refer to the 'computePCA' function for more details. The default is set to 'TRUE'.}

\item{scale}{Logical. Indicates whether to scale the data before applying principal component analysis.
Refer to the 'computePCA' function for more details. The default is 'FALSE'.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
'checkSeObj' will be applied. The default is set to 'TRUE'.}

\item{remove.na}{Character. Indicates whether to remove NA or missing values from either the 'assays', 'sample.annotation',
'both' or 'none'. If 'assays' is selected, the genes that contain NA or missing values will be excluded. If
'sample.annotation' is selected, the samples that contain NA or missing values for any 'bio.variables' or
'uv.variables' will be excluded. The default is set to 'none'.}

\item{save.se.obj}{Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
to output the result as a logical vector. The default is set to 'TRUE'. The file will be saved in
'se.obj@metadata$NCG$supervised$output.name".}

\item{output.name}{Character. A representation for the output's name. If set to 'NULL', the function will choose a name
automatically. In this case, the file name will be constructed as paste0(sum(ncg.selected),'|', paste0(bio.variables,
collapse = '&'), '|', paste0(uv.variables, collapse = '&'),'|TWAnova:', ncg.selection.method, '|', assay.name).}

\item{ncg.group}{Character. A symbol indicating the name of the group of NCG.}

\item{plot.output}{Character. Indicates whether to plot the output or not.}

\item{use.imf}{Logical. Indicates whether to use the intermediate file or not. The default is set to 'FALSE'.}

\item{save.imf}{Logical. Indicates whether to save the intermediate file in the SummarizedExperiment object or not.
If set to 'TRUE', the function saves the results of the two-way ANOVA. Subsequently, if users wish to adjust parameters
such as 'nb.ncg', 'ncg.selection.method', 'top.rank.bio.genes', and 'top.rank.uv.genes', the two-way ANOVA will not
be recalculated. This accelerates parameter tuning for NCG selection. The default value is 'FALSE'.}

\item{imf.name}{Character. Indicates the name to use when saving the intermediate file. If set to 'NULL', the function
will create a name. In this case, the file name will be constructed as
paste0(assay.name, '|TwoWayAnova|', ncg.selection.method). The default is 'NULL'.}

\item{verbose}{Logical. If 'TRUE', shows the messages of different steps of the function.}
}
\value{
Either the SummarizedExperiment object containing a set of negative control genes in the metadata or a
logical vector of the selected negative control genes.
}
\description{
This function contains three different functions including 'findNcgByTwoWayAnova', 'findNcgAcrossSamples' and
'findNcgPerBiologyPerBatch' to find a set of genes as negative control genes (NCG) for RUV-III-PRPS normalization. We
refer to each function for more details.
}
\author{
Ramyar Molania
}
