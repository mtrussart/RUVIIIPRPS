% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findMnn.R
\name{findMnn}
\alias{findMnn}
\title{Find mutual nearest neighbors in RNA-seq data.}
\usage{
findMnn(
  se.obj,
  assay.name,
  uv.variable,
  nb.mnn = 1,
  clustering.method = "kmeans",
  nb.clusters = 3,
  data.input = "expr",
  nb.pcs = 2,
  center = TRUE,
  scale = FALSE,
  svd.bsparam = bsparam(),
  normalization = "CPM",
  apply.cosine.norm = TRUE,
  regress.out.variables = NULL,
  hvg = NULL,
  apply.log = TRUE,
  pseudo.count = 1,
  mnn.bpparam = SerialParam(),
  mnn.nbparam = KmknnParam(),
  assess.se.obj = TRUE,
  remove.na = "both",
  plot.output = TRUE,
  output.name = NULL,
  prps.group = NULL,
  save.se.obj = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A summarized experiment object.}

\item{assay.name}{Symbol. A symbol that indicates the name of the assay in the SummarizedExperiment object. This data
should be the one that will be used as input for RUV-III normalization.}

\item{uv.variable}{Symbol. A symbol indicating the name of the column in the SummarizedExperiment object. The
'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
divided into 'nb.clusters' groups using the 'clustering.method' method.}

\item{clustering.method}{Symbol. A symbol that indicates the choice of clustering method for grouping the 'uv.variable'
if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.}

\item{nb.clusters}{Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
continuous variable. The default is 3.}

\item{normalization}{Symbol. A Symbol that indicates which normalization methods should be applied before finding the
knn. The option are 'CPM', 'TMM', 'VST'. The default is 'CPM'. Refer to the 'applyOtherNormalization' for more details.}

\item{apply.cosine.norm}{Logical. Idicates whether to apply cosine normalization on the data or not. The default is
set to 'FALSE'.}

\item{regress.out.variables}{Symbol. A symbol or vector of symbols indicating the column name(s) that contain
biological variable(s) in the SummarizedExperiment object. These variables will be regressed out from the data before
finding genes that are highly affected by unwanted variation variable. The default is set to 'NULL', indicates the
regression will not be applied.}

\item{hvg}{Vector. A vector of the names of the highly variable genes. These genes will be used to find mutual nearest
neighbors samples across the batches. The default is set to 'NULL'. The 'findBioGenes' function can be used for specify
a set of genes.}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.}

\item{pseudo.count}{Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.}

\item{mnn.bpparam}{Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
. The default is SerialParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.}

\item{mnn.nbparam}{Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
. The default is KmknnParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
function for more details.}

\item{remove.na}{Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
the checkSeObj function for more details.}

\item{plot.output}{Logical. If 'TRUE', the function plots the distribution of MNN across the batches.}

\item{output.name}{Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
based on "paste0(uv.variable, '||' , assay.name)".}

\item{prps.group}{Symbol. A symbol specifying the name of the PRPS group. If is 'NULL', the function will select a name
based on "paste0('prps|mnn|', uv.variable)".}

\item{save.se.obj}{Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
or to output the result as list. By default it is set to TRUE.}

\item{verbose}{Logical. If 'TRUE', shows the messages of different steps of the function.}

\item{mnn}{Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute. The default
is set 1.}
}
\description{
This function finds mutual nearest neighbors between all pairs of batches in RNA-seq data. The mutual nearest neighbors
will be used to find and create pseudo samples and eventually pseudo-replicates for RUV-III normalization.
}
\author{
Ramyar Molania
}
