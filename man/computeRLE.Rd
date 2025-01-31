% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeRLE.R
\name{computeRLE}
\alias{computeRLE}
\title{Calculates relative log expression (RLE) of RNA-seq data.}
\usage{
computeRLE(
  se.obj,
  assay.names = "all",
  apply.log = TRUE,
  pseudo.count = 1,
  outputs.to.return = "all",
  assess.se.obj = TRUE,
  remove.na = "assays",
  save.se.obj = TRUE,
  override.check = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object.}

\item{assay.names}{Character. A character string or a vector of character strings specifying the name(s) of the assay(s)
in the SummarizedExperiment object to calculate RLE data, medians, and interquartile ranges. The default is set to "all",
which indicates that all assays in the SummarizedExperiment object will be selected.}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data. The default is set to 'TRUE'.
Please note that any RNA-seq data (assays) must be in log scale before computing RLE.}

\item{pseudo.count}{Numeric. A value to be added as a pseudo count to all measurements of the assay(s) before applying
log transformation. This avoids -Inf for measurements equal to 0. The default is set to 1.}

\item{outputs.to.return}{Character. Specifies the type of RLE computations to be performed and the data to be returned.
Options include "all", "rle", "rle.med", "rle.iqr", and "rle.med.iqr". Selecting "all" returns RLE data along with
medians and interquartile ranges. Choosing "rle" returns only the RLE data for each assay. "rle.med" returns only the
RLE  medians. "rle.iqr" returns only the interquartile ranges of the RLE data. "rle.med.iqr" returns both the
RLE medians and interquartile ranges. The default is set to  'all'.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment object. If 'TRUE', the 'checkSeObj'
function will be applied. The default is set to 'TRUE'}

\item{remove.na}{Character. Indicates whether to remove NA or missing values from the assays. Options are 'assays' or
'none'. The default is "assays", meaning all NA or missing values in the assays will be removed before computing RLE.
Refer to the 'checkSeObj' function for more details.}

\item{save.se.obj}{Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object,
or output the result as a list. The default is set to 'TRUE'.}

\item{override.check}{Logical. When set to 'TRUE', the function verifies if the RLE has already been computed for the
current parameters in the SummarizedExperiment object. If already computed, the metric will not be recalculated. The
default is set to 'FALSE'.}

\item{verbose}{Logical. If 'TRUE', shows messages for each step of the function.}
}
\value{
Either a SummarizedExperiment object that contain the RLE results or a list containing the RLE data for each
individual assay in the SummarizedExperiment object.
}
\description{
This function calculates relative log expression (RLE) of the assay(s) in a SummarizedExperiment object. In addition,
the function returns the RLE medians and interquartile ranges (IQRs) of each sample for individual assay(s).
}
\details{
RLE plots are used to reveal trends, temporal clustering and other non-random patterns resulting from unwanted variation
in gene expression data. To generate RLE plots, we first form the log ratio log(yig/yg) of the raw count yig for
gene g in the sample labeled i relative to the median value yg of the counts for gene g taken across all samples. We
then generate a box plot from all the log ratios for sample i and plotted all such box plots along a line, where i
varies in a meaningful order, usually sample processing date. An ideal RLE plot should have its medians centered around
zero, and its box widths and their interquartile ranges (IQRs) should be similar in magnitude. Because of their
sensitivity to unwanted variation, we also examine the relationships between RLE medians and interquartile ranges with
potential sources of unwanted variation and individual gene expression levels in the datasets. In the absence of any
influence of unwanted variation in the data, we should see no such associations.
}
\references{
Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.

Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
2019.

Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
Nature Biotechnology, 2023
}
\author{
Ramyar Molania
}
