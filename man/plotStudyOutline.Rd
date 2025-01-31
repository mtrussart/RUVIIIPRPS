% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotStudyOutline.R
\name{plotStudyOutline}
\alias{plotStudyOutline}
\title{Generates study outline of a SummarizedExperiment object.}
\usage{
plotStudyOutline(
  se.obj,
  variables,
  variable.to.sort = NULL,
  assess.se.obj = TRUE,
  remove.na = "none",
  plot.output = TRUE,
  legend.font.size = 14,
  legend.ncol = 4,
  legend.direction = "horizontal",
  heatmap.legend.side = "bottom",
  column.names.rot = 25,
  save.se.obj = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object.}

\item{variables}{Character or character vector. The label(s) of variable(s) within the SummarizedExperiment
object. This can include a vector containing categorical, continuous, or a combination of both types of variables. This
cannot be empty or NULL.}

\item{variable.to.sort}{Character. The label of the variable used to sort the sample information. The default is 'NULL'.
This means the samples will be plotted as they are.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment object. See the 'checkSeObj()'
function for more details. The default is set to TRUE.}

\item{remove.na}{Character. Specifies whether to remove NA or missing values from the assays. The options are 'assays'
, 'sample.annotation', 'both' or 'none'. The default is set to "none". See the 'checkSeObj' function for more details.}

\item{plot.output}{Logical. Determines whether to plot the study outline. The default is set to 'TRUE'.}

\item{legend.font.size}{Numeric. A numeric value indicating the size of the font of the legend in the plot. The default
is set to 14.}

\item{legend.ncol}{Numeric. A numeric value indicating the number of columns in the legend of the heatmap. The default
is set to 4.}

\item{legend.direction}{Character. A string specifying the direction of the legend in the heatmap. The options are
'horizontal' or 'vertical'. The default is set to 'horizontal'.}

\item{heatmap.legend.side}{Character. A string indicating the side of the legend in the heatmap. The options are
'right', 'bottom', 'left', 'top'. The default is set to 'bottom'.}

\item{column.names.rot}{Numeric. A numeric value indicating the angle of the column labels in the heatmap. The default
is set to 25.}

\item{save.se.obj}{Logical. Indicates whether to save the study outline plot in the metadata of the SummarizedExperiment
object or output the result as a plot. The default is set to 'TRUE'. The plot will be saved in 'se.obj->metadata->StudyOutline'.}

\item{verbose}{Logical. If 'TRUE', shows the messages for different steps of the function.}
}
\value{
Either a SummarizedExperiment object containing a study outline plot or just the plot.
}
\description{
This function generates a heatmap of sample-level features such as batches, biological populations, library size, and
tumor purity. This plot is helpful for exploring how these factors are distributed across samples and for examining
visible unwanted variation in the data. Further, it can visually reveal confounder factors in the data. It can be initially
generated based on the known variables and updated with any estimated biological and unwanted variables during the
normalization process. We highly recommend generating the plot before performing downstream analyses.
}
\author{
Ramyar Molania
}
