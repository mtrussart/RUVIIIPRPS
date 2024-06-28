#' Add overall plot to the "metadata" of SummarizedExperiment object.

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param slot Symbol. A symbol indicating the name of the slots in the 'metadata' of the SummarizedExperiment object.
#' @param assessment.type Symbol. A symbol indicating the type of assessment. Options are 'gene.level' or 'global.level'.
#' @param assessment Symbol. A symbol indicating the name of the metric to be checked.
#' @param method Symbol. A symbol indicating the method used to calculate the metric.
#' @param variables Symbol. A symbol or a vector of symbols indicating the variables used to calculate the metric.
#' @param file.name Symbol. A symbol indicating the file name to which the results of the metric are assigned.
#' @param plot.data Symbol. A symbol indicating the name of the plots file for the metric.

addOverallPlotToSeObj <- function(
        se.obj,
        slot = 'Plots' ,
        assessment.type = 'gene.level',
        assessment,
        method,
        variables,
        file.name,
        plot.data
        ){
    if (!slot %in%  names(se.obj@metadata)) {
        se.obj@metadata[[slot]] <- list()
    }
    if (!assessment.type %in%  names(se.obj@metadata[[slot]]) ) {
        se.obj@metadata[[slot]][[assessment.type]] <- list()
    }
    if (!assessment %in%  names(se.obj@metadata[[slot]][[assessment.type]]) ) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]] <- list()
    }
    if (!method %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]] <- list()
    }
    if (!variables %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]] <- list()
    }
    if (!file.name %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- list()
    }
    se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- plot.data
    return(se.obj)
}


