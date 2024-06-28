#' Add results and plots of a metric to SummarizedExperiment object.

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param slot Symbol. A symbol indicating the name of the slots in the 'metadata' of the SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols specifying the name(s) of the assay(s) in the SummarizedExperiment
#' object. The default is set to 'all'.
#' @param assessment.type Symbol. A symbol indicating the type of assessment. Options are 'gene.level' or 'global.level'.
#' @param assessment Symbol. A symbol indicating the name of the metric to be checked.
#' @param method Symbol. A symbol indicating the method used to calculate the metric.
#' @param variables Symbol. A symbol or a vector of symbols indicating the variables used to calculate the metric.
#' @param file.name Symbol. A symbol indicating the file name to which the results of the metric are assigned.
#' @param results.data Symbol. A symbol indicating the name of the results file for the metric.

addMetricToSeObj <- function(
        se.obj,
        slot = 'Metrics',
        assay.names,
        assessment.type = 'gene.level',
        assessment,
        method,
        variables,
        file.name,
        results.data
        ){
    variables <- paste0(variables, collapse = '&')
    for(x in assay.names){
        if (length(se.obj@metadata) == 0) {
            se.obj@metadata[[slot]] <- list()
        }
        # Check the slot
        if (!slot %in% names(se.obj@metadata)) {
            se.obj@metadata[[slot]] <- list()
        }
        if (!x %in% names(se.obj@metadata[[slot]])) {
            se.obj@metadata[[slot]][[x]] <- list()
        }
        if (!assessment.type %in% names(se.obj@metadata[[slot]][[x]])) {
            se.obj@metadata[[slot]][[x]][[assessment.type]] <- list()
        }
        if (!assessment %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]])) {
            se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]] <- list()
        }
        if (!method %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]])) {
            se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]] <- list()
        }
        if (!variables %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]])) {
            se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]] <- list()
        }
        if (!file.name %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]]) ) {
            se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- list()
        }
        se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- results.data[[x]]
    }
    return(se.obj)
}
