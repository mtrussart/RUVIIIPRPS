#' Get a metric or plot from the "metadata" of SummarizedExperiment object.

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
#' @param sub.file.name Symbol. A symbol indicating the sub-file name to which the results of the metric are assigned.
#' @param required.function Symbol. A symbol indicating the name of function that is required to be applied before getting
#' the results of that metric.
#' @param message.to.print Symbol. A symbol to print.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

getMetricFromSeObj <- function(
        se.obj,
        slot = 'Metrics',
        assay.names,
        assessment.type = 'gene.level',
        assessment,
        method = NULL,
        variables,
        file.name,
        sub.file.name = NULL,
        required.function,
        message.to.print,
        verbose = TRUE
        ){
    if(is.null(required.function)){
        required.function <- 'required function'
    }
    # Create an empty list ####
    all.outputs <- list()

    # Check the metadata ####
    if (length(se.obj@metadata) == 0) {
        stop('The current SummarizedExperiment object does not contain "metadat".')
    }
    # Check the slot ####
    if (!slot %in% names(se.obj@metadata)) {
        stop(paste0('The "metadat" of the current SummarizedExperiment object does not contain ', slot ,' slot.'))
    }
    # Check the all ####
    for(x in assay.names){
        if (!x %in% names(se.obj@metadata[[slot]])) {
            stop(paste0('The "', slot, '" of in the "metadata" of the current SummarizedExperiment object does not contain any metrics for "', x ,'" data.'))
        }
        if (!assessment.type %in% names(se.obj@metadata[[slot]][[x]])) {
            stop(paste0('The metrics of for "', x , '" data does not contain any "', assessment.type ,'" assessments.'))
        }
        if (!assessment %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]])) {
            stop(paste0('The ', assessment.type , ' assessments of the "', x, '" data does not contain any ', assessment ,' data. ',
                        'Please run the "', required.function, '" function first.'))
        }
        if (!method %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]])) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data does not contain any data for the "', assessment ,'" method.',
                        'Please check the parameters of the "', required.function, '" function.'))
        }
        if (!variables %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]])) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ',
                        'does not contain any data for the "', variables ,'" variable.',
                        'Please check the parameters of the "', required.function, '" function.')
                 )
        }
        if (!file.name %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]]) ) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ', 'for the "',
                        variables ,'" variable.', 'does not contain any ', file.name, 'file.',
                        'Please check the parameters of the "', required.function, '" function.')
                 )
        }
        if(!is.null(sub.file.name)){
            if (!sub.file.name %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]]) ) {
                stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ', 'for the "',
                            variables ,'" variable.', ' with the ', file.name, 'file ', 'does not contain ', sub.file.name, 'file.',
                            'Please check the parameters of the "', required.function, '" function.')
                )
            }

            all.outputs[[x]] <- se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]][[sub.file.name]]
        } else {
            all.outputs[[x]] <- se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]]
        }
        printColoredMessage(
            message = paste0('- Obtain the ', message.to.print, ' the "', x , '" data.'),
            color = 'blue',
            verbose = verbose
        )
    }
    return(all.outputs)
}
