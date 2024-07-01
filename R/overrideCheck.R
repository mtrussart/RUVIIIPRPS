#' Check whether to override a current metric or not.

#' @author Ramyar Molania

#' @description
#' This function checks whether to re-calculate the specified metric or not.

#' @param se.obj A SummarizedExperiment object.
#' @param slot Symbol. A symbol indicating the name of the slots in the 'metadata' of the SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols specifying the name(s) of the assay(s) in the SummarizedExperiment
#' object. The default is set to 'all'.
#' @param assessment.type Symbol. A symbol indicating the type of assessment. Options are 'gene.level' or 'global.level'.
#' @param assessment Symbol. A symbol indicating the name of the metric to be checked.
#' @param method Symbol. A symbol indicating the method used to calculate the metric.
#' @param variable Symbol. A symbol indicating the variable used to calculate the metric.
#' @param file.name Symbol. A symbol indicating the file name to which the results of the metric are assigned.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assays
#' @export

overrideCheck <- function(
        se.obj,
        slot = 'Metrics',
        assay.names,
        assessment.type,
        assessment,
        method,
        variable,
        file.name,
        verbose = TRUE
        ){
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    printColoredMessage(
        message = paste('-- Check to override the current', assessment, ' analysis or not:'),
        color = 'magenta',
        verbose = verbose
    )
    compute.metric <- FALSE
    if (length(se.obj@metadata) == 0 | !'Metrics' %in% names(se.obj@metadata)) {
        printColoredMessage(
            message = paste0('- The ',  assessment, ' data cannot be found in the SummarizedExperiment, the this will be computed for all the assay(s).'),
            color = 'blue',
            verbose = verbose
        )
        compute.metric <- TRUE
        selected.assays <- assay.names
    } else {
        selected.assays <- c()
        for (x in levels(assay.names)) {
            if (!x %in% names(se.obj@metadata[['Metrics']])) {
                selected.assays <- c(selected.assays, x)
            } else if (!assessment.type %in% names(se.obj@metadata[['Metrics']][[x]])){
                selected.assays <- c(selected.assays, x)
            } else if (!assessment %in% names(se.obj@metadata[['Metrics']][[x]][[assessment.type]])) {
                selected.assays <- c(selected.assays, x)
            } else if (!method %in% names(se.obj@metadata[['Metrics']][[x]][[assessment.type]][[assessment]])) {
                selected.assays <- c(selected.assays, x)
            } else if (!variable %in% names(se.obj@metadata[['Metrics']][[x]][[assessment.type]][[assessment]][[method]])) {
                selected.assays <- c(selected.assays, x)
            } else if (!file.name %in% names(se.obj@metadata[['Metrics']][[x]][[assessment.type]][[assessment]][[method]][[variable]])) {
                selected.assays <- c(selected.assays, x)
            }
        }
        if(length(selected.assays) == 0){
            printColoredMessage(
                message = paste0(
                    '- The ', assessment , ' data have been already computed for the "',
                    variable, '" variable for all the assay(s).'),
                color = 'blue',
                verbose = verbose
            )
            compute.metric <- FALSE
        } else if (length(selected.assays) == length(assay.names)){
            printColoredMessage(
                message = paste0(
                    '- The ', assessment , ' data will be computed for the "', variable, '" variable for all the assay(s).'),
                color = 'blue',
                verbose = verbose
            )
            compute.metric <- TRUE
            assay.names <- assay.names
        } else if (length(selected.assays) < length(assay.names)){
            printColoredMessage(
                message = paste0(
                    '- The ', assessment , ' data have been already computed for the "', variable, '" variable for all the ',
                    paste0(assay.names[!assay.names %in% selected.assays], collapse = ', '), ' assays.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0('- Then, the ', assessment , ' data will be computed for only the ',
                                 paste0(selected.assays, collapse = ', '),
                                 ' assay(s).'),
                color = 'blue',
                verbose = verbose
            )
            compute.metric <- TRUE
            assay.names <- droplevels(assay.names[assay.names %in% selected.assays])
        }
    }
    if(isFALSE(compute.metric)){
        return(compute.metric)
    } else
        return(list(
            compute.metric = compute.metric,
            selected.assays = selected.assays)
        )
}

