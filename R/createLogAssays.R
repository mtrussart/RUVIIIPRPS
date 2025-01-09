#' Apply log-transformation to the assay(s) of a SummarizedExperiment object.

#' @author Marie Trussart

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A single symbol or a vector of symbols specifying the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is "all," indicating that all assays in the SummarizedExperiment object will
#' be selected.
#' @param pseudo.count Numeric. A value to be added to all measurements in the assay(s) before log-transformation to prevent
#' -Inf values for measurements equal to 0. The default is 1.
#' @param replace.assays Logical. If 'TRUE', the original assay(s) will be replaced with the log-transformed values.
#' @param new.name Character. A symbol to append to the name(s) of the assay(s) for the log-transformed data. The default is 'log.'.
#' @param apply.round Logical. If 'TRUE', the measurements in the individual assays will be rounded to two decimal places.
#' The default is 'TRUE'.
#' @param verbose Logical. If 'TRUE', the function will display messages about the different steps being executed.

#' @importFrom SummarizedExperiment assays assay colData


createLogAssays <- function(
        se.obj,
        assay.names = 'all',
        pseudo.count = 1,
        replace.assays = FALSE,
        new.name = 'Log.',
        apply.round = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The createLogAssays function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if(!is.vector(assay.names))
        stop('The "assay.names" must be a vector of the assay name(s) in the "se.obj" object.')

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names))
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')

    # apply log transformation and add to the SummarizedExperiment object ####
    for (x in  levels(assay.names)) {
        printColoredMessage(message = '-- Data transformation:',
            color = 'magenta',
            verbose = verbose)
        ## log transformation ####
        if (!is.null(pseudo.count)){
            printColoredMessage(
                message = paste0('Applying log2 + ', pseudo.count,' (pseudo.count) on the ', x,' data.'),
                color = 'blue',
                verbose = verbose
                )
            temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
        } else {
            printColoredMessage(
                message = paste0('Applying log2 on the', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- log2(assay(x = se.obj, i = x))
        }
        ## round the data ####
        if(isTRUE(apply.round))
            temp.data <- round(x = temp.data, digits = 2)
        ## save the data ####
        new.assay.name <- paste0(new.name, x)
        if(isTRUE(replace.assays)){
            se.obj@assays@data[x] <- NULL
            se.obj@assays@data[[new.assay.name]] <- temp.data
        } else se.obj@assays@data[[new.assay.name]] <- temp.data
    }
    printColoredMessage(message = '------------The createLogAssays function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}


