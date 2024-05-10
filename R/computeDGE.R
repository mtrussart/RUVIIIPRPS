#' Perform differential gene expression analysis using the Wilcoxon test.

#' @author Ramyar Molania

#' @description
#' This function performs differential gene expression analysis using Wilcoxon test between all possible pairs of a
#' categorical variable in a SummarizedExperiment object.

#' @details
#' DE analyses is performed using the Wilcoxon signed-rank test with log-transformed data e.g. raw counts, normalized data, ....
#' To evaluate the effects of the different sources of unwanted variation on the data, DE analyses is performed across
#' batches. In the absence of any batch effects, the histogram of the resulting unadjusted P values should be uniformly
#' distributed

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment object to compute the differential gene expression analysis. By default all the assays of the
#' SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol that indicates a column name of the SummarizedExperiment object that contains a
#' categorical variable such as batches. If the variable has more than two levels, the function perform DEG between all
#' possible pariwise groups.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before compuritn DGE. The
#' default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na Symbol. Indicates whether to remove missing/NA values from either the 'assays', 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains missing/NA values will be excluded. If 'sample.annotation'
#' is selected, the samples that contains NA or missing values for each 'variables' will be excluded. By default, it is
#' set to 'both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list containing all the stats of the Wilcoxon test and if requested the
#' associated p-value histograms.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export

computeDGE <- function(
        se.obj,
        assay.names = 'all',
        variable,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The computeDGE function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical variable.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop('The "variable" must have at least two levels (factors).')
    }
    if(!is.logical(apply.log)){
        stop('The "apply.log" must be "FALSE" or "TRUE".')
    }
    if (isTRUE(apply.log)){
        if (pseudo.count < 0)
            stop('The value of "pseudo.count" cannot be negative.')
    }
    if(!is.logical(assess.se.obj)){
        stop('The "assess.se.obj" must be "FALSE" or "TRUE".')
    } else if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be "FALSE" or "TRUE".')
    } else if (!is.logical(verbose)){
        stop('The "verbose" must be "FALSE" or "TRUE".')
    }

    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # Data log transformation ####
    printColoredMessage( message ='-- Data log transformation:',
                         color = 'magenta',
                         verbose = verbose)
    all.log.data <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (isTRUE(apply.log) & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('- Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- Apply log2 on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = x))
            } else if (isFALSE(apply.log)){
                printColoredMessage(
                    message = paste0('- The ', x, ' data will be used without log transformation.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = '- Please note, the assay should be in log scale before computing RLE.',
                    color = 'red',
                    verbose = verbose)
                expr.data <- assay(x = se.obj, i = x)
            }
        })
    names(all.log.data) <- levels(assay.names)

    # Apply Wilcoxon test ####
    printColoredMessage(
        message = '-- Perform Wilcoxon test:',
        color = 'magenta',
        verbose = verbose)
    all.contrasts <- combn(
        x = unique(colData(se.obj)[[variable]]),
        m = 2)
    all.wilcoxon.tests <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Apply the Wilcoxon test on the ', x, ' data:'),
                color = 'blue',
                verbose = verbose)
            de.results <- lapply(
                1:ncol(all.contrasts),
                function(i){
                    printColoredMessage(
                        message = paste0('- Wilcoxon test between the ', all.contrasts[1 , i], ' and ', all.contrasts[2 , i],'.'),
                        color = 'blue',
                        verbose = verbose)
                    data1 <- all.log.data[[x]][ , colData(se.obj)[[variable]] == all.contrasts[1 , i] ]
                    data2 <- all.log.data[[x]][ , colData(se.obj)[[variable]] == all.contrasts[2 , i] ]
                    de.table <- matrixTests::row_wilcoxon_twosample(data1, data2)[ , c('obs.x', 'obs.y', 'pvalue')]
                })
            names(de.results) <- sapply(
                1:ncol(all.contrasts),
                function (x)
                    paste(all.contrasts[1 , x], all.contrasts[2 , x], sep = '&'))
            de.results
        })
    names(all.wilcoxon.tests) <- levels(assay.names)

    # Save the results ####
    printColoredMessage(
        message = '-- Save the Wilcoxon test results:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'DGE' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['DGE']] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['DGE']][[variable]][['p.values']] <- all.wilcoxon.tests[[x]]
        }
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved to metadata@metric.',
            color = 'blue',
            verbose = verbose)

        printColoredMessage(message = '------------The computeDGE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## return the results as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeDGE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.wilcoxon.tests = all.wilcoxon.tests)
    }

}
