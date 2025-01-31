#' Computes the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function calculates the the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment object and a categorical variable (i.e. batch). Then, the functions generates a line-dot plot
#' between the first cumulative PCs and the correlation coefficient to see the relationship between different PCs with the
#' variable. An ideal normalization should results a low correlation with unwanted variation variables and high correlation
#' with known biology.

#' @details
#' We use the Rozeboom squared vector correlation to quantify the strength of (linear) relationships between two sets
#' of variables, such as the first k PCs (that is 1 ≤ k ≤ 10) and dummy variables representing time, batches, plates and
#' biological variables. Not only does this quantity summarize the full set of canonical correlations, but it also reduces
#' to the familiar R2 from multiple regression when one of the variable sets contains just one element.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or a vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object for calculating the vector correlation. The default is set to "all", which indicates
#' that all assays of the SummarizedExperiment object will be selected.
#' @param variable Character. A character string indicating the name of the column in the sample annotation of the SummarizedExperiment
#' object. The variable must be categorical variable such as batches or sample subtypes.
#' @param fast.pca Logical. Indicates whether to use the fast PCA or PCA results computed by the 'computePCA' function.
#' The default is set to 'TRUE'. If 'FALSE', the 'fast.pca' argument in the 'computePCA' function should be set to 'FALSE'.
#' @param nb.pcs Numeric. A numeric value indicating the number of the first principal components (PCs) to use for
#' calculating the vector correlation. The default is set to 10. This number cannot exceed the number of PCs calculated
#' by the 'computePCA' function.
#' @param save.se.obj Logical. Indicates whether to save the results, correlation coefficients, in the metadata of the
#' SummarizedExperiment object or to output the result as a list. By default, it is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', displays the messages for the different steps of the function.

#' @return A SummarizedExperiment object or a list containing the vector correlation plots of individual assay(s) for
#' a categorical variable.


#' @importFrom SummarizedExperiment assays assay
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cancor
#' @import ggplot2
#' @export

computePCVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The computePCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # Checking the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    }
    if (is.null(variable) | is.logical(variable)) {
        stop('The "variable" cannot be empty or logical.')
    }
    if (length(variable) > 1){
        stop('The "variable" must be the name of a column in the sample annotation of the SummarizedExperiment object.')
    }
    if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    }
    if (class(se.obj@colData[[variable]]) %in% c('numeric', 'integer')) {
        stop(paste0('The "', variable, '" must be a categorical varible.'))
    }
    if (length(unique(se.obj@colData[[variable]])) < 2) {
        stop('The "variable" must have at least two levels.')
    }
    if (!is.logical(fast.pca)){
        stop('The "fast.pca" must be logical(TRUE or FALSE)')
    }
    if (length(fast.pca) > 1){
        stop('The "nb.pcs" must be a numerical value.')
    }
    if (!is.numeric(nb.pcs)){
        stop('The "nb.pcs" must be a numerical value.')
    }
    if (is.null(nb.pcs)) {
        stop('The number of PCs (nb.pcs) must be specified.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA. Re-run the computePCA with "remove.na = both"'))
    }
    if(!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical (TRUE or FALSE)')
    }
    if(!is.logical(verbose)){
        stop('The "verbose" must be logical (TRUE or FALSE)')
    }
    # Checking the assays names ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Creating dummy variables ####
    printColoredMessage(
        message =  paste0(
            '-- Creating dummy variables for the "',
            variable,
            '" variable:'),
        color = 'magenta',
        verbose = verbose
        )
    catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[variable]])
    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
    printColoredMessage(
        message =  paste0(
            '- A data frame with ',
            ncol(catvar.dummies),
            ' binary columns is created.'),
        color = 'blue',
        verbose = verbose
        )

    # Computing vector correlation ####
    printColoredMessage(
        message =  paste0(
            '-- Computing vector correlation between the first ',
            nb.pcs,
            ' PCs (cumulatively) and the "',
            variable,
            '" variable:' ),
        color = 'magenta',
        verbose = verbose
        )

    ## obtaining computed PCs from the SummarizedExperiment object ####
    printColoredMessage(
        message = paste0(
            '- Obtaining the first ',
            nb.pcs,
            ' computed PCs from the SummarizedExperiment object.'),
        color = 'blue',
        verbose = verbose
    )
    if (isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    all.pca.data <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = levels(assay.names),
        slot = 'Metrics',
        assessment = 'PCA',
        assessment.type = 'global.level',
        method = method,
        variables = 'general',
        file.name = 'data',
        sub.file.name = 'svd',
        required.function = 'computePCA',
        message.to.print = 'PCs'
    )
    ## calculate the vector correlation ####
    printColoredMessage(
        message = '- calculating the vector correlation.',
        color = 'blue',
        verbose = verbose
        )
    all.vec.corr <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Computing the vector correlation for the "',
                    x,
                    '" data:'),
                color = 'orange',
                verbose = verbose)
            vector.corr <- sapply(
                1:nb.pcs,
                function(y) {
                    cca <- cancor(x = all.pca.data[[x]]$u[, 1:y, drop = FALSE], y = catvar.dummies)
                    1 - prod(1 - cca$cor ^ 2)
                })
            return(vector.corr)
        })
    names(all.vec.corr) <- levels(assay.names)

    # Save the results ####
    printColoredMessage(
        message = '-- Saving all the vector correlation results:',
        color = 'magenta',
        verbose = verbose
        )
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the vector correlations to the "metadata" in the SummarizedExperiment object.',
            color = 'orange',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'VCA',
            method = method,
            variables = variable,
            file.name = 'vector.correlations',
            results.data = all.vec.corr
            )
        printColoredMessage(
            message = paste0(
                '* The vector correlation for the individal assay(s) are saved to the',
                ' "se.obj@metadata$metric$AssayName$VCA" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computePCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## Return only the correlation result
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- The vector correlation for the individal assay(s) are outputed as list.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computePCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.vec.corr = all.vec.corr)
    }
}
