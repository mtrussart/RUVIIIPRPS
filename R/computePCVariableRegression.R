#' compute the linear regression.

#' @author Ramyar Molania

#' @description
#' This functions calculates the linear regression between the the first cumulative PCs of the gene expression (assays)
#' of a SummarizedExperiment object and a continuous variable (i.e. library size).

#' @details
#' R2 values of fitted linear models are used to quantity the strength of the (linear) relationships between a single
#' quantitative source of unwanted variation, such as sample (log) library size or tumor purity, and global sample
#' summary statistics, such as the first k PCs (1 ≤ k ≤ 10).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate regression analysis. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol indicating the name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable must be a continuous variable such as library size, tumor purity, ... .
#' @param fast.pca Logical. Indicates whether to use the computed fast PCA or PCA results computed by the computePCA function. The
#' default is 'TRUE'.
#' @param nb.pcs Numeric. A numeric value specifying the number of first PCs to be used to compute the vector correlation.
#' The default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the vector correlation plots to the 'metadata' of the
#' SummarizedExperiment object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object containing the computed regression for
#' the continuous variable and if requested the associated plot.

#' @importFrom stats lm var
#' @import ggplot2
#' @export

computePCVariableRegression <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computePCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (var(se.obj[[variable]]) == 0) {
        stop('The "variable" must have some variation.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a continuous varible.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
                    ' Run the checkSeObj function with "remove.na = both"',
                    ', then "computePCA"-->"computePCVariableRegression".'))
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Compute the regression on all assays ####
    printColoredMessage(
        message =  paste0('-- Compute linear regression between the first ',
                          nb.pcs, ' PCs (cumulatively) and the "', variable, 'variable:' ),
        color = 'magenta',
        verbose = verbose)
    ## obtain computed PCs from the SummarizedExperiment object ####
    printColoredMessage(
        message = paste0('- Obtain the first ', nb.pcs, ' computed PCs from the SummarizedExperiment object.'),
        color = 'blue',
        verbose = verbose
    )
    if(isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    all.pca.data <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = assay.names,
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
    ## obtain linear regression and obtain Rseq ####
    all.r.squared <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('- Compute the linear regression for the "', x , '" data:'),
                color = 'orange',
                verbose = verbose)
            printColoredMessage(
                message = '* compute the R squared of the linear regression analysis.',
                color = 'blue',
                verbose = verbose)
            r.squared <- sapply(
                1:nb.pcs,
                function(y) lm.ls <- summary(lm(se.obj@colData[[variable]] ~ all.pca.data[[x]]$u[, 1:y]))$r.squared)
            return(r.squared)
        })
    names(all.r.squared) <- levels(assay.names)

    # Save the results ####
    printColoredMessage(
        message = '-- Save all the R squared of the linear regression analysis:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save all the R squared of the linear regression analysis for each assay(s) the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'global.level',
            assessment = 'LRA',
            method = method,
            variables = variable,
            file.name = 'r.squared',
            results.data = all.r.squared
            )
        printColoredMessage(
            message = paste0('- The the R squared for the individal assay(s) are saved to the',
                             ' "se.obj@metadata$metric$AssayName$LRA$Varible$rseq" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computePCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    ## save the  regression r squared as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- The regression r squared valuse for individual assay(s) are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computePCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.r.squared = all.r.squared)
    }
}
