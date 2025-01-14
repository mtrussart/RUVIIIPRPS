#' Computes the average silhouette width using PCs. -- Finalized

#' @author Ramyar Molania

#' @description
#' This function calculates the average silhouette width for each categorical variables such as sample biological subtypes,
#' batches, etc.. The distance matrix is based on the principal components.

#' @details
#' The Silhouette coefficient analysis is used to assess the separation of categorical samples groups. We use first few
#' principal components to calculate both the similarity between one sample and the other samples in each cluster and the
#' separation between samples in different clusters. A better normalization method will lead to higher and lower silhouette
#' coefficients for biological and batch variables, respectively.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or a vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object for which to calculate the Silhouette coefficient. The default is set to 'all', which
#' indicates that all the assays of the SummarizedExperiment object will be selected.
#' @param variable Character. A character string indicating the column name that contains the categorical variable in the
#' SummarizedExperiment object. The variable can represent a biological or unwanted variable.
#' @param dist.measure Character. A character string indicating which distance measure to use for calculating the Silhouette
#' coefficient. The options are 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', or 'minkowski'. The default is
#' set to 'euclidean'. Refer to the 'dist' function from the 'stats' R package for more details.
#' @param fast.pca Logical. Indicates whether to use the PCs calculated using fast PCA. The default is set to 'TRUE'.
#' Fast PCA and ordinary PCA do not affect the Silhouette coefficient calculation.
#' @param nb.pcs Numeric. A numeric value indicating the number of principal components (PCs) to use when calculating the
#' distances between samples. The default is set to 3.
#' @param save.se.obj Logical. Indicates whether to save the silhouette coefficient  results in the metadata of the
#' SummarizedExperiment object or to output the result. The default is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return Either the SummarizedExperiment object containing the computed Silhouette coefficient stored in the metadata
#' or a list of Silhouette coefficients for all specified assay(s).


#' @importFrom SummarizedExperiment assays assay
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @import ggplot2
#' @export

computeSilhouette <- function(
        se.obj,
        assay.names = 'all',
        variable,
        dist.measure = 'euclidian',
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking  the function inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (is.null(variable) | is.logical(variable)){
        stop('The "variable" cannot be empty or logical.')
    }
    if (length(variable) > 1){
        stop('The "variable" must be a categorical variable.')
    }
    if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    }
    if (class(se.obj@colData[[variable]]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical variable.')
    }
    if (length(unique(se.obj@colData[, variable])) < 2) {
        stop('The "variable" must have at least two levels.')
    }
    if (!is.logical(fast.pca)){
        stop('The "fast.pca" must be logical(TRUE or FALSE)')
    }
    if (length(nb.pcs) > 1){
        stop('The "nb.pcs" must be a numerical value.')
    }
    if (!is.numeric(nb.pcs)){
        stop('The "nb.pcs" must be a numerical value.')
    }
    if (is.null(nb.pcs)) {
        stop('The number of PCs (nb.pcs) must be specified.')
    }
    if (!dist.measure %in% c('euclidian',
                              'maximum',
                              'manhattan',
                              'canberra',
                              'binary',
                              'minkowski')) {
        stop("The dist.measure should be one of the: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.")
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical(TRUE or FALSE)')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical(TRUE or FALSE)')
    }

    # Checking assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Compute silhouette coefficients on all assays ####
    printColoredMessage(
        message = paste0(
            '-- Computing silhouette coefficient using the first ',
            nb.pcs,
            ' PCS and the ',
            variable,
            ' variable.') ,
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'

    ## obtaining the PCA data from SummarizedExperiment ####
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
    ## computing average silhouette coefficients for all assay(s)  ####
    all.sil.coef <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '- Computing silhouette coefficient for the "',
                    x,
                    '" data:'),
                color = 'orange',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '* obtaining the first ',
                    nb.pcs,
                    ' computed PCs.'),
                color = 'blue',
                verbose = verbose
                )
            pca.data <- all.pca.data[[x]]$u
            if (ncol(pca.data) < nb.pcs){
                printColoredMessage(
                    message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                    color = 'blue',
                    verbose = verbose)
                stop(paste0(
                    'The number of PCs of the assay ',
                    x,
                    ' are less than',
                    nb.pcs,
                    '.',
                    'Re-run the computePCA function with nb.pcs = ',
                    nb.pcs,
                    '.'))
            }
            pca.data <- pca.data[ , seq_len(nb.pcs), drop = FALSE]
            if (!all.equal(row.names(pca.data), colnames(se.obj))){
                stop('The column names of the SummarizedExperiment object is not the same as row names of the PCA data.')
            }
            printColoredMessage(
                message = '* calculating the distance matrix on the PCs.',
                color = 'blue',
                verbose = verbose
                )
            d.matrix <- as.matrix(dist(pca.data[, seq_len(nb.pcs)], method = dist.measure))
            printColoredMessage(
                message = '* calculating the average Silhouette coefficient.',
                color = 'blue',
                verbose = verbose
                )
            avg.width <- summary(silhouette(as.numeric(as.factor(se.obj@colData[, variable])), d.matrix))$avg.width
            return(avg.width)
        })
    names(all.sil.coef) <- levels(assay.names)

    # Save the results ####
    printColoredMessage(
        message = '-- Saving all the Silhouette coefficients results :',
        color = 'magenta',
        verbose = verbose
        )
    ## adding results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving the silhouette coefficients of each assay(s) to the "metadata" in the SummarizedExperiment object:',
            color = 'orange',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'Silhouette',
            method = paste0('sil.', dist.measure),
            variables = variable,
            file.name = 'silhouette.coeff',
            results.data = all.sil.coef
        )
        printColoredMessage(
            message = paste0(
                '* The silhouette coefficients of the induvial assay(s) is saved to the .',
                ' ".se.obj@metadata$metric$RawCount$Silhouette" in the SummarizedExperiment objec.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    ## returning only the correlation result ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.sil.coef = all.sil.coef)
    }
}
