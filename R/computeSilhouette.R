#' Compute the average silhouette width.

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
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate Silhouette coefficient. The default is set to 'all', which indicates all the
#' assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol that indicates the column name that contain categorical variable in the
#' SummarizedExperiment object. The variable can be biological or unwanted variable.
#' @param dist.measure Symbol. A symbol that indicates which ditsance measure to be used to calculate silhouette coefficient.
#' The options are 'euclidean', 'maximum', manhattan', 'canberra', 'binary' or 'minkowski'. The default is set to 'euclidean'.
#' Refer to the function 'dist' from the 'stats' R package for more details.
#' @param fast.pca Logical. Indicates whether to use the PC calculated using fast PCA or not. The default is set to 'TRUE'.
#' The fast PCA and ordinary PCA do not affect the silhouette coefficient calculation.
#' @param nb.pcs Numeric. A numeric value indicates the number of first PCs to be used to calculated the distances between
#' samples. The default is set to 3.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object
#' or to output the result. By default it is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return Either the SummarizedExperiment object containing the computed silhouette coefficient in the metadata or a list
#' of silhouette coefficient for all specified assay(s)

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
) {
    printColoredMessage(message = '------------The computeSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)){
        stop('The "variable" cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be categorical variable.')
    }
    if(fast.pca){
        if (is.null(nb.pcs)) {
            stop('The number of PCs (nb.pcs) must be specified.')
        }
    }
    if (!dist.measure %in% c('euclidian',
                              'maximum',
                              'manhattan',
                              'canberra',
                              'binary',
                              'minkowski')) {
        stop("The dist.measure should be one of the: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.")
    }

    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Silhouette coefficients on all assays ####
    printColoredMessage(
        message = '-- Compute silhouette coefficient:',
        color = 'magenta',
        verbose = verbose)
    sil.coef <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('- Compute silhouette coefficient for the ', x, ' data.'),
                color = 'blue',
                verbose = verbose)

            printColoredMessage(
                message = paste0('- Obtain the first ', nb.pcs,' PCs.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fast.pca' %in% names(se.obj@metadata[['metric']][[x]][['PCA']]))
                    stop('To compute the Silhouette coefficient, the fast PCA must be computed first on the assay ', x, '.' )
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']][['fast.pca']][['pca.data']]$svd$u
            } else {
                if (!'pca' %in% names(se.obj@metadata[['metric']][[x]][['pca']]))
                    stop('To compute the Silhouette coefficient, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']][['pca']][['pca.data']]$svd$u
            }
            if(ncol(pca.data) < nb.pcs){
                printColoredMessage(
                    message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                    color = 'blue',
                    verbose = verbose)
                stop(paste0('The number of PCs of the assay ', x, ' are less than', nb.pcs, '.',
                            'Re-run the computePCA function with nb.pcs = ', nb.pcs, '.'))
            }
            pca.data <- pca.data[ , seq_len(nb.pcs)]
            if(!all.equal(row.names(pca.data), colnames(se.obj))){
                stop('The column names of the SummarizedExperiment object is not the same as row names of the PCA data.')
            }
            printColoredMessage(
                message = '- Calculate the distance matrix on the PCs.',
                color = 'blue',
                verbose = verbose)
            d.matrix <- as.matrix(dist(pca.data[, seq_len(nb.pcs)], method = dist.measure))
            printColoredMessage(
                message = '- Calculate the average Silhouette coefficient.',
                color = 'blue',
                verbose = verbose)
            avg.width <- summary(silhouette(as.numeric(as.factor(se.obj@colData[, variable])), d.matrix))$avg.width
            return(avg.width)
        })
    names(sil.coef) <- levels(assay.names)

    # Save the results ####
    printColoredMessage(
        message = '-- Save all the Silhouette:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save the silhouette coefficients to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)

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
            if (!'Silhouette' %in% names(se.obj@metadata[['metric']][[x]] )) {
                se.obj@metadata[['metric']][[x]][['Silhouette']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!paste0('sil.', dist.measure) %in% names(se.obj@metadata[['metric']][[x]][['Silhouette']])) {
                se.obj@metadata[['metric']][[x]][['Silhouette']][[paste0('sil.', dist.measure)]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['Silhouette']][[paste0('sil.', dist.measure)]][[variable]]$sil.coef <- sil.coef[[x]]
        }
        printColoredMessage(
            message = 'The silhouette coefficients for individual assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

    }
    ## return only the correlation result ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(sil.coef = sil.coef)
    }
}
