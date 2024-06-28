#' Compute adjusted rand index (ARI).

#' @author Ramyar Molania

#' @description
#' This functions computes the adjusted rand index for given a categorical variable using the first PCs of the assay(s)
#' in a SummarizedExperiment object.

#' @details
#' The ARI is the corrected-for-chance version of the Rand index. The ARI measures the percentage of matches between
#' two label lists. We use the ARI to assess the performance of normalization methods in terms of sample subtype
#' separation and batch mixing. We first calculate PCs and use the first PCs to perform ARI.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute ARI. The default is set to 'all', so all the assays of the SummarizedExperiment
#' object will be selected.
#' @param variable Symbol. Indicates the column name that contain categorical variable in the SummarizedExperiment object.
#' The variable can be biological or unwanted variable.
#' @param clustering.method Symbol. A symbol that indicates which clustering methods should be applied on the PCs to
#' calculate the ARI. The options are 'mclust' or 'hclust' methods. The default is set to 'hclust'.
#' @param hclust.method Symbol. Indicate the agglomeration method to be used when the  'clustering.method' is 'hclust'
#' method. The options are: 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA),
#' 'median' (= WPGMC) or centroid' (= UPGMC). Refer to the 'hclust' function from the stats R package for more details.
#' @param hclust.dist.measure Symbol. Indicates the distance measure to be used in the 'dist' function. The options are
#' 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'. Refer to the 'dist' function from the stats
#' R package for more details.
#' @param fast.pca Logical. Indicates whether to use the PC calculated using fast PCA or not. The default is set to 'TRUE'.
#' The fast PCA and ordinary PCA do not affect the silhouette coefficient calculation. Refer to detail for more information.
#' @param nb.pcs Numeric. The number of first PCs to be used to calculated the distances between samples. The default is
#' set to 3.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @importFrom stats cutree hclust dist
#' @import ggplot2
#' @export

computeARI <- function(
        se.obj,
        assay.names = 'all',
        variable,
        clustering.method = 'hclust',
        hclust.method = 'complete',
        hclust.dist.measure = 'euclidian',
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeARI function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s) or "assay.names = all".')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (length(unique(se.obj[[variable]])) == 1) {
        stop('The "variable" must have at least two levels.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical varible.')
    }
    if (is.null(nb.pcs)) {
        stop('The number of PCs (nb.pcs) must be specified')
    }
    if(length(clustering.method) > 1){
        stop('The "clustering.method" must be only one of the "mclust" or "hclust".')
    }
    if(is.null(clustering.method)){
        stop('The "clustering.method" cannot be empty.')
    } else if(!clustering.method %in% c('mclust', 'hclust')){
        stop('The "clustering.method" method must be one of the "mclust" or "hclust".')
    }
    if(clustering.method == 'hclust'){
        if(is.null(hclust.dist.measure)){
            stop('The "hclust.dist.measure" cannot be empty when the "clustering.method = hclust".')
        } else if (!hclust.dist.measure %in% c('euclidian',
                                               'maximum',
                                               'manhattan',
                                               'canberra',
                                               'binary',
                                               'minkowski')) {
            stop('The "hclust.dist.measure" should be one of the:"euclidean","maximum","manhattan","canberra","binary" or "minkowski".')
        }
        if(is.null(hclust.method)){
            stop('The "hclust.method" cannot be when the "clustering.method = hclust".')
        } else if (!hclust.method %in% c('complete',
                                         'ward.D',
                                         'ward.D2',
                                         'single',
                                         'average',
                                         'mcquitty',
                                         'median',
                                         'centroid')) {
            stop('The hclust.method should be one of the:"complete","ward.D","ward.D2","single","average", "mcquitty", "median" or "centroid".')
        }
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
                    ' Run the checkSeObj function with "remove.na = both"',
                    ', then "computePCA"-->"computeARI".'))
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Compute ARI on all the assays ####
    printColoredMessage(
        message = paste0('-- Compute adjusted rand index (ARI) using the first ',
                         nb.pcs, ' PCS for the ', variable, ' variable.') ,
        color = 'magenta',
        verbose = verbose)
    if(isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'

    ## obtain the PCA data from SummarizedExperiment ####
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

    ## compute ARI for all assay(s)  ####
    all.ari <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('- Compute the ARI for the "', x, '" data:'),
                color = 'orange',
                verbose = verbose)
            printColoredMessage(
                message = paste0('* obtain the first ', nb.pcs, ' computed PCs.'),
                color = 'blue',
                verbose = verbose)
            pca.data <- all.pca.data[[x]]$u
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

            if(clustering.method == 'mclust'){
                printColoredMessage(
                    message = '* cluster the PCs using the mclust function.',
                    color = 'blue',
                    verbose = verbose)
                bic <- mclustBIC(data = pca.data)
                mod <- Mclust(data = pca.data, x = bic, G = length(unique(se.obj@colData[[variable]])) )
                printColoredMessage(
                    message = '* calculate the adjusted rand index.',
                    color = 'blue',
                    verbose = verbose)
                ari <- adjustedRandIndex(mod$classification, se.obj@colData[, variable])
            } else {
                printColoredMessage(
                    message = '* cluster the PCs using the hclust function.',
                    color = 'blue',
                    verbose = verbose)
                clusters <- cutree(
                    tree = hclust(d = dist(x = pca.data, method = hclust.dist.measure), method = hclust.method),
                    k = length(unique(se.obj@colData[, variable])))
                printColoredMessage(
                    message = '* calculate the adjusted rand index.',
                    color = 'blue',
                    verbose = verbose)
                ari <- adjustedRandIndex(clusters, se.obj@colData[, variable])
            }
            return(ari)
        })
    names(all.ari) <- levels(assay.names)

    # save the results ####
    printColoredMessage(
        message = '-- Save all the ARI:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save the ARIs of each assay(s) to the "metadata" in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose)
        if(clustering.method == 'mclust'){
            method <- 'mclust'
        } else method <- paste0('hclust.', hclust.method, '.', hclust.dist.measure)
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'global.level',
            assessment = 'ARI',
            method = method,
            variables = variable,
            file.name = 'ari',
            results.data = all.ari
            )
        printColoredMessage(
            message = paste0('- The ARI of induvial assay(s) is saved to the .',
                             ' ".se.obj@metadata$metric$RawCount$ARI" in the SummarizedExperiment objec.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## return only the adjusted rand index results ####
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '-Save the ARIs as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.ari = all.ari)
    }
}
