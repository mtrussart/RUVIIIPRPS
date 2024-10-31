#' Perform principal component analysis (PCA) using singular value decomposition (SVD).

#' @author Ramyar Molania

#' @description
#' This function uses singular value decomposition to perform principal component on the assay(s) in a SummarizedExperiment
#' object. The function provides fast singular value decomposition using the BiocSingular R package.

#' @details
#' The PCs (in this context also called singular vectors) of the sample × transcript array of log counts are the linear
#' combinations of the transcript measurements having the largest, second largest, third largest, etc., variation,
#' standardized to be of unit length and orthogonal to the preceding components. Each will give a single value for
#' each sample.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute PCA. The default is "all, which indicates all the assays of the
#' SummarizedExperiment object will be selected.
#' @param fast.pca Logical. Indicates whether to calculate a specific number of left singular vectors instead of the
#' full possible vectors to speed up the process. The default is 'TRUE'.
#' @param nb.pcs Numeric. The number of first left singular vectors to be calculated for the fast PCA process. The
#' default is 10.
#' @param center Logical. Indicates whether to scale the data or not. If center is 'TRUE', then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before computing the SVD. The
#' default is 'TRUE'. The data must be in log scale before computing the SVD.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before applying log transformation.
#' The default is 1. This argument cannot be NULL or negative.
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package for more details.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object. The default is TRUE.
#' We refer to the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the SVD results in the metadata of the SummarizedExperiment object
#' or to output the results as list. By default it is set to 'TRUE'.
#' @param remove.na To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' @param override.check Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
#' determine if the PCA has already been computed for the current parameters. If it has, the metric will not be recalculated.
#' The default is set to FALSE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object or a list that containing the singular value decomposition results and the
#' percentage variation of each PCs.

#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular runSVD bsparam
#' @import ggplot2
#' @export

computePCA <- function(
        se.obj,
        assay.names = 'all',
        fast.pca = TRUE,
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        svd.bsparam = bsparam(),
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computePCA function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if(isTRUE(override.check)){
        if(isTRUE(fast.pca)){
            method <- 'fast.svd'
        } else method <- 'ordinary.svd'
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'global.level',
            assessment = 'PCA',
            method = method,
            variable = 'general',
            file.name = 'data',
            verbose = verbose
        )
        if(is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if(isTRUE(compute.metric)){
        # Check the inputs ####
        if (is.null(assay.names)) {
            stop('The "assay.names" cannot be empty.')
        }
        if (isTRUE(fast.pca) & is.null(nb.pcs)) {
            stop('To perform fast PCA, the number of PCs (left singular vectors) must be specified.')
        } else if (isTRUE(fast.pca) & nb.pcs == 0) {
            stop('To perform fast PCA, the number of PCs (left singular vectors) must be specified.')
        }
        if (isTRUE(apply.log)){
            if(pseudo.count < 0){
                stop('The value of "pseudo.count" cannot be negative.')
            }
            if (is.null(pseudo.count)){
                stop('A value for the "pseudo.count" must be specified.')
            }
        }
        if(!remove.na %in% c('assays','none')){
            stop('The "remove.na" must be on of the "assays" or "none"')
        }
        if (isTRUE(scale)) {
            printColoredMessage(
                message = 'Note: highly recommend not to scale the data before computing the PCA.',
                color = 'red',
                verbose = verbose)
        }

        # Check the assays ####
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
        if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }
        # Assess the SummarizedExperiment object ####
        if (isTRUE(assess.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = NULL,
                remove.na = remove.na,
                verbose = verbose)
        }
        # Data transformation ####
        printColoredMessage(
            message = '-- Data log transformation:',
            color = 'magenta',
            verbose = verbose)
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assessment = 'computing "PCA"',
            verbose = verbose
        )
        # Compute SVD ####
        printColoredMessage(
            message = '-- Perform singular value decomposition (SVD):',
            color = 'magenta',
            verbose = verbose
        )
        ## compute fast SVD ####
        if (isTRUE(fast.pca)) {
            printColoredMessage(
                message = paste0(
                    '- Perform "fast" singular value decomposition with scale = ',
                    scale, ' and center = ', center, '.'),
                color = 'orange',
                verbose = verbose)
            if (is.null(svd.bsparam))
                svd.bsparam <- bsparam()
            all.sv.decomposition <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0('* perform fast SVD on the "', x , '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    sv.dec <- BiocSingular::runSVD(
                        x = t(all.assays[[x]]),
                        k = nb.pcs,
                        BSPARAM = svd.bsparam,
                        center = center,
                        scale = scale)
                    rownames(sv.dec$u) <- colnames(se.obj)
                    rownames(sv.dec$v) <- row.names(se.obj)
                    percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                    percentage <- sapply(
                        seq_along(percentage),
                        function(i) round(percentage [i], 1))
                    return(list(svd = sv.dec, percentage.variation = percentage))
                })
            printColoredMessage(
                message = paste0(
                    '- Note: in the fast svd analysis, the percentage of variation of PCs will be ',
                    'computed proportional to the highest selected number of PCs (left singular vectors), not on all the PCs.'),
                color = 'red',
                verbose = verbose)
            names(all.sv.decomposition) <- levels(assay.names)
        }
        ## compute ordinary SVD ####
        if (isFALSE(fast.pca)) {
            printColoredMessage(
                message = paste0(
                    '- Perform "ordinary" singular value decomposition with scale = ',
                    scale, ' and center = ', center, '.'),
                color = 'orange',
                verbose = verbose
            )
            all.sv.decomposition <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0('* Perform singular value decomposition on the "', x , '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    sv.dec <- svd(scale(
                        x = t(all.assays[[x]]),
                        center = center,
                        scale = scale)
                    )
                    rownames(sv.dec$u) <- colnames(se.obj)
                    rownames(sv.dec$v) <- row.names(se.obj)
                    percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                    percentage <- sapply(
                        seq_along(percentage),
                        function(i) round(percentage [i], 1))
                    return(list(svd = sv.dec, percentage.variation = percentage))
                })
            names(all.sv.decomposition) <- levels(assay.names)
        }

        # Save all the results ####
        printColoredMessage(
            message = '-- Save the SVD results:',
            color = 'magenta',
            verbose = verbose)
        ## add the results to the SummarizedExperiment object ####
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- Save all the SVD results to the "metadata" of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if(isTRUE(fast.pca)){
                method <- 'fast.svd'
            } else method <- 'ordinary.svd'
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'global.level',
                assessment = 'PCA',
                method = method,
                file.name = 'data',
                variables = 'general',
                results.data = all.sv.decomposition
            )
            printColoredMessage(
                message = paste0('- The SVD results of individual assay (s) are saved to the',
                                 ' "se.obj@metadata$metric$AssayName$global.level$PCA$', method, '$data" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)

            printColoredMessage(message = '------------The computePCA function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj)
        }

        if (isFALSE(save.se.obj)) {
            ## return a list ####
            printColoredMessage(
                message = '- The SVD results of individual assays are outputed as a list.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(message = '------------The computePCA function finished.',
                                color = 'white',
                                verbose = verbose)
            return(all.sv.decompositions = all.sv.decomposition)
        }
    } else {
        printColoredMessage(message = '------------The computePCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)

    }
}
