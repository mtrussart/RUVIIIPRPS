#' Calculates relative log expression (RLE) of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function calculates relative log expression (RLE) of the assay(s) in a SummarizedExperiment object. In addition,
#' the function returns the RLE medians and interquartile ranges (IQRs) of each sample for individual assay(s).

#' @details
#' RLE plots are used to reveal trends, temporal clustering and other non-random patterns resulting from unwanted variation
#' in gene expression data. To generate RLE plots, we first form the log ratio log(yig/yg) of the raw count yig for
#' gene g in the sample labeled i relative to the median value yg of the counts for gene g taken across all samples. We
#' then generate a box plot from all the log ratios for sample i and plotted all such box plots along a line, where i
#' varies in a meaningful order, usually sample processing date. An ideal RLE plot should have its medians centered around
#' zero, and its box widths and their interquartile ranges (IQRs) should be similar in magnitude. Because of their
#' sensitivity to unwanted variation, we also examine the relationships between RLE medians and interquartile ranges with
#' potential sources of unwanted variation and individual gene expression levels in the datasets. In the absence of any
#' influence of unwanted variation in the data, we should see no such associations.

#' @references
#' Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#'
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#'
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or a vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object to calculate RLE data, medians, and interquartile ranges. The default is set to "all",
#' which indicates that all assays in the SummarizedExperiment object will be selected.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. The default is set to 'TRUE'.
#' Please note that any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value to be added as a pseudo count to all measurements of the assay(s) before applying
#' log transformation. This avoids -Inf for measurements equal to 0. The default is set to 1.
#' @param outputs.to.return Character. Specifies the type of RLE computations to be performed and the data to be returned.
#' Options include "all", "rle", "rle.med", "rle.iqr", and "rle.med.iqr". Selecting "all" returns RLE data along with
#' medians and interquartile ranges. Choosing "rle" returns only the RLE data for each assay. "rle.med" returns only the
#' RLE  medians. "rle.iqr" returns only the interquartile ranges of the RLE data. "rle.med.iqr" returns both the
#' RLE medians and interquartile ranges. The default is set to  'all'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. If 'TRUE', the 'checkSeObj'
#' function will be applied. The default is set to 'TRUE'
#' @param remove.na Character. Indicates whether to remove NA or missing values from the assays. Options are 'assays' or
#' 'none'. The default is "assays", meaning all NA or missing values in the assays will be removed before computing RLE.
#' Refer to the 'checkSeObj' function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object,
#' or output the result as a list. The default is set to 'TRUE'.
#' @param override.check Logical. When set to 'TRUE', the function verifies if the RLE has already been computed for the
#' current parameters in the SummarizedExperiment object. If already computed, the metric will not be recalculated. The
#' default is set to 'FALSE'.
#' @param verbose Logical. If 'TRUE', shows messages for each step of the function.
#'
#' @return Either a SummarizedExperiment object that contain the RLE results or a list containing the RLE data for each
#' individual assay in the SummarizedExperiment object.

#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats median
#' @import ggplot2
#' @export

computeRLE <- function(
        se.obj,
        assay.names = "all",
        apply.log = TRUE,
        pseudo.count = 1,
        outputs.to.return = 'all',
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeRLE function starts:',
                        color = 'white',
                        verbose = verbose)

    # Checking to override or not ####
    if (!is.logical(override.check)){
        stop('The "override.check" must be logical (TRUE or FALSE)')
    }
    if(isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'global.level',
            assessment = 'RLE',
            method = 'gene.median.center',
            variable = 'general',
            file.name = 'data'
        )
        if(is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if(isTRUE(compute.metric)){
        # Checking inputs ####
        if(is.list(assay.names) | is.logical(assay.names)){
            stop('The "assay.names" must be a vector of assay names(s) or set to "all".')
        }
        if (isFALSE(is.logical(apply.log))) {
            stop('The "apply.log" must be "TRUE" or "FALSE".')
        }
        if (isTRUE(apply.log)){
            if (pseudo.count < 0)
                stop('The value of "pseudo.count" cannot be negative.')
        }
        if (is.logical(outputs.to.return)) {
            stop('The "outputs.to.return" must be one of the "all", "rle", "rle.med", "rle.iqr" or "rle.med.iqr".')
        }
        if (!outputs.to.return %in% c('all', 'rle', 'rle.med', 'rle.iqr', 'rle.med.iqr')) {
            stop('The "outputs.to.return" must be one of the "all", "rle", "rle.med", "rle.iqr" or "rle.med.iqr".')
        }
        if (!remove.na %in% c('assays', 'none')) {
            stop('The "remove.na" must be one of the "assays" or "none".')
        }
        if (!is.logical(assess.se.obj)){
            stop('The "save.se.obj" must be logical (TRUE or FALSE)')
        }
        if (!is.logical(save.se.obj)){
            stop('The "save.se.obj" must be logical (TRUE or FALSE)')
        }
        if (!is.logical(verbose)){
            stop('The "verbose" must be logical (TRUE or FALSE)')
        }

        # Assess the SummarizedExperiment object ####
        if (isTRUE(assess.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = assay.names,
                variables = NULL,
                remove.na = remove.na,
                verbose = verbose
                )
        }

        # Checking the assays ####
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
        if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }

        # Data log transformation ####
        printColoredMessage( message ='-- Data log transformation:',
                             color = 'magenta',
                             verbose = verbose)
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            pseudo.count = pseudo.count,
            assessment = 'RLE',
            verbose = verbose
            )
        names(all.assays) <- levels(assay.names)

        # Compute RLE for each assay ####
        printColoredMessage(
            message = '-- Computing the RLE data for individual assay(s):',
            color = 'magenta',
            verbose = verbose
            )
        all.assays <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('- Computing the RLE for the "', x, '" data:'),
                    color = 'orange',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - matrixStats::rowMedians(all.assays[[x]])
                if (outputs.to.return == 'all'){
                    printColoredMessage(
                        message = '* obtaining the RLE data.',
                        color = 'blue',
                        verbose = verbose
                        )
                    printColoredMessage(
                        message = '* obtaining the RLE medians and interquartile ranges.',
                        color = 'grey',
                        verbose = verbose
                        )
                    rle.med <- matrixStats::colMedians(rle.data)
                    rle.iqr <- matrixStats::colIQRs(rle.data)
                    rle <- list(
                        rle.data = rle.data,
                        rle.med = rle.med,
                        rle.iqr = rle.iqr
                        )
                } else if (outputs.to.return == 'rle.data'){
                    printColoredMessage(
                        message = '* obtaining the RLE data.',
                        color = 'blue',
                        verbose = verbose)
                    rle <- list(rle.data = rle.data)
                } else if (outputs.to.return == 'rle.med'){
                    printColoredMessage(
                        message = '* obtaining the RLE data.',
                        color = 'blue',
                        verbose = verbose
                        )
                    printColoredMessage(
                        message = '* obtaining the RLE medians.',
                        color = 'blue',
                        verbose = verbose
                        )
                    rle.med <- matrixStats::colMedians(rle.data)
                    rle <- list(rle.data = rle.data, rle.med = rle.med)
                } else if (outputs.to.return == 'rle.iqr'){
                    printColoredMessage(
                        message = '* obtaining the RLE data.',
                        color = 'blue',
                        verbose = verbose)
                    printColoredMessage(
                        message = '* obtaining the interquartile ranges.',
                        color = 'blue',
                        verbose = verbose)
                    rle.iqr <- matrixStats::colIQRs(rle.data)
                    rle <- list(rle.data = rle.data, rle.iqr = rle.iqr)
                } else if (outputs.to.return == 'rle.med.iqr'){
                    printColoredMessage(
                        message = '* obtaining the RLE data.',
                        color = 'blue',
                        verbose = verbose)
                    printColoredMessage(
                        message = '* obtaining the RLE medians and interquartile ranges.',
                        color = 'blue',
                        verbose = verbose)
                    rle.med <- matrixStats::colMedians(rle.data)
                    rle.iqr <- matrixStats::colIQRs(rle.data)
                    rle <- list(rle.med = rle.med, rle.iqr = rle.iqr)
                }
                return(rle)
            })
        names(all.assays) <- levels(assay.names)

        # Save the results ####
        ## add results to the SummarizedExperiment object ####
        printColoredMessage(
            message = '-- Saving the results of the RLE computation:',
            color = 'magenta',
            verbose = verbose
            )
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- Saving all the RLE data to the "metadata" of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
                )
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                slot = 'Metrics',
                assessment = 'RLE',
                assessment.type = 'global.level',
                method = 'gene.median.center',
                file.name = 'data',
                variables = 'general',
                results.data = all.assays
                )
            printColoredMessage(
                message = paste0(
                    '* The RLE data, RLE medians, and RLE interquartile of individual assay(s) are saved to ',
                    'the "se.obj@metadata$metric$AssayName$RLE$rle.data" in SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(
                message = '------------The computeRLE function finished.',
                color = 'white',
                verbose = verbose)
            return(se.obj = se.obj)
        }

        ## return a list ####
        if (isFALSE(save.se.obj)) {
            printColoredMessage(
                message = '- The RLE data, RLE medians and RLE interquartile of the individual assay is outputed as a list.',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(
                message = '------------The computeRLE function finished.',
                color = 'white',
                verbose = verbose)
            return(rle.data = all.assays)
        }
    }
    if(isFALSE(compute.metric)){
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
}
