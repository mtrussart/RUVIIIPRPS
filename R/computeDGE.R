#' Perform differential gene expression analysis.

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
#' categorical variable such as batches. If the variable has more than two levels, the function perform DGE between all
#' possible pairwise groups.
#' @param method Symbol. A symbol that indicates which DE method should be used. Options are 'limma' and 'Wilcoxon'. The
#' defalut is set to 'limma'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before compuritn DGE. The
#' default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na Symbol. Indicates whether to remove missing/NA values from either the 'assays', 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains missing/NA values will be excluded. If 'sample.annotation'
#' is selected, the samples that contains NA or missing values for each 'variables' will be excluded. By default, it is
#' set to 'both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param override.check Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
#' determine if the DGE has already been computed for the current parameters. If it has, the metric will not be recalculated.
#' The default is set to FALSE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list containing all the stats of the Wilcoxon test and if requested the
#' associated p-value histograms.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_wilcoxon_twosample
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes
#' @importFrom stats model.matrix
#' @export

computeDGE <- function(
        se.obj,
        assay.names = 'all',
        method = 'limma',
        variable,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeDGE function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if(isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'DGE',
            method = 'Wilcoxon',
            variable = variable,
            file.name = 'p.values'
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
        } else if (is.null(variable)) {
            stop('The "variable" cannot be empty.')
        } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
            stop('The "variable" must be a categorical variable.')
        } else if (length(unique(se.obj[[variable]])) < 2) {
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
                variables = variable,
                remove.na = remove.na,
                verbose = verbose
                )
        }

        # Data log transformation ####
        printColoredMessage(
            message = '-- Data log transformation:',
            color = 'magenta',
            verbose = verbose
            )
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assessment = 'DGE',
            verbose = verbose
            )
        # Apply Wilcoxon test ####
        if(method == 'wilcoxon'){
            printColoredMessage(
                message = paste0(
                    '-- Perform Wilcoxon test between all possible contrasts of the ',
                    variable,
                    ' variable.') ,
                color = 'magenta',
                verbose = verbose
            )
            all.contrasts <- combn(
                x = unique(colData(se.obj)[[variable]]),
                m = 2
            )
            all.wilcoxon.tests <- lapply(
                levels(assay.names),
                function(x){
                    printColoredMessage(
                        message = paste0('- Apply the Wilcoxon test on the "', x, '" data:'),
                        color = 'orange',
                        verbose = verbose
                    )
                    de.results <- lapply(
                        1:ncol(all.contrasts),
                        function(i){
                            printColoredMessage(
                                message = paste0(
                                    '* Wilcoxon test between the ', all.contrasts[1 , i], ' and ', all.contrasts[2 , i],'.'),
                                color = 'blue',
                                verbose = verbose
                            )
                            data.x <- all.assays[[x]][ , colData(se.obj)[[variable]] == all.contrasts[1 , i] ]
                            data.y <- all.assays[[x]][ , colData(se.obj)[[variable]] == all.contrasts[2 , i] ]
                            de.table <- matrixTests::row_wilcoxon_twosample(data.x, data.y)[ , c('obs.x', 'obs.y', 'pvalue')]
                        })
                    names(de.results) <- sapply(
                        1:ncol(all.contrasts),
                        function (x)
                            paste(all.contrasts[1 , x], all.contrasts[2 , x], sep = '&')
                    )
                    return(de.results)
                })
            names(all.wilcoxon.tests) <- levels(assay.names)
        }
        if(method == 'limma'){
            printColoredMessage(
                message = paste0(
                    '-- Perform limma between all possible contrasts of the ',
                    variable,
                    ' variable.') ,
                color = 'magenta',
                verbose = verbose
            )
            de.design <- model.matrix(~ 0 + se.obj[[variable]])
            colnames(de.design) <- levels(factor(se.obj[[variable]]))
            all.possible.contrasts <- combn(x = colnames(de.design), m = 2)
            all.possible.contrasts <- sapply(
                1:ncol(all.possible.contrasts),
                function(x) paste0(all.possible.contrasts[,x], collapse = '-')
                )
            all.possible.contrasts <- makeContrasts(
                contrasts = all.possible.contrasts,
                levels = colnames(de.design)
            )
            all.limma.tests <- lapply(
                levels(assay.names),
                function(x){
                    printColoredMessage(
                        message = paste0('- Apply limma test on the "', x, '" data:'),
                        color = 'orange',
                        verbose = verbose
                    )
                    vfit <- lmFit(object = all.assays[[x]] , design = de.design)
                    vfit  <- contrasts.fit(fit = vfit , contrasts = all.possible.contrasts)
                    efit <- eBayes(fit = vfit)
                    p.values <- as.data.frame(efit$p.value)
                    p.values.list <- lapply(
                        colnames(p.values),
                        function(x){
                           d <- p.values[ , x, drop = FALSE]
                           colnames(d) <- 'pvalue'
                           d
                           })
                    names(p.values.list) <- gsub('-', '&', colnames(p.values))
                    p.values.list
                    })
            names(all.limma.tests) <- levels(assay.names)
        }


        # Save the results ####
        printColoredMessage(
            message = '-- Save the Wilcoxon test results:',
            color = 'magenta',
            verbose = verbose
            )
        ## add results to the SummarizedExperiment object ####
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- Save all the Wilcoxon test results in the "metadata" in the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if(method == 'wilcoxon'){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = levels(assay.names),
                    assessment.type = 'gene.level',
                    assessment = 'DGE',
                    method = 'Wilcoxon',
                    variables = variable,
                    file.name = 'p.values',
                    results.data = all.wilcoxon.tests
                )
            }
            if(method == 'limma'){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = levels(assay.names),
                    assessment.type = 'gene.level',
                    assessment = 'DGE',
                    method = 'limma',
                    variables = variable,
                    file.name = 'p.values',
                    results.data = all.limma.tests
                )
            }
            printColoredMessage(
                message = paste0('- All the Wilcoxon results for indiviaul assay(s) are saved to the .',
                                 ' "se.obj@metadata$metric$AssayName$DGE" in the SummarizedExperiment object.'),
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
                verbose = verbose
                )
            printColoredMessage(message = '------------The computeDGE function finished.',
                                color = 'white',
                                verbose = verbose)
            if(method == 'wilcoxon' ){
                return(all.wilcoxon.tests = all.wilcoxon.tests)
            } else {
                return(all.limma.tests = all.limma.tests)
            }
        }
    }
    return(se.obj)
    printColoredMessage(message = '------------The computeDGE function finished.',
                        color = 'white',
                        verbose = verbose)
}







