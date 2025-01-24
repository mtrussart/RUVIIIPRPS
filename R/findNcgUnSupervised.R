#' Finds negative control genes (NCG) using an unsupervised approach.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function identifies a set of genes as negative control genes in cases where no biological variations are known.
#' The function uses gene-level correlation, ANOVA, and MAD analyses across and between sample groups to find a set
#' of genes as negative control genes for the RUV-III-PRPS normalization.
#'
#' @details
#' The function initially utilizes gene-level correlation and ANOVA to identify genes significantly influenced by continuous
#' and categorical sources of variation, respectively. Subsequently, it conducts Median Absolute Deviation (MAD) analysis
#' on each gene within homogeneous sample groups, considering unwanted variables, to pinpoint genes highly variable due
#' to biological factors. Lastly, various methods are employed to consolidate the statistical findings, ultimately determining
#' an appropriate set of genes as negative control genes for RUV-III-PRPS normalization.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. The name of the assay in the SummarizedExperiment object to be used for RUV-III-PRPS
#' normalization.
#' @param uv.variables Character. A character vector indicating the names of the columns in the SummarizedExperiment object
#' that contain unwanted variables. If all unwanted variation is unknown, the 'identifyUnknownUV' function can estimate them.
#' @param nb.ncg Numeric. Specifies the number of genes to be chosen as negative control genes (NCG) when the '
#' ncg.selection.method' is set to 'auto'. This value corresponds to a fraction of the total genes in the SummarizedExperiment
#'  object. The default is 0.1.
#' @param ncg.selection.method Character. Specifies the method used to select negative control genes (NCG). The available options
#' are 'prod', 'sum', 'average', 'non.overlap', 'quantile', and 'auto'. The default method is 'non.overlap'.
#' @param top.rank.bio.genes Numeric. Specifies the fraction of top-ranked genes that are highly affected by biological variation.
#' These genes have the highest MAD values. The default is 0.5.
#' @param top.rank.uv.genes Numeric. Specifies the fraction of top-ranked genes that are highly affected by unwanted variation
#' variables. These genes have the highest F-statistics and correlation coefficients. The default is 0.5.
#' @param uv.percentile Numeric. The percentile cut-off for selecting genes that are highly affected by unwanted variation.
#' The default is 0.8.
#' @param bio.percentile Numeric. The percentile cut-off for selecting genes that are highly affected by biological variation.
#' The default is 0.8.
#' @param grid.nb Numeric. Specifies the number of genes to consider during the grid search when the 'ncg.selection.method'
#' is set to 'auto'. The default is 20.
#' @param grid.group Character. Specifies whether the grid search should focus on biological, unwanted, or both factors. The
#' options are 'bio', 'uv', or 'both'.
#' @param grid.direction Character. Specifies whether the grid search should proceed in decreasing or increasing order.
#' The options are 'increase' or 'decrease'.
#' @param normalization Character. Specifies the normalization method used for library size adjustment before identifying
#' genes that are highly affected by biological variation. The default is 'CPM'. Refer to the 'applyOtherNormalization'
#' function for further details.
#' @param regress.out.variables Character. Specifies which variables to regress out from the data prior to analysis.
#' @param min.sample.for.mad Numeric. The minimum number of samples required to perform MAD analysis on each gene within homogeneous
#' sample groups. The default is 3.
#' @param min.sample.for.aov Numeric. The minimum number of samples required for ANOVA analysis between categorical sources
#' of variation with individual gene expression. The default is 3.
#' @param min.sample.for.correlation Numeric. The minimum number of samples required for correlation analyses between continuous
#' sources of unwanted variation and individual gene expression. The default is 10.
#' @param corr.method Character. The correlation method to use for the analysis. Options are 'pearson' or 'spearman'. The default
#' is 'spearman'.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation analysis. The default is 0.05.
#' @param rho Numeric. The hypothesized correlation value to be used in the hypothesis testing. The default is 0.
#' @param anova.method Character. The ANOVA method to use. Options are 'aov' or 'welch'. The default is 'aov'. Refer to the function
#' 'computeGenesVariableAnova' for more details.
#' @param clustering.method Character. Specifies the clustering method to use for grouping continuous sources of unwanted variation.
#' The available options are 'kmeans', 'cut', and 'quantile'. The default is 'kmeans'.
#' @param nb.clusters Numeric. The number of clusters for grouping continuous sources of unwanted variation. The default is 3.
#' @param assess.ncg Logical. Specifies whether to assess the performance of the selected NCGs. This involves principal component
#' analysis on the selected NCGs and exploration of R^2 or vector correlation between the first 'nb.pcs' principal components
#' and the biological and unwanted variables.
#' @param variables.to.assess.ncg Character. The column names of the SummarizedExperiment object that contain variables to be
#' assessed for their association with the selected NCGs. By default, all variables in 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Specifies the number of principal components from the selected NCGs to be used for performance assessment.
#' The default is 5.
#' @param center Logical. Specifies whether to center the data before applying Singular Value Decomposition (SVD). The
#' default is set to 'TRUE'.
#' @param scale Logical. Specifies whether to scale the data before applying SVD. The default is FALSE.
#' @param apply.log Logical. Specifies whether to apply a log-transformation to the data. The default is TRUE.
#' @param pseudo.count Numeric. The pseudo count to be added to all measurements before log transformation.
#' @param assess.se.obj Logical. Specifies whether to assess the SummarizedExperiment object. The default is TRUE.
#' @param remove.na Character. Specifies whether to remove missing values (NA) from the 'assays', the 'sample.annotation',
#' 'both', or 'none'. The default is set to 'both'.
#' @param save.se.obj Logical. Specifies whether to save the result in the metadata of the SummarizedExperiment object or
#' output the result. The default is set to 'TRUE'.
#' @param output.name Character. Specifies the name for the output in the metadata of the SummarizedExperiment object. The
#' default is set to NULL, meaning the function will generate a name based on the specified argument.
#' @param ncg.group Character. Specifies the name of the group of NCGs. If NULL, the function will generate a name using
#' "ncg|unsupervised".
#' @param plot.output Logical. If TRUE, a plot of the NCG assessment will be displayed while the function is running.
#' @param use.imf Logical. Specifies whether to use the intermediate file. The default is FALSE.
#' @param save.imf Logical. Specifies whether to save the intermediate file. If TRUE, the function will save the results
#' of the ANOVA. If users want to change parameters such as 'nb.ncg', 'ncg.selection.method', 'top.rank.bio.genes', or
#' 'top.rank.uv.genes', the analysis will not be recalculated.
#' @param imf.name Character. Specifies the name for saving the intermediate file. If NULL, the function will generate a name.
#' @param verbose Logical. If TRUE, process messages will be displayed.
#'
#' @return A SummarizedExperiment object containing the selected negative control genes or a logical vector indicating
#' the selected genes.

#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr progress_estimated
#' @importFrom fastDummies dummy_cols
#' @importFrom BiocSingular bsparam
#' @importFrom BiocSingular runSVD
#' @importFrom tidyr pivot_longer
#' @importFrom Rfast correls
#' @importFrom stats aov
#' @import ggplot2
#' @export

findNcgsUnSupervised <- function(
        se.obj,
        assay.name,
        uv.variables,
        nb.ncg = 0.1,
        ncg.selection.method = 'non.overlap',
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        normalization = 'CPM',
        regress.out.variables = NULL,
        min.sample.for.mad = 3,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        scale = FALSE,
        center = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        output.name = NULL,
        ncg.group = NULL,
        plot.output = TRUE,
        use.imf = FALSE,
        save.imf = FALSE,
        imf.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgsUnSupervised function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check functions inputs ####
    if (length(assay.name) > 1 | is.logical(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
        stop('Some or all the "uv.variables" cannot be found in the SummarizedExperiment object.')
    }
    if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
    }
    if (!ncg.selection.method %in% c('auto', 'non.overlap')){
        stop('The "ncg.selection.method" must be one of "auto" or "non.overlap".')
    }
    if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" msut be a positve value  0 < top.rank.bio.genes < 1.')
    }
    if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    }
    if (grid.nb < 1 | grid.nb > nrow(se.obj)){
        stop(paste0('The "grid.nb" must be a positve value  0 < grid.nb < ', nrow(se.obj), '.'))
    }
    if (!grid.group %in% c('bio', 'uv', 'both')){
        stop('The "grid.group" must be one of "bio", "uv" or "non.overlap".')
    }
    if (!grid.direction %in% c('increase', 'decrease', 'auto')){
        stop('The "grid.direction" must be one of "increase", "decrease" or "auto".')
    }
    if (is.null(min.sample.for.aov)){
        stop('The "min.sample.for.aov" cannot be empty.')
    }
    if (min.sample.for.aov <= 2){
        stop('The "min.sample.for.aov" should be at least 3.')
    }
    if (is.null(min.sample.for.correlation)){
        stop('The min.sample.for.correlation cannot be empty.')
    }
    if (min.sample.for.correlation >= ncol(se.obj) | min.sample.for.correlation < 3){
        stop('The "min.sample.for.correlation" msut be more than 2 and less than the total number of samples in the data.')
    }
    if (!anova.method %in% c('aov', 'welch')){
        stop('The anova.method must be one of the "aov" or "welch".')
    }
    if (isFALSE(is.logical(assess.ncg))){
        stop('The "assess.ncg" must be "TRUE" or "FALSE.')
    }
    if (length(nb.pcs) > 1){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE.')
    }
    if (length(pseudo.count) > 1){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (pseudo.count < 0){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (isFALSE(is.logical(assess.se.obj))) {
        stop('The "assess.se.obj" must be "TRUE" or "FALSE.')
    }
    if (is.null(assess.se.obj)) {
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        } else if (!is.null(variables.to.assess.ncg)) {
            if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg))) {
                stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.null(regress.out.variables)){
        if (isTRUE(sum(regress.out.variables %in% colnames(colData(se.obj))) != length(regress.out.variables))) {
            stop('All or some of "regress.out.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (isTRUE(ncg.selection.method == 'quantile')){
        if(is.null(bio.percentile) | is.null(uv.percentile))
            stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
        if(bio.percentile > 1 | bio.percentile < 0)
            stop('The "bio.percentile" must be a postive value between 0 and 1.')
        if(uv.percentile > 1 | uv.percentile < 0)
            stop('The "uv.percentile" must be a postive value between 0 and 1.')
    }

    # Checking the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(uv.variables, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Data transformation and normalization ####
    printColoredMessage(
        message = '-- Applying data normalization and transformation:',
        color = 'magenta',
        verbose = verbose
        )
    ## normalization ####
    if (!is.null(normalization)){
        expr.data.nor <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
    }
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 + ',
                pseudo.count,' (pseudo.count) on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without any transformation.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- assay(x = se.obj, i = assay.name)
    }

    # Finding negative control genes unsupervised ####
    if (isFALSE(use.imf)){
        printColoredMessage(
            message = '-- Finding a subset of genes as negative control genes:',
            color = 'magenta',
            verbose = verbose
            )
        ## identify genes that are highly affected by unwanted variation ####
        printColoredMessage(
            message = '-- Finding genes that are highly affected by each specified source(s) of unwanted variation:',
            color = 'orange',
            verbose = verbose
            )

        ### find classes of variables ####
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]]))
            )
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]

        ### anova between genes and categorical sources of unwanted variation ####
        if (isTRUE(length(categorical.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Performing ANOVA between individual gene-level ',
                    'expression and each categorical source of unwanted variation: ',
                    paste0(categorical.uv, collapse = ' & '),
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            anova.genes.uv <- lapply(
                categorical.uv,
                function(x) {
                    keep.samples <- findRepeatingPatterns(
                        vec = colData(se.obj)[[x]],
                        n.repeat = min.sample.for.aov
                        )
                    if( isTRUE(length(keep.samples) == 0)){
                        stop(paste0(
                            'There are not enough samples to perform ANOVA between individual gene expression and the ',
                            x,
                            ' variable. Possible solutions is to lower min.sample.for.aov or remove the "',
                            x,
                            '" from the uv.variables and re-run the function.'))
                    } else if (isTRUE(length(keep.samples) == 1)){
                        stop(paste0(
                            'There is only a single batch from in the ',
                            x,
                            ' variable that have enough samples ',
                            min.sample.for.aov,
                            ' (min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove the "',
                            x,
                            '" from the uv.variables and re-run the function.'
                        ))
                    } else if (isTRUE(length(keep.samples) != length(unique(colData(se.obj)[[x]]))) ){
                        not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                        printColoredMessage(
                            message = paste0(
                                'Note, the ',
                                paste0(not.coverd, collapse = '&'),
                                ' batches do not have enough samples for the ANOVA analysis.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                    if (anova.method == 'aov'){
                        anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples])
                            )
                    } else if (anova.method == 'welch'){
                        anova.gene.batch <- as.data.frame(row_oneway_welch(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples])
                            )
                    }
                    set.seed(2233)
                    anova.gene.batch$ranked.genes <- rank(
                        -anova.gene.batch[ , 'statistic'],
                        ties.method = 'random'
                        )
                    anova.gene.batch
                })
            names(anova.genes.uv) <- categorical.uv
        } else anova.genes.uv <- NULL

        ### correlation between genes and continuous sources of unwanted variation ####
        if (isTRUE(length(continuous.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Performing ',
                    corr.method,
                    ' correlation between individual gene-level ',
                    'expression and each continuous source of unwanted variations: ',
                    paste0(continuous.uv, collapse = '&'),
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            if(isTRUE(ncol(se.obj) <= min.sample.for.correlation)){
                stop(paste0(
                    'There are not enough samples (min.sample.for.correlation:',
                    min.sample.for.correlation,
                    ') to perform correlation analysis.',
                    ' A possible solution is to lower the "min.sample.for.correlation" value.'))
            }
            corr.genes.uv <- lapply(
                continuous.uv,
                function(x) {
                    corr.genes.var <- as.data.frame(correls(
                        y = se.obj@colData[, x],
                        x = t(expr.data),
                        type = corr.method,
                        a = a ,
                        rho = rho)
                        )
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 3),
                        corr.genes.var[, 5, drop = FALSE]
                        )
                    set.seed(2233)
                    colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                    corr.genes.var$ranked.genes <- rank(
                        -abs(corr.genes.var[, 'statistic']),
                        ties.method = 'random'
                        )
                    row.names(corr.genes.var) <- row.names(expr.data)
                    corr.genes.var
                })
            names(corr.genes.uv) <- continuous.uv
        } else corr.genes.uv <- NULL

        ## finding genes that are highly affected by possible biological variation ####
        printColoredMessage(
            message = '-- Finding genes that are potentially highly affected by biological variation:',
            color = 'orange',
            verbose = verbose
        )
        if (!is.null(normalization)) {
            data.to.use <- expr.data.nor
        } else data.to.use <- expr.data

        #### regress out variables ####
        if(!is.null(regress.out.variables)){
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(regress.out.variables, collapse = ' & '),
                    ' variables will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose)
            data.to.use <- t(data.to.use)
            lm.formula <- paste('se.obj', regress.out.variables, sep = '$')
            adjusted.data <- lm(as.formula(
                paste('data.to.use', paste0(lm.formula, collapse = '+') , sep = '~'))
                )
            data.to.use <- t(adjusted.data$residuals)
            colnames(data.to.use) <- colnames(se.obj)
            row.names(data.to.use) <- row.names(se.obj)
        }
        ### apply mad within each homogeneous sample groups with respect to the unwanted variable ####
        printColoredMessage(
            message = paste0(
                '- Performing MAD on individual gene expression',
                ' within each homogeneous sample groups with respect to the unwanted variables.'),
            color = 'orange',
            verbose = verbose
            )
        #### find all possible sample groups with respect to the unwanted variables ####
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            clustering.method =  clustering.method,
            nb.clusters = nb.clusters,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            verbose = verbose
            )
        #### apply mad  ####
        selected.homo.uv.groups <- findRepeatingPatterns(
            vec = homo.uv.groups,
            n.repeat = min.sample.for.mad
            )
        if (isTRUE(length(selected.homo.uv.groups) > 0)){
            bio.genes <- sapply(
                selected.homo.uv.groups,
                function(x){
                    index.samples <- homo.uv.groups == x
                    matrixStats::rowMads(x = data.to.use[ , index.samples, drop = FALSE])
                })
            bio.genes <- matrixStats::rowMedians(bio.genes)
            bio.genes <- data.frame(bio.mad = bio.genes)
            set.seed(3322)
            bio.genes$bio.ranks <- rank(x = bio.genes$bio.mad, ties.method = 'random')
            } else if (isTRUE(length(selected.homo.uv.groups) == 0))
            stop(paste0(
                'There is no any homogenous sample groups with at least ',
                min.sample.for.mad,
                ' samples to perform MAD.')
                )
    }

    # Intermediate file ####
    ## read intermediate file ####
    if (isTRUE(use.imf)){
        if(is.null(imf.name)){
            imf.name <- paste0(assay.name, '|un.supervised|', ncg.selection.method)
        }
        if(is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
        bio.genes <- all.tests$bio.genes
        anova.genes.uv <- all.tests$anova.genes.uv
        corr.genes.uv <- all.tests$corr.genes.uv
    }

    ## save intermediate file ####
    if(isTRUE(save.imf)){
        if(length(se.obj@metadata$IMF) == 0 ) {
            se.obj@metadata[['IMF']] <- list()
        }
        if(!'NCG' %in% names(se.obj@metadata[['IMF']])){
            se.obj@metadata[['IMF']][['NCG']] <- list()
        }
        if(is.null(imf.name)){
            imf.name <- paste0(assay.name, '|un.supervised|', ncg.selection.method)
        }
        if(!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
        }
        se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list(
            bio.genes = bio.genes,
            anova.genes.uv = anova.genes.uv,
            corr.genes.uv = corr.genes.uv)
    }

    # Selection of NCG ####
    printColoredMessage(
        message = '-- Summarizing the statistical results to select a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose
        )
    ## product, sum or average of ranks ####
    if (ncg.selection.method %in% c('prod', 'sum', 'average')) {
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        all.uv.ranks <- lapply(
            all.uv.tests,
            function(x){
                temp <- get(x)
                if(length(names(temp))!=0){
                    ranks.data <- lapply(
                        names(temp),
                        function(y) temp[[y]]$ranked.genes)
                    ranks.data <- do.call(cbind, ranks.data)
                    colnames(ranks.data) <- names(temp)
                    ranks.data}
            })
        all.uv.ranks <- do.call(cbind, all.uv.ranks)
        row.names(all.uv.ranks) <- row.names(se.obj)
        all.ranks <- cbind(all.uv.ranks, bio.genes[ , 'bio.ranks', drop = FALSE])
        ### product of ranks ####
        if (ncg.selection.method == 'prod'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowProds(as.matrix(all.ranks))
            if(sum(is.infinite(stat.summary)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## average of ranks ####
        if (ncg.selection.method == 'sum'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowSums(as.matrix(all.ranks))
        }
        ## sum of ranks ####
        if (ncg.selection.method == 'average'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowMeans(as.matrix(all.ranks))
        }
        ## select top genes as NCGS ####
        all.ranks$stat.summary <- stat.summary
        set.seed(112233)
        all.ranks$rank.stat.summary <- rank(x = all.ranks$stat.summary, ties.method = 'random')
        all.ranks <- all.ranks[order(all.ranks$rank.stat.summary, decreasing = FALSE) , ]
        ncg.selected <- row.names(all.ranks[1:round(nb.ncg* nrow(se.obj), digits = 0) , ])
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## non-overlap approach ####
    if (ncg.selection.method == 'non.overlap'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs using the "non.overlap" approach.',
            color = 'orange',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '*1: selecting top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by possible bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
        )
        ### select genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
        printColoredMessage(
            message = paste0(
                '- ',
                length(top.bio.genes),
                ' genes are selected.'),
            color = 'blue',
            verbose = verbose
        )

        ## select genes affected by unwanted variation ####
        printColoredMessage(
            message = paste0(
                '*2: selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by each unwanted variation.'),
            color = 'blue',
            verbose = verbose
        )
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                            row.names(temp.data[[y]])[index] })
                        ))
                }
            })))
        printColoredMessage(
            message = paste0(
                '- ',
                length(top.uv.genes),
                ' genes are selected.'),
            color = 'blue',
            verbose = verbose
        )
        ## select of NCGS ####
        printColoredMessage(
            message = '- all genes found in 1 will be excluded from ones found in 2.',
            color = 'blue',
            verbose = verbose
        )
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## quantile approach ####
    if (isTRUE(ncg.selection.method == 'quantile')){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "quantile" approach..',
            color = 'blue',
            verbose = verbose)
        ### find biological percentile ####
        bio.quan <- quantile(x = bio.genes$bio.mad, probs = bio.percentile)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.mad > bio.quan]

        ## find UV percentile ####
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$statistic < quantile(x = temp.data[[y]]$statistic , probs = uv.percentile)
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        printColoredMessage(
            message = paste0(
                '- Select ', length(top.uv.genes), ' genes with the UV F-statistics higher than ',
                ' (' , uv.percentile* 100, '% percentile), and exclude any genes presents in ',
                length(top.bio.genes),
                ' genes with the biological F-statistics higher than ',
                ' (' , bio.percentile* 100, '% percentile).'),
            color = 'blue',
            verbose = verbose)
        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
    }

    ## auto approach ####
    if (ncg.selection.method == 'auto'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs using the "auto" approach.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of the most affected genes by unwanted variation, and then exclude any genes that are also in the top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### select genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]

        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        ## select NCG ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        printColoredMessage(
            message = paste0('- ', length(ncg.selected), ' genes are found.'),
            color = 'blue',
            verbose = verbose
            )
        ## assess the need for grid search ####
        nb.ncg <- round(c(nb.ncg * nrow(se.obj)), digits = 0)
        ncg.ranges <- round(x = 0.01 *nb.ncg, digits = 0)
        if(length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges)) {
            if(isTRUE(nb.ncg > length(ncg.selected))){
                con <- parse(text = paste0("nb.ncg", ">", "length(ncg.selected)"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ', length(ncg.selected),
                        ' is less than the number (', nb.ncg ,') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            if(isTRUE(nb.ncg < length(ncg.selected))){
                con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ', length(ncg.selected),
                        ' is larger than the number (', nb.ncg ,') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ## grid search ####
            ### grid group: both bio and uv variable ####
            if(grid.group == 'both'){
                printColoredMessage(
                    message = '- The grid search will be applied on both biological and unwanted factors. ',
                    color = 'blue',
                    verbose = verbose
                    )
                #### increasing order ####
                if(grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the values of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                        )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if(top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ### check selection ####
                if(length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                ### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 2)
                if(top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if(top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Update the selection. Select top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
                ncg.selected <- row.names(se.obj) %in% ncg.selected
            }
            ##### grid group: bio ####
            if (grid.group == 'bio'){
                printColoredMessage(
                    message = '- The grid search will be applied on biological factor. ',
                    color = 'blue',
                    verbose = verbose
                    )
                ###### increasing order ####
                if(grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- nrow(se.obj) - top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if(top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### check selection ####
                if(length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 0)
                if(top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Update the selection. Select top ',
                        top.rank.uv.genes * 100,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose
                    )
            }
            ##### grid group: uv ####
            if (grid.group == 'uv'){
                printColoredMessage(
                    message = '- The grid search will be applied on genes that are affected by unwanted variation. ',
                    color = 'blue',
                    verbose = verbose)
                ###### increasing order ####
                if(grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- nrow(se.obj) - top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### check selection ####
                if(length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if(top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Update the selection. Select top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes * 100,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }else {
                printColoredMessage(
                    message = paste0('- ', length(ncg.selected), ' genes are selected as NCGs.'),
                    color = 'blue',
                    verbose = verbose)
            }
        } else ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    printColoredMessage(
        message = paste0(
            'A set of ',
            sum(ncg.selected),
            ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose
        )

    # Plotting ####
    all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
    all.test.res <- lapply(
        all.uv.tests,
        function(x){
            if (!is.null(x)){
                temp.data <- get(x)
                all.test.res <- lapply(
                    names(temp.data),
                    function(y){
                        test.res <- temp.data[[y]][ , 'ranked.genes', drop = FALSE]
                        test.res$group <- rep(y , nrow(test.res))
                        test.res
                    })
                all.test.res <- do.call(cbind, all.test.res)
            }
            all.test.res
        })
    all.test.res <- do.call(cbind, all.test.res)
    temp.data <- lapply(
        seq(1, ncol(all.test.res), 2),
        function(x){
            temp.data <- all.test.res[ , x, drop = FALSE]
            colnames(temp.data) <- all.test.res[ , x+1][1]
            temp.data
        })
    temp.data <- do.call(cbind, temp.data)
    temp.data$bio <- bio.genes$bio.ranks
    temp.data$NCG <- ncg.selected
    ha <- ComplexHeatmap::rowAnnotation(
        NCG = temp.data$NCG,
        col = list(NCG = c('TRUE' = 'gray10', 'FALSE' = 'gray'))
    )
    ncg.plot <- ComplexHeatmap::Heatmap(
        temp.data[ , seq_len(ncol(temp.data) - 1)],
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        right_annotation = ha,
        column_names_rot = 45,
        heatmap_legend_param = list(
            title = 'Ranks',
            title_gp = grid::gpar(fontsize = 14),
            by_row = TRUE,
            ncol = 1)
    )
    if (isTRUE(plot.output)) print(ncg.plot)

    # Performance assessment of the selected NCG ####
    ## pca ####
    if (isTRUE(assess.ncg)){
        printColoredMessage(
            message = '-- Assessing the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Performing PCA on only the selected genes as NCG.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Exploring the association of the first ',
                nb.pcs,
                ' with the ',
                paste0(uv.variables, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose
            )
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u
        if (is.null(variables.to.assess.ncg))
            variables.to.assess.ncg <- uv.variables
        ## regression and vector correlations ####
        all.corr <- lapply(
            variables.to.assess.ncg,
            function(x){
                if (class(se.obj[[x]]) %in% c('numeric', 'integer')){
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared)
                } else if(class(se.obj[[x]]) %in% c('factor', 'character')){
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y){ cca <- cancor(
                            x = pca.data[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)})
                }
            })
        names(all.corr) <- variables.to.assess.ncg
        pcs <- Groups <- NULL
        pca.ncg <- as.data.frame(do.call(cbind, all.corr)) %>%
            dplyr::mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'Groups',
                values_to = 'ls'
                )
        p.assess.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), size = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab ('Correlations') +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
            )
        if (isTRUE(plot.output)) print(p.assess.ncg)
    }
    # Save the NCGs ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Saving the selected NCGs:',
        color = 'magenta',
        verbose = verbose
        )
    if(is.null(ncg.group)){
        ncg.group <- paste0('ncg|unsupervised')
    }
    if(is.null(output.name)){
        output.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|AnoCorrMad:',
            ncg.selection.method,
            '|',
            assay.name
            )
    }

    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Saving the selected set of NCG to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        ## Check
        if (!'NCG' %in% names(se.obj@metadata)) {
            se.obj@metadata[['NCG']] <- list()
        }
        ## check
        if (!'un.supervised' %in% names(se.obj@metadata[['NCG']])) {
            se.obj@metadata[['NCG']][['un.supervised']] <- list()
        }
        ## check
        if (!ncg.group %in% names(se.obj@metadata[['NCG']][['un.supervised']])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]] <- list()
        }
        ## check
        if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['ncg.set']] <- list()
        }
        ## check
        if (!output.name %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['ncg.set']])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['ncg.set']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['ncg.set']][[output.name]] <- ncg.selected

        if (isTRUE(assess.ncg)){
            ## check
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]])) {
                se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['assessment.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['assessment.plot']])) {
                se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['assessment.plot']][[output.name]] <- list()
            }
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group]][['assessment.plot']][[output.name]] <- p.assess.ncg
        }

        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## export output as vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '- The NCGs are outpputed as a logical vector.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose
            )
        return(ncg.selected)
    }
}
