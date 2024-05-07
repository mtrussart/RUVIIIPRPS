#' Find a set of negative control genes (NCG) using unsupervised approach.

#' @author Ramyar Molania

#' @description
#' This function identifies a set of genes as negative control genes in cases where no biological variations are known.
#' The function uses the gene-level correlation, ANOVA and MAD  analyses across and between sample groups to find a set
#' of genes as negative control genes for the RUV-III-PRPS normalization.

#' @details
#' The function initially utilizes gene-level correlation and ANOVA to identify genes significantly influenced by continuous
#' and categorical sources of variation, respectively. Subsequently, it conducts Median Absolute Deviation (MAD) analysis
#' on each gene within homogeneous sample groups, considering unwanted variables, to pinpoint genes highly variable due
#' to biological factors. Lastly, various methods are employed to consolidate the statistical findings, ultimately determining
#' an appropriate set of genes as negative control genes for RUV-III-PRPS normalization.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object. The selected assay should
#' be the one that will be used for the RUV-III-PRPS normalization.
#' @param uv.variables Symbol. A symbol or symbols indicating name(s) of the columns in the sample annotation that contains
#' unwanted variables in the SummarizedExperiment object. If there are unknown, the 'identifyUnknownUV' function can
#' estimate them.
#' @param nb.ncg Numeric. Specifies the number of genes to be chosen as negative control genes (NCG) when the
#' ncg.selection.method parameter is set to 'auto'. This value corresponds to a fraction, 'nb.ncg', of the total genes
#' in the SummarizedExperiment object. By default, this fraction is set to 0.1.
#' @param ncg.selection.method Symbol. Indicates how to select a set genes as negative control genes (NCG). For individual
#' genes F-statistics (from ANOVA analysis) , correlation coefficient and Median Absolute Deviation (MAD) is calculated.
#' Ideal control genes have should have possible high F-statistics and correlation coefficient,  and possible lowest MAD
#' values. To summarize these statistic there two options: 'non-overlap' and 'auto'. If 'non-overlap' is selected, the possible
#' non-overlap genes between 'top.rank.uv.genes' and 'top.rank.uv.genes' genes will be selected as NCG. If 'auto' selected,
#' the function finds possible 'nb.ncg' genes as NCGs by changes the 'top.rank.uv.genes' and 'top.rank.uv.genes' automatically.
#' @param top.rank.bio.genes Numeric. Indicates the fraction of top ranked genes that are highly affected by the biological
#' variation. These genes have highest MAD values. The default is 0.5.
#' @param top.rank.uv.genes Numeric. Indicates the fraction of top ranked genes that are highly affected by the unwanted
#' variation variables. These genes have highest F-statistics and correlation coefficients. The default is 0.5.
#' @param uv.percentile Numeric. The percentile cut off to select genes that are highly affected with the unwanted
#' variation. The default is 0.8.
#' @param bio.percentile Numeric. The percentile cut off to select genes that are highly affected with the unwanted
#' variation. The default is 0.8.
#' @param grid.nb Numeric. Indicates the number of genes for grid search when the 'ncg.selection.method' is auto'. The
#' default is 20.
#' @param grid.group Symbol. Indicates whether the grid search should be performed on biological, unwanted or both factors.
#' The options are 'bio', 'uv' or 'both'.
#' @param grid.direction Symbol. Indicates the grid search should be performed in decreasing or increasing order. The
#' option are 'increase' or 'decrease'.
#' @param normalization Symbols. Indicates which normalization method to use for library size adjustment before fining
#' genes that are highly affected by biological variation. The default is CPM. Refer to the 'applyOtherNormalization'
#' function for further details.
#' @param regress.out.variables Symbols. TTTT
#' @param min.sample.for.mad Numeric. Indicates the minimum number of samples to perform Median Absolute Deviation on
#' each gene within homogeneous sample groups, considering unwanted variables. The default is 3.
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to perform ANOVA analyses between categorical
#' sources of variation with individual gene expression. The default is 3.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to perform correlation analyses
#' between continuous sources of unwanted variation with individual gene expression. The default is 10.
#' @param corr.method Numeric. Indicates which correlation methods should be used for the correlation analyses. The options
#' are pearson' or 'spearman'. The default is 'spearman'.
#' @param a Numeric. the significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param anova.method Symbols. Indicates which ANOVA methods to use. The options are 'aov' or 'welch'. The default in
#' 'aov'. Refer to the function 'computeGenesVariableAnova' for more details.
#' @param clustering.method Symbol. A symbol specifying the clustering method to be applied for grouping each continuous
#' source of unwanted variables. Options include 'kmeans', 'cut', and 'quantile'. The default is 'kmeans' clustering.
#' @param nb.clusters Numeric. A value indicating the number of groups for continuous sources of unwanted variation.
#' The default is 3. This implies that each continuous source will be split into 3 groups using the specified
#' 'clustering.method' method.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not. This involves principal
#' component analysis on the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs' first principal
#' components and with biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected genes as NCG needs to be evaluated. The default is NULL. This
#' means all the variables in the 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the
#' performance of NCGs. The defaut is 5.
#' @param center Logical. Indicates whether to scale the data or not before applying SVD. If 'center' is TRUE, then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param assess.variable Logical.
#' @param cat.cor.coef Vector. Containing two numerical values that specifies the cut-off for the correlation coefficient
#' between each pair of categorical variables. The first value applies to correlations between each pair of 'uv.variables',
#' and the second value applies to correlations between each pair of 'bio.variables'. Correlation is computed using the
#' ContCoef function from the DescTools R package. If the correlation between a pair of variables exceeds the cut-off,
#' only the variable with the highest number of factors will be retained, and the other will be excluded from further
#' analysis. By default, both values are set to 0.9.
#' @param cont.cor.coef Vector. Containing two numerical values that specifies the cut-off for the correlation coefficient
#' between each pair of continuous variables. The first value applies to correlations between each pair of continuous
#' 'uv.variables', and the second value applies to correlations between each pair of continuous 'bio.variables'.
#' Correlation is computed using the ContCoef function from the DescTools R package. If the correlation between a pair of
#' variables exceeds the cut-off, only the variable with the highest variance will be retained, and the other will be
#' excluded from further analysis. By default, both values are set to 0.9.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment object or
#' to output the result. The default is TRUE.
#' @param output.name Symbol. Indicates the name of the output in the meta data of the SummarizedExperiment object. The
#' default is 'NULL'. This means the function will create a name based the specified argument.
#' @param use.imf Logical. Indicates whether to use the intermediate file or not. The default is set to 'FALSE'.
#' @param save.imf Logical. Indicates whether to save the intermediate file or not. If 'TRUE', the function save the results
#' of two-way ANOVA. Then, if users want to change the parameters including 'nb.ncg', 'ncg.selection.method',
#' 'top.rank.bio.genes' and 'top.rank.uv.genes', the two-way ANOVA will not be re-calculated.
#' @param imf.name Symbol. A name to save the intermediate file. If is 'NULL', the function creates a name.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return Either the SummarizedExperiment object containing the a set of negative control genes or a logical vector of
#' the selected negative control genes
#'
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
        grid.nb = 20,
        grid.group = 'uv',
        grid.direction = 'decrease',
        normalization = 'CPM',
        regress.out.variables = NULL,
        min.sample.for.mad = 3,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        clustering.method = 'kmeans',
        nb.clusters = 3,
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        scale = FALSE,
        center = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        assess.variable = FALSE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        save.se.obj = TRUE,
        remove.na = 'both',
        output.name = NULL,
        use.imf = FALSE,
        save.imf = FALSE,
        imf.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgsUnSupervised function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check functions inputs ####
    if(length(assay.name) > 1 | is.logical(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    } else if(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
        stop('Some or all the "uv.variables" cannot be found in the SummarizedExperiment object.')
    } else if(nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
    } else if (!ncg.selection.method %in% c('auto', 'non.overlap')){
        stop('The "ncg.selection.method" must be one of "auto" or "non.overlap".')
    } else if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" msut be a positve value  0 < top.rank.bio.genes < 1.')
    } else if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    } else if( grid.nb < 1 | grid.nb > nrow(se.obj)){
        stop(paste0('The "grid.nb" must be a positve value  0 < grid.nb < ', nrow(se.obj), '.'))
    } else if(!grid.group %in% c('bio', 'uv', 'both')){
        stop('The "grid.group" must be one of "bio", "uv" or "non.overlap".')
    } else if(!grid.direction %in% c('increase', 'decrease', 'auto')){
        stop('The "grid.direction" must be one of "increase", "decrease" or "auto".')
    } else if (is.null(min.sample.for.aov)){
        stop('The "min.sample.for.aov" cannot be empty.')
    } else if (min.sample.for.aov <= 2){
        stop('The "min.sample.for.aov" should be at least 3.')
    } else if (is.null(min.sample.for.correlation)){
        stop('The min.sample.for.correlation cannot be empty.')
    } else if (min.sample.for.correlation >= ncol(se.obj) | min.sample.for.correlation < 3){
        stop('The "min.sample.for.correlation" msut be more than 2 and less than the total number of samples in the data.')
    } else if (!anova.method %in% c('aov', 'welch')){
        stop('The anova.method must be one of the "aov" or "welch".')
    } else if (isFALSE(is.logical(assess.ncg))){
        stop('The "assess.ncg" must be "TRUE" or "FALSE.')
    } else if (length(nb.pcs) > 1){
        stop('The "nb.pcs" must be a postive integer value.')
    } else if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    } else if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE.')
    } else if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE.')
    } else if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE.')
    } else if(length(pseudo.count) > 1){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    } else if(pseudo.count < 0){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    } else if (isFALSE(is.logical(assess.se.obj))) {
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
    if(!is.null(regress.out.variables)){
        if (isTRUE(sum(regress.out.variables %in% colnames(colData(se.obj))) != length(regress.out.variables))) {
            stop('All or some of "regress.out.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if(isTRUE(ncg.selection.method == 'quantile')){
        if(is.null(bio.percentile) | is.null(uv.percentile))
            stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
        if(bio.percentile > 1 | bio.percentile < 0)
            stop('The "bio.percentile" must be a postive value between 0 and 1.')
        if(uv.percentile > 1 | uv.percentile < 0)
            stop('The "uv.percentile" must be a postive value between 0 and 1.')
    }

    # Check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(uv.variables, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose)
    }
    # Data transformation and normalization ####
    printColoredMessage(
        message = '-- Data transformation and normalization:',
        color = 'magenta',
        verbose = verbose)
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0('Applying log2 + ', pseudo.count,' (pseudo.count) on the ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0('Applying log2 on the ', assay.name,' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without any transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(x = se.obj, i = assay.name)
    }

    ## normalization ####
    if(!is.null(normalization)){
        expr.data.nor <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
    }

    # Finding negative control genes ####
    if(isFALSE(use.imf)){
        printColoredMessage(
            message = '-- Find negative control genes',
            color = 'magenta',
            verbose = verbose)
        ## find genes that are highly affected by unwanted variation ####
        printColoredMessage(
            message = '-- Find genes that are highly affected by each sources of unwnated variation:',
            color = 'blue',
            verbose = verbose)
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]])))
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
        ### anova between genes and categorical sources of unwanted variation ####
        if(isTRUE(length(categorical.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Perform ANOVA between individual gene-level ',
                    'expression and each categorical source of unwanted variation: ',
                    paste0(categorical.uv, collapse = ' & '), '.'),
                color = 'blue',
                verbose = verbose)
            anova.genes.uv <- lapply(
                categorical.uv,
                function(x) {
                    keep.samples <- findRepeatingPatterns(
                        vec = colData(se.obj)[[x]],
                        n.repeat = min.sample.for.aov)
                    if(isTRUE(length(keep.samples) == 0)){
                        stop(paste0(
                            'There are not enough samples to perfrom ANOVA between individual gene expression and the ',
                            x,
                            ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                            x,
                            ' from the uv.variables and re-run the function.'))
                    } else if(isTRUE(length(keep.samples) == 1)){
                        stop(paste0(
                            'There is only a single batch from in the ',
                            x,
                            ' variable that have enough samples ',
                            min.sample.for.aov,
                            ' (min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove',
                            x,
                            ' from the uv.variables and re-run the function.'
                        ))
                    } else if(isTRUE(length(keep.samples) != length(unique(colData(se.obj)[[x]])))){
                        not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                        printColoredMessage(
                            message = paste0('Note, the ',
                                             paste0(not.coverd, collapse = '&'),
                                             ' batches do not have enough samples for the ANOVA analysis.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                    if(anova.method == 'aov'){
                        anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples]))
                    } else if(anova.method == 'welch'){
                        anova.gene.batch <- as.data.frame(row_oneway_welch(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples]))
                    }
                    set.seed(2233)
                    anova.gene.batch$ranked.genes <- rank(
                        -anova.gene.batch[, 'statistic'],
                        ties.method = 'random')
                    anova.gene.batch
                })
            names(anova.genes.uv) <- categorical.uv
        } else anova.genes.uv <- NULL
        ### correlation between genes and continuous sources of unwanted variation ####
        if(isTRUE(length(continuous.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Perform ',
                    corr.method,
                    ' correlation between individual gene-level ',
                    'expression and each continuous source of unwanted variations: ',
                    paste0(continuous.uv, collapse = '&'), '.'),
                color = 'blue',
                verbose = verbose)
            if(isTRUE(ncol(se.obj) <= min.sample.for.correlation)){
                stop(paste0('There are not enough samples (min.sample.for.correlation:',
                            min.sample.for.correlation,
                            ') to perform correlation analysis.',
                            ' A possible soultion in to lower min.sample.for.correlation.'))
            }
            corr.genes.uv <- lapply(
                continuous.uv,
                function(x) {
                    corr.genes.var <- as.data.frame(correls(
                        y = se.obj@colData[, x],
                        x = t(expr.data),
                        type = corr.method,
                        a = a ,
                        rho = rho))
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 3),
                        corr.genes.var[, 5, drop = FALSE])
                    set.seed(2233)
                    colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                    corr.genes.var$ranked.genes <- rank(-abs(corr.genes.var[, 'statistic']), ties.method = 'random')
                    row.names(corr.genes.var) <- row.names(expr.data)
                    corr.genes.var
                })
            names(corr.genes.uv) <- continuous.uv
        } else corr.genes.uv <- NULL

        ## find genes that are not highly affected by biology ####
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
            adjusted.data <- lm(as.formula(paste('data.to.use', paste0(lm.formula, collapse = '+') , sep = '~')))
            data.to.use <- t(adjusted.data$residuals)
            colnames(data.to.use) <- colnames(se.obj)
            row.names(data.to.use) <- row.names(se.obj)
        }
        ### apply mad within each homogeneous sample groups with respect to the unwanted variable ####
        printColoredMessage(
            message = paste0(
                '- Perform MAD on individual gene expression',
                ' within each homogeneous sample groups with respect to the unwanted variables.'),
            color = 'blue',
            verbose = verbose)
        #### find all possible sample groups with respect to the unwanted variables ####
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            clustering.method =  clustering.method,
            nb.clusters = nb.clusters,
            assess.variables = assess.variable,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            verbose = verbose)
        #### apply mad  ####
        homo.uv.groups <- findRepeatingPatterns(
            vec = homo.uv.groups,
            n.repeat = min.sample.for.mad)
        groups <- unique(homo.uv.groups)
        if(isTRUE(length(groups) > 0)){
            bio.genes <- sapply(
                groups,
                function(x){
                    index.samples <- homo.uv.groups == x
                    matrixStats::rowMads(x = data.to.use[ , index.samples, drop = FALSE])
                })
            bio.genes <- matrixStats::rowMaxs(bio.genes)
            bio.genes <- data.frame(bio.mad = bio.genes)
            set.seed(3322)
            bio.genes$bio.ranks <- rank(x = bio.genes$bio.mad, ties.method = 'random')
        } else if (isTRUE(length(groups) == 0))
            stop('There is no any homogenous sample groups with enough samples to perform MAD.')
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
    # if(!is.null(anova.genes.uv)){
    #     anova.genes.uv.rank <- lapply(
    #         names(anova.genes.uv),
    #         function(x){
    #             temp <- anova.genes.uv[[x]][ , c('ranked.genes'), drop = FALSE]
    #             temp$ranked.genes1 <- rank(x = -temp$ranked.genes, ties.method = 'random')
    #             colnames(temp) <- x
    #             temp
    #         })
    #     anova.genes.uv.rank <- do.call(cbind, anova.genes.uv.rank)
    # } else anova.genes.uv <- NULL
    # if(!is.null(corr.genes.uv)){
    #     corr.genes.uv.rank <- lapply(
    #         names(corr.genes.uv),
    #         function(x){
    #             temp <- corr.genes.uv[[x]][ , c('ranked.genes'), drop = FALSE]
    #             temp$ranked.genes <- rank(x = -temp$ranked.genes, ties.method = 'random')
    #             colnames(temp) <- x
    #             temp
    #         })
    #     corr.genes.uv.rank <- do.call(cbind, corr.genes.uv.rank)
    # } else corr.genes.uv.rank <- NULL
    # all.ranks <- cbind(bio.genes.rank, corr.genes.uv.rank, anova.genes.uv.rank)
    # all.ranks <- all.ranks[ , c('biology', uv.variables )]
    # all.ranks <- all.ranks[order(all.ranks[['biology']]) , ]
    # p <- ComplexHeatmap::Heatmap(
    #     matrix = all.ranks,
    #     cluster_rows = FALSE,
    #     cluster_columns = FALSE,
    #     show_row_names = FALSE,
    #     col = c('grey90', 'grey50', 'cyan', 'darkgreen', 'orange', 'red')
    #     )

    # Selection of NCG ####
    printColoredMessage(message = '-- Selection a set of genes as NCG:',
                        color = 'magenta',
                        verbose = verbose)
    ## product, sum or average of ranks ####
    if (ncg.selection.method %in% c('prod', 'sum', 'average')) {
        all.uv.tests <- c(
            'anova.genes.uv',
            'corr.genes.uv')
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
        if(ncg.selection.method == 'prod'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.summary <- rowProds(as.matrix(all.ranks))
            if(sum(is.infinite(stat.summary)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## average of ranks ####
        if(ncg.selection.method == 'sum'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.summary <- rowSums(as.matrix(all.ranks))
        }
        ## sum of ranks ####
        if(ncg.selection.method == 'average'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
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

    ## non.overlap approach ####
    if (ncg.selection.method == 'non.overlap'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                '- Select top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude all top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose)
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
        ## select of NCGS ####
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
            message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                '- Select top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose)
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
            verbose = verbose)

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
                    verbose = verbose)
                #### increasing order ####
                if(grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the values of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(nrow(se.obj) - top.rank.uv.genes.nb,
                              top.rank.bio.genes.nb)
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
                    verbose = verbose)
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
                        verbose = verbose)
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
                    verbose = verbose)
            }
            ##### grid group: uv ####
            if (grid.group == 'uv'){
                printColoredMessage(
                    message = '- The grid search will be applied on unwanted factor. ',
                    color = 'blue',
                    verbose = verbose)
                ###### increasing order ####
                if(grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose)
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
                        verbose = verbose)
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
        }
    }

    printColoredMessage(
        message = paste0('A set of ', sum(ncg.selected), ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose)

    # Performance assessment of the selected NCG ####
    ## pca ####
    if(assess.ncg){
        printColoredMessage(
            message = '-- Assess the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = '- Perform PCA on only the selected genes as NCG.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                '- Explore the association of the first ',
                nb.pcs,
                '  with the ',
                paste0(uv.variables, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u
        if(is.null(variables.to.assess.ncg)) variables.to.assess.ncg <- uv.variables
        ## regression and vector correlations ####
        all.corr <- lapply(
            variables.to.assess.ncg,
            function(x){
                if(class(se.obj[[x]]) %in% c('numeric', 'integer')){
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
                values_to = 'ls')
        pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
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
        if(verbose) print(pca.ncg)
    }
    # Save the NCGs ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save the selected NCGs:',
        color = 'magenta',
        verbose = verbose)
    if(is.null(output.name)){
        output.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|AnoCorrMad:',
            ncg.selection.method,
            '|',
            assay.name)
    }

    if(isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Save the selected set of NCG to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if (!'NCG' %in% names(se.obj@metadata)) {
            se.obj@metadata[['NCG']] <- list()
        }
        ## check if supervised already exists in the PRPS slot
        if (!'un.supervised' %in% names(se.obj@metadata[['NCG']])) {
            se.obj@metadata[['NCG']][['un.supervised']] <- list()
        }
        ## check if prps.set.name already exists in the PRPS$supervised slot
        if (!output.name %in% names(se.obj@metadata[['NCG']][['un.supervised']])) {
            se.obj@metadata[['NCG']][['un.supervised']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['un.supervised']][[output.name]] <- ncg.selected
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)

        ## export output as vector ####
    }
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = '- The NCGs are outpputed as a logical vector.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose)
        return(ncg.selected)
    }
}




