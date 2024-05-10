#' Find a set of negative control genes (NCG) using ANOVA and correlation across all samples.

#' @author Ramyar Molania

#' @description
#' This function employs correlation and ANOVA analyses across all samples to identify a set of genes suitable as negative
#' control genes for RUV-III-PRPS normalization. Correlation analysis is utilized to identify genes highly affected by
#' continuous sources of variation, while ANOVA is used to identify genes affected by categorical sources of variation.
#' The function selects genes as negative control genes (NCG) based on high correlation coefficients and F-statistics for
#' unwanted sources of variation, and low correlation coefficients and F-statistics for biological sources of variation.
#' Various approaches are employed for the final gene selection; please refer to the details for more information.


#' @details
#' The function uses 5 ways to summarize two gene-level F-statistics obtained for the biological and unwanted variation
#' . The function uses either the values or the ranks of F-statistics for NCGs selection. The function ranks the
#' negative of F-statistics values for unwanted variation. The lower the ranks, the greater the impact of unwanted
#' variation on genes. The function ranks the F-statistics for biological variation. The higher the ranks, the greater
#' the impact of biological variation on genes. The options are 'prod', sum', 'average', 'auto' or 'non.overlap' and
#' 'quantile'.
#'
#' If 'prod', 'sum' and 'average' is set:
#'
#' The product, sum or average of ranks of F-statistics is calculated. Then, the function selects 'nb.ncg' 'numbers of
#' genes as negative control genes that have the lowest ranks.
#'
#' If 'non.overlap' is selected:
#' \enumerate{
#'    \item The function selects the top ‘top.rank.bio.genes’ genes that have the highest ranks of F-statistics
#'    for biological variation.
#'    \item The function selects the top ‘top.rank.uv.genes’ genes that have the lowest ranks of F-statistics for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1. This will be a set of genes as
#'    negative control genes.
#' }
#'
#' If 'auto' is selected:
#' \enumerate{
#'    \item The function selects the top ‘top.rank.bio.genes’ genes that have the highest ranks of F-statistics for
#'    biological variation.
#'    \item  The function selects the top ‘top.rank.uv.genes’ genes that have the lowest ranks of F-statistics  for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#'    \item If the number of selected genes is larger or smaller than the specified ‘nb.ncg’, the function applies an
#'    auto search to find approximate ‘nb.ncg’ of genes as negative control genes as follow. The auto search will either
#'    decrease or increase the values of either ‘top.rank.bio.genes’ or ‘top.rank.uv.genes’ or both till to find
#'    approximate ‘nb.ncg’ of genes as negative control genes.
#' }
#' If 'quantile' is selected:
#' \enumerate{
#'    \item The function selects the ‘bio.percentile’ percentile of F-statistics for biological variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function selects the ‘uv.percentile’ percentile of F-statistics for unwanted variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#' }


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object. The selected assay should
#' be the one that will be used for the RUV-III-PRPS normalization. We recommend to use raw data without any transformation.
#' @param nb.ncg Numeric. Specifies the number of genes to be chosen as negative control genes (NCG) when the
#' 'ncg.selection.method' parameter is set to 'auto'. This value ,'nb.ncg', corresponds to a fraction of the total genes
#' in the SummarizedExperiment object. By default, this fraction is set to 0.1.
#' @param bio.variables Symbols. Indicates the columns name(s) that contain biological variable(s) in the SummarizedExperiment
#' object. This cannot be set to 'NULL'.
#' @param uv.variables Symbols. Indicates the column name(s) that contain unwanted variable(s) in the SummarizedExperiment
#' object. This cannot be set to 'NULL'.
#' @param ncg.selection.method Symbol. Indicates how to summarize different statistics and select a set genes as negative
#' control genes. The options are 'prod', 'average', 'sum', 'non.overlap', 'auto' and 'quantile'. The default is set to
#' 'non.overlap'. For more information, refer to the details of the function.
#' @param top.rank.bio.genes Numeric. Indicates the percentage of top ranked genes that are highly affected by the biological
#' variation. This is required to be specified when the 'ncg.selection.method' is either 'auto' or 'non.overlap'. The default
#' is set to 0.2.
#' @param top.rank.uv.genes Numeric. Indicates the percentage of top ranked genes that are highly affected by the unwanted
#' variables. This is required to be specified when the 'ncg.selection.method' is either 'auto' or 'non.overlap'. The default
#' is set to 0.2.
#' @param bio.percentile Numeric. The percentile cut off to select genes that are highly affected with the biological
#' variation. The default is 0.8.
#' @param uv.percentile Numeric. The percentile cut off to select genes that are highly affected with the unwanted
#' variation. The default is 0.8.
#' @param bio.percentile Numeric. The percentile cut off of F-statistics and correlation coefficients to select genes that
#' are highly affected with the biological variation. The default is set to 0.8.
#' @param bio.percentile Numeric. The percentile cut off of F-statistics and correlation coefficients to select genes that
#' are highly affected with the unwanted variation. The default is set to 0.8.
#' @param grid.nb Numeric. Indicates the number of genes for grid search when the 'ncg.selection.method' is auto'. In the
#' 'auto' approach, the grid search starts with the initial top.rank.uv.genes' and 'top.rank.bio.genes' values and the add
#' or drop  the 'grid.nb' in each loop to find 'nb.ncg' of genes as negative control genes.
#' @param grid.group Symbol. Indicates whether the grid search should be performed on biological, unwanted or both factors
#' when the ncg.selection.method' is set to 'auto'. The options are 'bio', 'uv' or 'both'. The default is set to 'uv'.
#' For more details, refer to the details of the function.
#' @param grid.direction Symbol. Indicates the grid search should be performed in decreasing or increasing order, when the
#' 'ncg.selection.method' is set to 'auto'.The options are 'increase' or 'decrease'. The default is set to 'decrease'.
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to perform correlation analyses between
#' continuous sources of variation (biological and unwanted variation) with individual gene expression. The default is set
#' to 3. The minimum value can be 3.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to perform correlation analyses
#' between continuous sources of variation (biological and unwanted variation) with individual gene expression.The default
#' is set to 10. The minimum value can be 3.
#' @param regress.out.bio.variables Symbol. Indicates the column names of biological variables in the SummarizedExperiment
#' object that will be regressed out from the data before performing correlation and ANOVA. Regressing out biological variables
#' might help to better identify genes that are highly affected by unwanted variation. The default is 'NULL'.
#' @param regress.out.uv.variables Symbol. Indicates the column names of unwanted variables the SummarizedExperiment object
#' that will be regressed out from the data before performing correlation and ANOVA. Regressing out unwanted variables
#' might help to better identify genes that are highly affected by biological variation.The default is 'NULL'.
#' @param normalization Symbols.I ndicates winch normalization method to use for library size adjustment before fining genes
#' that are highly affected by biological variation. The default is 'CPM'. Refer to the 'applyOtherNormalization' function
#' for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before any statistical analyses.
#' The default is set to 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before applyinglog transformation.
#' The default is set to 1.
#' @param corr.method Numeric. Indicates which correlation methods should be used for the correlation analyses. The options
#' are 'pearson' or 'spearman'. The default is set to 'spearman'.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation, by default it is set
#' to 0.05. Refer to the function 'correls' from the Rfast R package for more details.
#' @param rho Numeric. The value of the hypothesised correlation to be used in the hypothesis testing. The default it is
#' set to 0. Refer to the function 'correls' from the Rfast R package for more details.
#' @param anova.method Symbols. Indicates which ANOVA methods to use. The options are 'aov' or 'welch'. The default is
#' 'aov'. Refer to the function 'row_oneway_equalvar' or 'row_oneway_welch' from the R package matrixTests for more details.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected genes as negative control or not.
#' This analysis involves principal component analysis on the selected genes and then explore the R^2 or vector correlation
#' between the 'nb.pcs' first principal components and with the biological and unwanted variables. The default is set to 'TRUE'.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected genes as NCG needs to be evaluated. The default is 'NULL'. This
#' means all the variables specified in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected negative control genes to be
#' used to assess the performance of them. The default is set to 5.
#' @param center Logical. Indicates whether to scale the data or not before applying SVD. If center is 'TRUE', then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is 'FALSE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not. If is 'TRUE',
#' the function 'checkSeObj' will be used. The default is 'TRUE'.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If
#' 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and 'uv.variables'
#'  will be excluded. By default, it is set to 'non'.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment
#' object or to output the result. The default is 'TRUE'.
#' @param output.name Symbol. A representation for the output's name. If set to 'NULL', the function will choose a name
#' automatically.
#' @param save.imf Logical. Indicates whether to save the intermediate file or not. If 'TRUE', the function save the results
#' of the statistical analyses in the metadata of the SummarizedExperiment object. Then, if users want to change the parameters
#' including 'nb.ncg', 'ncg.selection.method', top.rank.bio.genes' and 'top.rank.uv.genes', the analyses  will not be re-calculated.
#' The default is set to 'FALSE'.
#' @param use.imf Logical. Indicates whether to use the intermediate file or not. The default is set to 'FALSE'.
#' @param imf.name Symbol. A name to save the intermediate file. If is 'NULL', the function creates a name.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return Either the SummarizedExperiment object containing the a set of negative control genes in the metadata or a
#' logical vector of the selected negative control genes.

#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom matrixTests row_oneway_welch row_oneway_equalvar
#' @importFrom BiocSingular runSVD bsparam
#' @importFrom fastDummies dummy_cols
#' @importFrom matrixStats rowProds
#' @importFrom tidyr pivot_longer
#' @importFrom stats quantile
#' @importFrom Rfast correls
#' @importFrom dplyr mutate
#' @import ggplot2
#' @export

findNcgAcrossSamples <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.2,
        top.rank.uv.genes = 0.2,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        regress.out.bio.variables = NULL,
        regress.out.uv.variables = NULL,
        normalization = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        output.name = NULL,
        save.imf = FALSE,
        imf.name = NULL,
        use.imf = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgAcrossSamples function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####
    if(!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    } else if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty or "NULL".')
    } else if (is.null(uv.variables)){
        stop('The "uv.variables" cannot be empty or "NULL".')
    } else if(!is.vector(bio.variables) | !is.vector(uv.variables) ){
        stop('The "uv.variables" and "bio.variables" must be a vector of variables names in the SummarizedExperiment object.')
    } else if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('Individual specified variable must be either in "bio.variables" or "uv.variables".')
    } else if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" must be a positve value 0 < nb.ncg < 1.')
    } else if (!ncg.selection.method %in% c('prod', 'sum', 'average', 'auto', 'non.overlap', 'quantile')){
        stop('The "ncg.selection.method" muat be one of "prod", "sum", "average", "auto", "non.overlap" or "quantile".')
    } else if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" must be a positve value  0 < top.rank.bio.genes < 1.')
    } else if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    } else if (is.null(min.sample.for.aov)){
        stop('The "min.sample.for.aov" cannot be empty.')
    } else if (min.sample.for.aov <= 2){
        stop('The "min.sample.for.aov" should be at least 3.')
    } else if (is.null(min.sample.for.correlation)){
        stop('The min.sample.for.correlation cannot be empty.')
    } else if (min.sample.for.correlation >= ncol(se.obj) | min.sample.for.correlation < 3){
        stop('The "min.sample.for.correlation" msut be more than 2 and less than the total number of samples in the data.')
    } else if (!anova.method %in% c('aov', 'welch')){
        stop('The "anova.method" must be one of the "aov" or "welch".')
    } else if (isFALSE(is.logical(assess.ncg))){
        stop('The "assess.ncg" must be "TRUE" or "FALSE.')
    } else if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    } else if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(assess.se.obj))) {
        stop('The "assess.se.obj" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(save.se.obj))) {
        stop('The "save.se.obj" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(use.imf))) {
        stop('The "use.imf" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(save.imf))) {
        stop('The "save.imf" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(verbose))) {
        stop('The "verbose" must be "TRUE" or "FALSE".')
    }
    if(!is.null(regress.out.bio.variables) | !is.null(regress.out.uv.variables)){
        if(is.logical(regress.out.bio.variables) | is.logical(regress.out.uv.variables))
            stop(paste0('The "regress.out.bio.variables" or "regress.out.bio.variables" ',
                        'must names of columns in the the SummarizedExperiment object.'))
    }
    if(isTRUE(ncg.selection.method == 'auto')){
        if(isFALSE(is.numeric(grid.nb))){
            stop('The "grid.nb" must be a postive integer value.')
        } else if(grid.nb < 0 | length(grid.nb) > 1 ){
            stop('The "grid.nb" must be a postive integer value.')
        } else if(isTRUE(is.logical(grid.group))){
            stop('The "grid.group" must be on of the "both", "uv" or "bio".')
        } else if(isTRUE(length(grid.group) > 1)){
            stop('The "grid.group" must be on of the "both", "uv" or "bio".')
        } else if (isTRUE(!grid.group %in% c('both', 'uv', 'bio'))){
            stop('The "grid.group" must be on of the "both", "uv" or "bio".')
        } else if(isTRUE(is.logical(grid.direction))){
            stop('The "grid.direction" must be on of the "decrease",or "increase".')
        } else if(isTRUE(length(grid.direction) > 1)){
            stop('The "grid.direction" must be on of the "decrease" or "increase".')
        } else if (isTRUE(!grid.direction %in% c('decrease', 'increase'))){
            stop('The "grid.direction" must be on of the "decrease" or "increase".')
        }
    }
    if(!is.null(normalization)){
        if(!is.null(regress.out.uv.variables))
            printColoredMessage(
                message = paste0('Both normalization and regress.out.uv.variables are selected.',
                'The function will perfom normalization first and the regression the UV variables.'),
                color = 'magenta',
                verbose = verbose)
    }
    if(isTRUE(apply.log)){
        if(length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
            stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if(isTRUE(ncg.selection.method == 'quantile')){
        if(is.null(bio.percentile) | is.null(uv.percentile))
            stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
        if(bio.percentile > 1 | bio.percentile < 0)
            stop('The "bio.percentile" must be a postive value between 0 and 1.')
        if(uv.percentile > 1 | uv.percentile < 0)
            stop('The "uv.percentile" must be a postive value between 0 and 1.')
    }
    if (is.null(assess.se.obj)) {
        if (isTRUE(sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables))) {
            stop('All or some of "bio.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (!is.null(variables.to.assess.ncg)) {
            if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg)))
                stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
        }
    }

    # Check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, variables.to.assess.ncg)),
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
            message = paste0('Applying log2 + ', pseudo.count, ' (pseudo.count) on the ', assay.name, ' data.'),
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
            message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
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
    ## regress out unwanted variables ####
    if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Note, we do not recommend regressing out the ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' if they are largely associated with the ',
                paste0(bio.variables, collapse = ' & '),
                '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.uv <- t(expr.data.nor)
        uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
        expr.data.reg.uv <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
    }
    if(!is.null(regress.out.uv.variables) & is.null(normalization)){
        printColoredMessage(
            message = paste0(
                'The',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Note: we do not recommend regressing out ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                'if they are largely associated with the ',
                paste0(bio.variables, collapse = ' & '), '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.uv <- t(expr.data)
        uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
        expr.data.reg.uv <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~'
        )))
        expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
    }

    ## regress out biological variables ####
    if(!is.null(regress.out.bio.variables)){
        printColoredMessage(
            message = paste0(
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out the ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                'if they are largely associated with the ',
                paste0(uv.variables, collapse = ' & '), '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.bio <- t(expr.data)
        bio.variables.all <- paste('se.obj', regress.out.bio.variables, sep = '$')
        expr.data.reg.bio <- lm(as.formula(paste(
            'expr.data.reg.bio',
            paste0(bio.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data.reg.bio <- t(expr.data.reg.bio$residuals)
        colnames(expr.data.reg.bio) <- colnames(se.obj)
        row.names(expr.data.reg.bio) <- row.names(se.obj)
    }

    # Gene-level ANOVA and correlation analyses  ####
    if(isFALSE(use.imf)){
        ## select genes that are highly affected by unwanted variation ####
        printColoredMessage(
            message = '-- Find genes that are highly affected by each sources of unwnated variation:',
            color = 'magenta',
            verbose = verbose)
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]])))
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
        ### anova between genes and categorical sources of unwanted variation ####
        if(length(categorical.uv) > 0 ){
            if(!is.null(regress.out.bio.variables)){
                data.to.use <- expr.data.reg.bio
            } else data.to.use <- expr.data
            printColoredMessage(
                message = paste0(
                    '- Perform ANOVA between individual gene-level ',
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
                        n.repeat = min.sample.for.aov)
                    if(length(keep.samples) == 0){
                        stop(paste0(
                            'There are not enough samples to perfrom ANOVA between individual gene expression and the ',
                            x,
                            ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                            x,
                            'from the uv.variables and re-run the function.'))
                    } else if(length(keep.samples) == 1 ){
                        stop(paste0(
                            'There is only a single batch from ',
                            x,
                            ' that have enough samples ',
                            min.sample.for.aov,
                            '(min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove'))
                    } else if(length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
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
                    if(anova.method == 'aov'){
                        anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                    } else if(anova.method == 'welch.correction'){
                        anova.gene.batch <- as.data.frame(row_oneway_welch(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                    }
                    set.seed(2233)
                    anova.gene.batch$ranked.genes <- rank(
                        x = -anova.gene.batch[, 'statistic'],
                        ties.method = 'random')
                    anova.gene.batch
                })
            names(anova.genes.uv) <- categorical.uv
            rm(data.to.use)
        } else anova.genes.uv <- NULL
        ### correlation between genes and categorical sources of unwanted variation ####
        if(length(continuous.uv) > 0 ){
            if(!is.null(regress.out.bio.variables)){
                data.to.use <- expr.data.reg.bio
            } else data.to.use <- expr.data
            printColoredMessage(
                message = paste0(
                    '- Perform ',
                    corr.method,
                    ' correlation between individual gene-level ',
                    'expression and each continuous source of unwanted variations: ',
                    paste0(continuous.uv, collapse = '&'),
                    '.'),
                color = 'blue',
                verbose = verbose
            )
            if(ncol(se.obj) <= min.sample.for.correlation){
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
                        x = t(data.to.use),
                        type = corr.method,
                        a = a ,
                        rho = rho))
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 3),
                        corr.genes.var[, 5, drop = FALSE])
                    set.seed(2233)
                    colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                    corr.genes.var$ranked.genes <- rank(-abs(corr.genes.var[, 'statistic']), ties.method = 'random')
                    row.names(corr.genes.var) <- row.names(data.to.use)
                    corr.genes.var
                })
            names(corr.genes.uv) <- continuous.uv
            rm(data.to.use)
        } else corr.genes.uv <- NULL

        ## select genes that are not highly affected by biology ####
        printColoredMessage(
            message = '-- Find genes that are not highly affected by sources of biological variation:',
            color = 'magenta',
            verbose = verbose)
        bio.var.class <- unlist(lapply(
            bio.variables,
            function(x) class(colData(se.obj)[[x]]) ))
        continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
        categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
        ### anova between genes and categorical sources of biological variation ####
        if(length(categorical.bio) > 0 ){
            if(!is.null(normalization) & is.null(regress.out.uv.variables)){
                data.to.use <- expr.data.nor
            } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data.reg.uv
            } else if(!is.null(regress.out.uv.variables) & is.null(normalization)){
                data.to.use <- expr.data.reg.uv
            } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data
            }
            printColoredMessage(
                message = paste0(
                    '- Perform ANOVA between individual gene-level ',
                    'expression and each categorical source of biological variation: ',
                    paste0(categorical.bio, collapse = ' & '),
                    '.'),
                color = 'blue',
                verbose = verbose
            )
            anova.genes.bio <- lapply(
                categorical.bio,
                function(x) {
                    keep.samples <- findRepeatingPatterns(
                        vec = colData(se.obj)[[x]],
                        n.repeat = min.sample.for.aov)
                    if(length(keep.samples) == 0){
                        stop(paste0(
                            'There are not enough samples to perfrom ANOVA between individual genes expression and the ',
                            x,
                            ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                            x,
                            'from the bio.variables and re-run the function.'))
                    } else if(length(keep.samples) == 1 ){
                        stop(paste0(
                            'There is only a single batch from ',
                            x,
                            ' that have enough samples ',
                            min.sample.for.aov,
                            '(min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove',
                            x,
                            'from the bio.variables and re-run the function'))
                    } else if(length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
                        not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                        printColoredMessage(
                            message = paste0(
                                'Note, the',
                                paste0(not.coverd, collapse = '&'),
                                ' groups do not have enough samples for the ANOVA analysis.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                    if(anova.method == 'aov'){
                        anova.genes <- as.data.frame(row_oneway_equalvar(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                    } else if(anova.method == 'welch.correction'){
                        anova.genes <- as.data.frame(row_oneway_equalvar(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                    }
                    set.seed(2233)
                    anova.genes$ranked.genes <- rank(anova.genes[, 'statistic'], ties.method = 'random')
                    anova.genes
                })
            names(anova.genes.bio) <- categorical.bio
            rm(data.to.use)
        } else anova.genes.bio <- NULL
        ### correlation between genes and continuous sources of biological variation ####
        if(length(continuous.bio) > 0 ){
            if(!is.null(normalization) & is.null(regress.out.uv.variables)){
                data.to.use <- expr.data.nor
            } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data.reg.uv
            } else if(!is.null(regress.out.uv.variables) & is.null(normalization)){
                data.to.use <- expr.data.reg.uv
            } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data
            }
            ### gene-batch anova
            printColoredMessage(
                message = paste0(
                    '- Perform ',
                    corr.method,
                    ' correlation between individual gene-level ',
                    'expression and each continuous sources of biological variation: ',
                    paste0(continuous.bio, collapse = '&'),
                    '.'),
                color = 'blue',
                verbose = verbose)
            if(ncol(se.obj) <= min.sample.for.correlation){
                stop(paste0('There are not enough samples (min.sample.for.correlation:',
                            min.sample.for.correlation,
                            ') to perform correlation analysis.',
                            ' A possible soultion in to lower min.sample.for.correlation.'))
            }
            corr.genes.bio <- lapply(
                continuous.bio,
                function(x) {
                    corr.genes.var <- as.data.frame(correls(
                        y = se.obj@colData[, x],
                        x = t(data.to.use),
                        type = corr.method,
                        a = a ,
                        rho = rho))
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 3),
                        corr.genes.var[, 5, drop = FALSE])
                    row.names(corr.genes.var) <- row.names(data.to.use)
                    set.seed(2233)
                    colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                    corr.genes.var$ranked.genes <- rank(
                        x = abs(corr.genes.var[, 'statistic']),
                        ties.method = 'random')
                    corr.genes.var
                })
            names(corr.genes.bio) <- continuous.bio
        } else corr.genes.bio <- NULL
    }
    # Read the intermediate file ####
    if (isTRUE(use.imf)){
        printColoredMessage(
            message = paste0('- Retrieve the results of ANOVA and correlations from the the SummarizedExperiment object .'),
            color = 'blue',
            verbose = verbose)
        if(is.null(imf.name)){
            imf.name <- paste0(assay.name, '|AcrossSamples|', ncg.selection.method)
        }
        if(is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
        anova.genes.bio <- all.tests$anova.genes.bio
        corr.genes.bio <- all.tests$corr.genes.bio
        anova.genes.uv <- all.tests$anova.genes.uv
        corr.genes.uv <- all.tests$corr.genes.uv
    }
    # Save the intermediate file ####
    if(isTRUE(save.imf)){
        printColoredMessage(
            message = '-- Save a intermediate file:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = '- The results of ANOVA and correlations are saved in the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        if(length(se.obj@metadata$IMF) == 0 ) {
            se.obj@metadata[['IMF']] <- list()
        }
        if(!'NCG' %in% names(se.obj@metadata[['IMF']])){
            se.obj@metadata[['IMF']][['NCG']] <- list()
        }
        if(is.null(imf.name)){
            imf.name <- paste0(assay.name, '|AcrossSamples|', ncg.selection.method)
        }
        if(!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
        }
        se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list(
            anova.genes.bio = anova.genes.bio,
            corr.genes.bio = corr.genes.bio,
            anova.genes.uv = anova.genes.uv,
            corr.genes.uv = corr.genes.uv)
    }

    # Selection of NCG ####
    printColoredMessage(
        message = '-- Selection of a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose)
    ## prod, sum or average of ranks ####
    if(ncg.selection.method %in% c('prod', 'sum', 'average')){
        all.tests <- c(
            'anova.genes.bio',
            'corr.genes.bio',
            'anova.genes.uv',
            'corr.genes.uv')
        all.stats <- lapply(
            all.tests,
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
        all.stats <- do.call(cbind, all.stats)
        row.names(all.stats) <- row.names(se.obj)
        ### product of ranks ####
        if(ncg.selection.method == 'prod'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.summary <- rowProds(all.stats)
            if(sum(is.infinite(stat.summary)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## average of ranks ####
        if(ncg.selection.method == 'sum'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.summary <- rowSums(all.stats)
        }
        ## sum of ranks ####
        if(ncg.selection.method == 'average'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.summary <- rowMeans(all.stats)
        }
        ## select top genes as NCGS ####
        all.stats <- as.data.frame(all.stats)
        row.names(all.stats) <- row.names(se.obj)
        all.stats$stat.summary <- stat.summary
        set.seed(112233)
        all.stats$rank.stat.summary <- rank(x = all.stats$stat.summary, ties.method = 'random')
        all.stats <- all.stats[order(all.stats$rank.stat.summary, decreasing = FALSE) , ]
        ncg.selected <- row.names(all.stats[1:round(nb.ncg* nrow(se.obj), digits = 0) , ])
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
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))

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
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$statistic > quantile(x = temp.data[[y]]$statistic , probs = bio.percentile)
                            row.names(temp.data[[y]])[index] })))
                }
            })))

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
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(c(top.rank.uv.genes * nrow(se.obj)), digits = 0)
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
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
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
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
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
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
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
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x){
                                if(!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
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
                }
            } else {
            printColoredMessage(
                message = paste0('- ', length(ncg.selected), ' genes are selected as NCGs.'),
                color = 'blue',
                verbose = verbose)
        }
    }

    printColoredMessage(
        message = paste0('A set of ', sum(ncg.selected), ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose)

    # Assessment of the selected set of NCG ####
    ### pca on the NCG ####
    if(assess.ncg){
        printColoredMessage(
            message = '-- Assess the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = paste0('- Perform PCA on only selected genes as NCG.'),
            color = 'blue',
            verbose = verbose)
        if(is.null(variables.to.assess.ncg))
            variables.to.assess.ncg <- c(bio.variables, uv.variables)
        printColoredMessage(
            message = paste0(
                '- Explore the association of the first ',
                nb.pcs,
                '  PCs with the ',
                paste0(variables.to.assess.ncg, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u
        ### regression and vector correlations ####
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
            ylab (expression("Correlations")) +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
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
    # Save results ####
    printColoredMessage(
        message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
        color = 'magenta',
        verbose = verbose)
    if(is.null(output.name)){
        output.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|AnoCorrAs:',
            ncg.selection.method,
            '|',
            assay.name)
    }
    ### add results to the SummarizedExperiment object ####
    if(isTRUE(save.se.obj)){
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if(!'supervised' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['supervised']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['NCG']][['supervised']])){
            se.obj@metadata[['NCG']][['supervised']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['supervised']][[output.name]] <- ncg.selected
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgAcrossSamples function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    }
    ### export results as logical vector ####
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The NCGs are outputed as a logical vector.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgAcrossSamples function finished.',
            color = 'white',
            verbose = verbose)
        return(ncg.selected)
    }
}
