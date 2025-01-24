#' Finds negative control genes using ANOVA and correlation per sample groups.

#' @author Ramyar Molania

#' @description
#' This function employs correlation and ANOVA analyses across within groups of samples that are homogeneous with respect
#' to biological or unwanted variation. Correlation analysis is utilized to identify genes highly affected by continuous
#' sources of variation, while ANOVA is used to identify genes affected by categorical sources of variation.
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
#'
#' The data pre-processing for the functions findNcgAcrossSamples() and findNcgPerBiologyPerBatch().
#' The input data for correlation and ANOVA analyses can be one or a combination of the following options:
#' \enumerate{
#'    \item item $log_2$ (raw data + pseudo-count).
#'    \item Apply some library size normalizations (if the input data is gene count expression).
#'    \item This will facilitate finding genes that are highly affected by biological variation. item Regress out
#'    biological variation.
#'    \item This will facilitate finding genes that are highly affected by unwanted variation.
#'    item Regress out unwanted variation.
#' }


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name A character string. The name of the assay in the SummarizedExperiment object to be used for selecting
#' negative control genes.
#' @param bio.variables A character vector. The column names that contain biological variables in the SummarizedExperiment
#' object. Both continuous and categorical variables can be provided.
#' @param uv.variables A character vector. The column names that contain unwanted variation (UV) variables in the
#' SummarizedExperiment object.
#' @param nb.ncg Numeric. The proportion of the total genes to be selected as a Negative Control Gene (NCG) set. The
#' default is .1.
#' @param ncg.selection.method A character string. Indicates how to select the genes as NCGs. For individual genes,
#' the two-way ANOVA calculates F-statistics for biological and unwanted variation factors separately. An ideal NCG set
#' should have high F-statistics for the unwanted variation variables and low F-statistics for the biological variables.
#' The function ranks the F-statistics obtained for the biological variables and the negative of the F-statistics obtained
#' for the unwanted variables. The function then offers 5 ways to summarize the ranks of the two F-statistics:
#' 'Prod' (product of the ranks), 'Sum' (sum of the ranks), 'Average' (average of the ranks), 'AbsNoneOverlap' (non-overlapping genes
#' between 'top.rank.uv.genes' and 'top.rank.bio.genes'), and 'NoneOverlap' (non-overlapping genes between 'top.rank.uv.genes'
#' and at least 'top.rank.bio.genes'). The F-statistics for biological and UV are first ranked.
#' @param grid.nb Numeric. The percentage for grid search when the 'ncg.selection.method' is 'NoneOverlap'. In the 'NoneOverlap'
#' approach, the grid search starts with the initial `top.rank.uv.genes` value and adds `grid.nb` in each loop to find
#' the desired 'nb.ncg'.
#' @param top.rank.bio.genes Numeric. The percentage of top-ranked genes that are highly affected by the biological variation.
#' This parameter is required when `ncg.selection.method` is either 'NoneOverlap' or 'AbsNoneOverlap'.
#' @param top.rank.uv.genes Numeric. The percentage of top-ranked genes that are highly affected by the unwanted variation variables.
#' This parameter is required when `ncg.selection.method` is either 'NoneOverlap' or 'Auto'.
#' @param bio.percentile Numeric. The percentile cut-off for selecting genes that are highly affected by the biological
#' variation. The default is set to 0.8.
#' @param uv.percentile Numeric. The percentile cut-off for selecting genes that are highly affected by the unwanted
#' variation. The default is set to 0.8.
#' @param bio.groups A character vector or a symbol. Indicates the column names that contain biological variables in the
#' SummarizedExperiment object.
#' If not NULL, bio.groups' will be used to group samples into different homogeneous biological groups.
#' @param grid.group A character string. Indicates whether the grid search should be performed on biological ('top.rank.bio.genes'),
#' unwanted ('top.rank.uv.genes') or both factors. The options are 'bio', 'uv', or 'both'. The default is set to 'uv'.
#' @param grid.direction A character string. Indicates the direction of the grid search. The options are 'increase' or 'decrease'.
#' The default is set to 'decrease'.
#' @param bio.clustering.method A character string. Indicates which clustering method should be used to group continuous
#' sources of biological variation if any are provided. The default is 'kmeans' clustering.
#' @param nb.bio.clusters Numeric. Indicates the number of clusters for each continuous source of biological variation.
#' The default is set to 3.
#' @param uv.groups A character vector or a symbol. Indicates the column names that contain unwanted variation variables in the
#' SummarizedExperiment object. If not NULL, `uv.groups` will be used to group samples into different homogeneous unwanted
#' variation groups.
#' @param uv.clustering.method A character string. Indicates which clustering method should be used to group continuous
#' sources of unwanted variation. The default is set to 'kmeans' clustering.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous source of unwanted variation (UV).
#' By default, it is set to 2.
#' @param normalization A character string. Indicates which normalization method should be applied to the data before
#' identifying genes that are affected by biological variation. The default is 'CPM'. Refer to the 'applyOtherNormalizations
#' function for more details.
#' @param regress.out.bio.variables A character vector. Indicates the column names of biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before identifying genes that are
#' highly affected by unwanted variation. The default is NULL.
#' @param regress.out.uv.variables A character vector. Indicates the column names of unwanted variation variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before identifying genes that are
#' highly affected by biological variation. The default is NULL.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing correlation and
#' ANOVA. The default is TRUE.
#' @param pseudo.count Numeric. A value to be added as a pseudo-count to all measurements before log-transformation.
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to be present in each group before applying
#' the ANOVA. The default is set to 3.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to be considered in each group before
#' applying correlation analysis. The default is set to 10.
#' @param corr.method A character string. Indicates which correlation method should be used to compute associations between
#' gene-level expression and a continuous variable. The default is set to 'spearman'.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation analysis. The default
#' is set to 0.05.
#' @param rho Numeric. The hypothesized correlation value to be used in the hypothesis testing. The default is set to 0.
#' @param anova.method A character string. Indicates which ANOVA method should be used to compute associations between
#' gene-level expression and a categorical variable. The default is 'aov'.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCGs. This analysis involves
#' principal component analysis (PCA) on the selected NCGs, followed by exploring the R-squared or vector correlation between
#' the first 'nb.pcs' principal components and biological and unwanted variation variables.
#' @param variables.to.assess.ncg A character vector. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected NCGs needs to be evaluated. The default is NULL, meaning all
#' variables in 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components of selected NCGs to be used to assess
#' the performance of NCGs.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis. The default
#'  is set to 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data before applying principal component analysis. The default is
#' set to 'FALSE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is set to  'TRUE'.
#' @param remove.na A character string. Indicates whether to remove NA or missing values from either the 'assays', 'sample.annotation',
#' 'both', or 'none'. If 'assays' is selected, genes containing NA or missing values will be excluded. If 'sample.annotation'
#' is selected, samples with NA or missing values for any `bio.variables` and `uv.variables` will be excluded. The default
#' is set to  'both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object
#' ('se.obj') or to output the result.
#' The default is TRUE.
#' @param output.name A character string. The name for the output. If set to NULL, the function will automatically generate a name.
#' @param ncg.group A character string. Indicates the name of the group of NCGs.
#' @param plot.output Logical. If TRUE, a plot of the NCG selection process will be generated.
#' @param save.imf Logical. Indicates whether to save the intermediate file in the SummarizedExperiment object. If set
#' to TRUE, the results of the two-way ANOVA
#' will be saved. If users wish to adjust parameters such as `nb.ncg`, `ncg.selection.method`, `top.rank.bio.genes`,
#' and `top.rank.uv.genes`, the ANOVA
#' will not be recalculated. This accelerates parameter tuning for NCG selection. The default is FALSE.
#' @param imf.name A character string. The name to use when saving the intermediate file. If set to NULL, the function
#' will automatically generate a name.
#' The default is NULL.
#' @param use.imf Logical. Indicates whether to use the intermediate file for subsequent steps. The default is FALSE.
#' @param verbose Logical. If TRUE, process messages will be displayed.
#'
#' @return Either the SummarizedExperiment object containing a set of negative control genes, or a logical vector of
#' the selected negative control genes.

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

findNcgPerBiologyPerBatch <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        regress.out.bio.variables = NULL,
        regress.out.uv.variables = NULL,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
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
        ncg.group = NULL,
        plot.output = TRUE,
        save.imf = FALSE,
        imf.name = NULL,
        use.imf = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgPerBiologyPerBatch function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check functions inputs ####
    if(length(assay.name) > 1 | is.logical(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if(nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
    }
    if (!ncg.selection.method %in% c('prod', 'sum', 'avergae', 'auto', 'non.overlap')){
        stop('The "ncg.selection.method" must be one of "prod", "sum", "avergae", "auto" or "non.overlap".')
    }
    if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" msut be a positve value  0 < top.rank.bio.genes < 1.')
    }
    if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    }
    if( grid.nb < 1 | grid.nb > nrow(se.obj)){
        stop(paste0('The "grid.nb" must be a positve value  0 < grid.nb < ', nrow(se.obj), '.'))
    }
    if(!grid.group %in% c('bio', 'uv', 'both')){
        stop('The "grid.group" must be one of "bio", "uv" or "non.overlap".')
    }
    if(!grid.direction %in% c('increase', 'decrease', 'auto')){
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
    if (isFALSE(assess.se.obj)) {
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        } else if (!is.null(variables.to.assess.ncg)) {
            if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg))) {
                stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.null(regress.out.uv.variables)){
        if (isTRUE(sum(regress.out.uv.variables %in% colnames(colData(se.obj))) != length(regress.out.uv.variables))) {
            stop('All or some of "regress.out.uv.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!is.null(regress.out.bio.variables)){
        if (isTRUE(sum(regress.out.bio.variables %in% colnames(colData(se.obj))) != length(regress.out.bio.variables))) {
            stop('All or some of "regress.out.bio.variables" cannot be found in the SummarizedExperiment object.')
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

    # Check the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, bio.groups, uv.groups, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose)
    }

    if (remove.na == 'none'){
        if (is.null(variables.to.assess.ncg))
            variables.to.assess.ncg <- c(bio.variables, uv.variables)
        mout <- lapply(
            variables.to.assess.ncg,
            function(x){
                if (sum(is.na(se.obj[[x]])) > 0)
                    stop('There are NA or missing values in the specified variables.')
            })
    }

    # Data transformation and normalization ####
    printColoredMessage(message = '-- Applying data transformation and normalization:',
        color = 'magenta',
        verbose = verbose
        )
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 + ',
                pseudo.count,
                ' (pseudo.count) on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(x = se.obj, i = assay.name)
    }

    ## normalization ####
    if (!is.null(normalization)) {
        printColoredMessage(
            message = '-- Applying data normalization:',
            color = 'magenta',
            verbose = verbose
            )
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
    # Regressing out variables ####
    if (!is.null(regress.out.uv.variables) | !is.null(regress.out.bio.variables)){
        printColoredMessage(
            message = '-- Regressing out unwanted or biological variables:',
            color = 'magenta',
            verbose = verbose
            )
    }
    ## regress out unwanted variables ####
    if(!is.null(regress.out.uv.variables)){
        printColoredMessage(
            message = '- Regressing out the specified unwanted variables:',
            color = 'blue',
            verbose = verbose
            )
        if(!is.null(normalization)){
            expr.data.reg.uv <- expr.data.nor
        } else expr.data.reg.uv <- expr.data
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' if they are largely associated with the ',
                paste0(bio.variables, collapse = ' & '),
                ' variables.'),
            color = 'red',
            verbose = verbose
            )
        expr.data.reg.uv <- t(expr.data.reg.uv)
        uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
        expr.data.reg.uv <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
    }

    ## regressing out biological variables ####
    if (!is.null(regress.out.bio.variables)){
        printColoredMessage(
            message = '- Regressing the specified biological variables:',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' if they are largely associated with the ',
                paste0(uv.variables, collapse = ' & '),
                ' variable(s).'),
            color = 'red',
            verbose = verbose
        )
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

    # Statistical analyses ####
    if (isFALSE(use.imf)){
        printColoredMessage(
            message = '-- Finding a subset of genes as negative control genes:',
            color = 'magenta',
            verbose = verbose
            )
        ## select genes that are highly affected by sources of unwanted variation ####
        printColoredMessage(
            message = '-- Selecting genes that are highly affected by each source(s) of unwanted variation:',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Note, this step will be performed within each homogeneous sample groups with respect to the biological variation.',
            color = 'blue',
            verbose = verbose
            )

        ### create all possible major homogeneous biological groups ####
        printColoredMessage(
            message = '- Creating all possible major homogeneous sample groups with respect to biological variables:',
            color = 'blue',
            verbose = verbose
            )
        if (is.null(bio.groups)){
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(bio.variables, collapse = ' & '),
                    ' variable(s) will be used to create all possible major homogeneous biological groups.'),
                color = 'blue',
                verbose = verbose
                )
            all.bio.groups <- createHomogeneousBioGroups(
                se.obj = se.obj,
                bio.variables = bio.variables,
                nb.clusters = nb.bio.clusters,
                clustering.method = bio.clustering.method,
                assess.se.obj = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )
        } else if (!is.null(bio.groups)) {
            printColoredMessage(
                message = paste0(
                    'The',
                    paste0(bio.groups, collapse = ' & '),
                    ' variables will be used to create all possible major homogeneous biological groups.'),
                color = 'blue',
                verbose = verbose)
            all.bio.groups <- createHomogeneousBioGroups(
                se.obj = se.obj,
                bio.variables = bio.groups,
                nb.clusters = nb.bio.clusters,
                clustering.method = bio.clustering.method,
                assess.se.obj = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)
        }

        ### correlation between gene expression and all continuous source of unwanted variation with each biological groups ####
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]]))
            )
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
        if (isTRUE(length(continuous.uv) > 0)) {
            printColoredMessage(
                message = paste0(
                    '-- Performing correlation analysis between gene expression and ',
                    'each specified continuous sources of unwanted variation:'),
                color = 'blue',
                verbose = verbose
            )
            selected.bio.groups <- findRepeatingPatterns(
                vec = all.bio.groups,
                n.repeat = min.sample.for.correlation
                )
            if (isTRUE(length(selected.bio.groups) > 0)) {
                if (is.null(regress.out.bio.variables)) {
                    data.to.use <- expr.data
                } else data.to.use <- expr.data.reg.bio
                corr.genes.uv <- lapply(
                    continuous.uv,
                    function(x) {
                        all.corr <- lapply(
                            selected.bio.groups,
                            function(y) {
                                selected.samples <- all.bio.groups == y
                                corr.genes <- as.data.frame(correls(
                                    x = t(data.to.use[, selected.samples]),
                                    y = se.obj@colData[[x]][selected.samples],
                                    type = corr.method,
                                    a = a ,
                                    rho = rho
                                    ))
                                corr.genes <- cbind(
                                    round(x = corr.genes[, 1:4], digits = 3),
                                    corr.genes[, 5, drop = FALSE]
                                    )
                                set.seed(2233)
                                corr.genes$ranked.genes <- rank(
                                    -abs(corr.genes[, 'correlation']),
                                    ties.method = 'random'
                                    )
                                row.names(corr.genes) <- row.names(data.to.use)
                                corr.genes
                            })
                        names(all.corr) <- selected.bio.groups
                        all.corr
                    })
                names(corr.genes.uv) <- continuous.uv
            } else if (isTRUE(length(selected.bio.groups) == 0))
                stop(paste0(
                    'There are not homogeneous biological groups that have at least ',
                    min.sample.for.correlation,
                    ' (min.sample.for.correlation) samples for correlation analysis.')
                    )
        } else corr.genes.uv <- NULL

        ### anova between gene expression and all categorical source of variation within each biological groups ####
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        if (isTRUE(length(categorical.uv) > 0)) {
            printColoredMessage(
                message = paste0(
                    '-- Performing ANOVA between individual gene expression and each ',
                    'specified categorical sources of unwanted variation:'),
                color = 'blue',
                verbose = verbose
                )
            anova.genes.uv <- lapply(
                categorical.uv,
                function(x) {
                    bio.batch <- table(all.bio.groups, colData(se.obj)[[x]])
                    cover.sample.groups <- rowSums(bio.batch >= min.sample.for.aov) == length(unique(se.obj[[x]]))
                    if (isTRUE(sum(cover.sample.groups) > 0)) {
                        printColoredMessage(
                            message = paste0(
                                sum(cover.sample.groups),
                                ' homogeneous biological group(s) have at least ',
                                min.sample.for.aov,
                                ' (min.sample.for.aov) samples within individual batches of the ',
                                x,
                                ' variable.'),
                            color = 'blue',
                            verbose = verbose)
                    }
                    if (isTRUE(sum(cover.sample.groups) == 0)){
                        printColoredMessage(
                            message = paste0(
                                'There are not homogeneous biological groups that have at least ',
                                min.sample.for.aov ,
                                ' (min.sample.for.aov) samples within each batches of the ',
                                x,
                                ' variable. This may result in unsatisfactory NCG selection.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    selected.bio.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                    if(isTRUE(length(selected.bio.groups) == 0)){
                        stop('There is not enough groups to perform ANOVA.')
                    }
                    if (is.null(regress.out.bio.variables)) {
                        data.to.use <- expr.data
                    } else data.to.use <- expr.data.reg.bio
                    all.anova <- lapply(
                        selected.bio.groups,
                        function(i) {
                            selected.samples <- all.bio.groups == i
                            if (anova.method == 'aov') {
                                anova.genes.batch <- as.data.frame(row_oneway_equalvar(
                                    x = data.to.use[, selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                            } else if (anova.method == 'welch.correction') {
                                anova.genes.batch <- as.data.frame(row_oneway_welch(
                                    x = data.to.use[, selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                            }
                            set.seed(2233)
                            anova.genes.batch$ranked.genes <- rank(-anova.genes.batch[, 'statistic'], ties.method = 'random')
                            anova.genes.batch
                        })
                    names(all.anova) <- selected.bio.groups
                    all.anova
                })
            names(anova.genes.uv) <- categorical.uv
        } else anova.genes.uv <- NULL

        ## select genes that are not highly affected by biological variation  ####
        printColoredMessage(
            message = '-- Selecting genes that are highly affected by each specified source(s) of biological variation:',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '- This step will be performed within each possible homogeneous unwanted groups.',
            color = 'blue',
            verbose = verbose
            )
        ### create all possible homogeneous uv groups ####
        printColoredMessage(
            message = '- Creating all possible homogeneous sample groups with respect to the specified unwanted variables:',
            color = 'blue',
            verbose = verbose
            )
        if (is.null(uv.groups)){
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(uv.variables, collapse = ' & '),
                    ' variables will be used as a major sources of unwanted variation',
                    ' to find all possible groups.'),
                color = 'blue',
                verbose = verbose)
            all.uv.groups <- createHomogeneousUVGroups(
                se.obj = se.obj,
                uv.variables = uv.variables,
                nb.clusters = nb.uv.clusters,
                clustering.method = uv.clustering.method,
                assess.se.obj = FALSE,
                save.se.obj = FALSE,
                verbose = verbose
                )
        } else if(!is.null(uv.groups)){
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(uv.groups, collapse = ' & '),
                    ' variables will be used as a major sources of unwanted variation',
                    ' to find all possible groups.'),
                color = 'blue',
                verbose = verbose
                )
            all.uv.groups <- createHomogeneousUVGroups(
                se.obj = se.obj,
                uv.variables = uv.groups,
                nb.clusters = nb.uv.clusters,
                clustering.method = uv.clustering.method,
                assess.se.obj = FALSE,
                save.se.obj = FALSE,
                verbose = verbose)
        }
        ### correlation between gene expression and all continuous source of biological variation with each uv groups ####
        bio.var.class <- unlist(lapply(
            bio.variables,
            function(x) class(colData(se.obj)[[x]])))
        continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
        if (length(continuous.bio) > 0) {
            printColoredMessage(
                message = '-- Correlation analyses:',
                color = 'magenta',
                verbose = verbose
                )
            selected.uv.groups <- findRepeatingPatterns(
                vec = all.uv.groups,
                n.repeat = min.sample.for.correlation)
            if (length(selected.uv.groups) > 0) {
                if (length(selected.uv.groups) == 1) {
                    group = 'group has'
                } else group = 'groups have'
                printColoredMessage(
                    message = paste0(
                        length(selected.uv.groups),
                        ' homogeneous groups with respect to the sources of unwanted variation ',
                        group,
                        ' at least ',
                        min.sample.for.correlation,
                        ' (min.sample.for.correlation) samples to pefrom correlation between gene-level',
                        'expression and all the continuous sources of bioloical variation.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (is.null(regress.out.uv.variables) & is.null(normalization)) {
                    data.to.use <- expr.data
                } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)) {
                    data.to.use <- expr.data.reg.uv
                } else if (is.null(regress.out.uv.variables) & !is.null(normalization)) {
                    data.to.use <- expr.data.nor
                }
                corr.genes.bio <- lapply(
                    continuous.bio,
                    function(x) {
                        all.corr <- lapply(
                            selected.uv.groups,
                            function(y) {
                                selected.samples <- all.uv.groups == y
                                corr.genes <- as.data.frame(correls(
                                    y = se.obj@colData[, x][selected.samples],
                                    x = t(data.to.use[, selected.samples]),
                                    type = corr.method,
                                    a = a ,
                                    rho = rho
                                    ))
                                corr.genes <- cbind(round(x = corr.genes[, 1:4], digits = 3),
                                                    corr.genes[, 5, drop = FALSE])
                                set.seed(2233)
                                corr.genes$ranked.genes <- rank(
                                    abs(corr.genes[, 'correlation']),
                                    ties.method = 'random'
                                    )
                                row.names(corr.genes) <- row.names(data.to.use)
                                corr.genes
                            })
                        names(all.corr) <- selected.uv.groups
                        all.corr
                    })
                names(corr.genes.bio) <- continuous.bio
            } else {
                stop(
                    paste0(
                        'There are not homogeneous groups with respect to sources of unwanted variation that have at least ',
                        min.sample.for.correlation,
                        ' (min.sample.for.correlation) samples for correlation analysis between gene-level expression and all',
                        ' the continuous sources of bioloical variation.'))
            }
        } else corr.genes.bio <- NULL

        ### anova between gene expression and all categorical source of biological variation with each uv groups ####
        categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
        if (length(categorical.bio) > 0) {
            printColoredMessage(
                message = '-- ANOV analyses:',
                color = 'magenta',
                verbose = verbose
                )
            anova.genes.bio <- lapply(
                categorical.bio,
                function(x) {
                    bio.batch <- table(all.uv.groups, colData(se.obj)[[x]])
                    cover.sample.groups <- rowSums(bio.batch >= min.sample.for.aov) == length(unique(se.obj[[x]]))
                    if (isTRUE(sum(cover.sample.groups) > 0)) {
                        printColoredMessage(
                            message = paste0(
                                sum(cover.sample.groups),
                                ' homogeneous unwanted group(s) have at least ',
                                min.sample.for.aov,
                                ' (min.sample.for.aov) samples within individual groups of the ', x, ' variable.'),
                            color = 'blue',
                            verbose = verbose)
                    }
                    if (isTRUE(sum(cover.sample.groups) == 0)){
                        printColoredMessage(
                            message = paste0(
                                'There are not homogeneous unwanted groups that have at least ',
                                min.sample.for.aov , ' (min.sample.for.aov) samples within each batches of the ',
                                x,
                                ' variable. This may result in unsatisfactory NCG selection.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    selected.uv.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                    if (length(selected.uv.groups) == 0) {
                        stop(paste0(
                            'It seems there is complete association between ',x,
                            ' homogeneous groups with respect to unwanted variation.'))
                    }
                    if (is.null(regress.out.uv.variables) &is.null(normalization)) {
                        data.to.use <- expr.data
                    } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)) {
                        data.to.use <- expr.data.reg.uv
                    } else if (is.null(regress.out.uv.variables) & !is.null(normalization)) {
                        data.to.use <- expr.data.nor
                    }
                    all.anova <- lapply(
                        selected.uv.groups,
                        function(i) {
                            selected.samples <- all.uv.groups == i
                            if (anova.method == 'aov') {
                                anova.gene.bio <- as.data.frame(row_oneway_equalvar(
                                    x = data.to.use[, selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                            } else if (anova.method == 'welch.correction') {
                                anova.gene.bio <- as.data.frame(row_oneway_welch(
                                    x = data.to.use[, selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                            }
                            set.seed(2233)
                            anova.gene.bio$ranked.genes <- rank(anova.gene.bio[, 'statistic'], ties.method = 'random')
                            anova.gene.bio
                        })
                    names(all.anova) <- selected.uv.groups
                    all.anova
                })
            names(anova.genes.bio) <- categorical.bio
        } else anova.genes.bio <- NULL
    }

    # Intermediate file ####
    ## read intermediate file ####
    if (isTRUE(use.imf)){
        if(is.null(imf.name)){
            imf.name <- paste0(assay.name, '|PerBiologyPerBatch|', ncg.selection.method)
        }
        if(is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
        anova.genes.bio <- all.tests$anova.genes.bio
        corr.genes.bio <- all.tests$corr.genes.bio
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
            imf.name <- paste0(assay.name, '|PerBiologyPerBatch|', ncg.selection.method)
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
    printColoredMessage(message = '-- Selection a set of genes as NCG:',
                        color = 'magenta',
                        verbose = verbose)

    ## product, sum or average of ranks ####
    if (ncg.selection.method %in% c('prod', 'sum', 'average')) {
        all.tests <-
            c('anova.genes.bio',
              'corr.genes.bio',
              'anova.genes.uv',
              'corr.genes.uv')
        all.stats <- lapply(
            all.tests,
            function(x) {
                if (!is.null(x)) {
                    temp.data <- get(x)
                    ranks.data <-
                        lapply(
                            names(temp.data),
                            function(y) {
                                all.ranks <- lapply(
                                    names(temp.data[[y]]),
                                    function(i) temp.data[[y]][[i]]$ranked.genes)
                                names(all.ranks) <- names(temp.data[[y]])
                                all.ranks <- do.call(cbind, all.ranks)
                                all.ranks
                            })
                    ranks.data <- do.call(cbind, ranks.data)
                    names(ranks.data) <- names(temp.data)
                    ranks.data
                }
            })
        all.stats <- as.data.frame(do.call(cbind, all.stats))
        row.names(all.stats) <- row.names(se.obj)
        ### product of ranks ####
        if (ncg.selection.method == 'prod') {
            printColoredMessage(
                message = 'A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose
                )
            all.stats$all.rank <- 10 ^ (rowSums(log(all.stats, base = 10)))
            if (sum(is.infinite(all.stats$all.rank)) > 0) {
                stop('The product of ranks results in infinity values.')
            }
        }
        ## sum of ranks ####
        if (ncg.selection.method == 'sum') {
            printColoredMessage(
                message = 'A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            all.stats$all.rank <- rowSums(x = all.stats, na.rm = TRUE)
        }
        ## average of ranks ####
        if (ncg.selection.method == 'average') {
            printColoredMessage(
                message = 'A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose
                )
            all.stats$all.rank <- rowMeans(x = all.stats, na.rm = TRUE)
        }

        all.stats <- all.stats[order(all.stats$all.rank, decreasing = FALSE), ]
        ncg.selected <- row.names(all.stats)[1:round(nb.ncg * nrow(se.obj), digits = 0)]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## non.overlap approach ####
    if (ncg.selection.method == 'non.overlap') {
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude all top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### select genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x) {
                if (isTRUE(!is.null(x))) {
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y) {
                            all.ranks <- sapply(
                                names(temp.data[[y]]),
                                function(z) temp.data[[y]][[z]]$ranked.genes)
                            set.seed(2233)
                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                        })))
                }
            })))
        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x) {
                temp <- get(x)
                if (length(names(temp)) != 0) {
                    ranks.data <- lapply(
                        names(temp),
                        function(y) {
                            unlist(lapply(names(temp[[y]]),
                                          function(z) {
                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                              row.names(temp[[y]][[z]])[index]
                                          }))
                        })
                }
            })))
        ## select of NCGS ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## auto approach ####
    if (isTRUE(ncg.selection.method == 'auto')){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### find highly affected genes by biology ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x) {
                if (isTRUE(!is.null(x))) {
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y) {
                            all.ranks <- sapply(
                                names(temp.data[[y]]),
                                function(z) temp.data[[y]][[z]]$ranked.genes)
                            set.seed(2233)
                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                        })))
                }
            })))
        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(c(top.rank.uv.genes * nrow(se.obj)), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x) {
                temp <- get(x)
                if (length(names(temp)) != 0) {
                    ranks.data <- lapply(
                        names(temp),
                        function(y) {
                            unlist(lapply(names(temp[[y]]),
                                          function(z) {
                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                              row.names(temp[[y]][[z]])[index]
                                          }))
                        })
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

        if (length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges)) {
            if(isTRUE(nb.ncg > length(ncg.selected))){
                con <- parse(text = paste0("nb.ncg", ">", "length(ncg.selected)"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is less than the number (',
                        nb.ncg ,
                        ') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose
                    )
            }
            if (isTRUE(nb.ncg < length(ncg.selected))){
                con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is larger than the number (',
                        nb.ncg ,
                        ') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ## grid search ####
            ### grid group: both bio and uv variable ####
            if (grid.group == 'both'){
                #### increasing order ####
                if (grid.direction == 'increase'){
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                        )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x) {
                                temp <- get(x)
                                if (length(names(temp)) != 0) {
                                    ranks.data <- lapply(
                                        names(temp),
                                        function(y) {
                                            unlist(lapply(names(temp[[y]]),
                                                          function(z) {
                                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                              row.names(temp[[y]][[z]])[index]
                                                          }))
                                        })
                                }
                            })))
                        # bio
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 0) top.rank.bio.genes.nb = 1
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x) {
                                if (isTRUE(!is.null(x))) {
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y) {
                                            all.ranks <- sapply(
                                                names(temp.data[[y]]),
                                                function(z) temp.data[[y]][[z]]$ranked.genes)
                                            set.seed(2233)
                                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                        })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if(length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                ### decreasing order ####
                if (grid.direction == 'decrease'){
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 0 ) top.rank.uv.genes.nb = 1
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x) {
                                temp <- get(x)
                                if (length(names(temp)) != 0) {
                                    ranks.data <- lapply(
                                        names(temp),
                                        function(y) {
                                            unlist(lapply(names(temp[[y]]),
                                                          function(z) {
                                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                              row.names(temp[[y]][[z]])[index]
                                                          }))
                                        })
                                }
                            })))
                        # bio
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if(top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x) {
                                if (isTRUE(!is.null(x))) {
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y) {
                                            all.ranks <- sapply(
                                                names(temp.data[[y]]),
                                                function(z) temp.data[[y]][[z]]$ranked.genes)
                                            set.seed(2233)
                                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                        })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if(length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
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
                        '- Updating the selection. Select top ',
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
                ###### increasing order ####
                if(grid.direction == 'increase'){
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x) {
                                if (isTRUE(!is.null(x))) {
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y) {
                                            all.ranks <- sapply(
                                                names(temp.data[[y]]),
                                                function(z) temp.data[[y]][[z]]$ranked.genes)
                                            set.seed(2233)
                                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                        })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        ncg.selected
                    }
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    lo <- nrow(se.obj) - top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if(top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                        top.bio.genes <- unique(unlist(lapply(
                            all.bio.tests,
                            function(x) {
                                if (isTRUE(!is.null(x))) {
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y) {
                                            all.ranks <- sapply(
                                                names(temp.data[[y]]),
                                                function(z) temp.data[[y]][[z]]$ranked.genes)
                                            set.seed(2233)
                                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                        })))
                                }
                            })))
                    }
                    ##### check selection ####
                    if(length(ncg.selected) == 0) stop('NCGs cannot be found based on the current parameters.')
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
            }
            ##### grid group: bio ####
            if (grid.group == 'uv'){
                printColoredMessage(
                    message = '- The grid search will be applied on unwanted factor. ',
                    color = 'blue',
                    verbose = verbose)
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- nrow(se.obj) - top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x) {
                                temp <- get(x)
                                if (length(names(temp)) != 0) {
                                    ranks.data <- lapply(
                                        names(temp),
                                        function(y) {
                                            unlist(lapply(names(temp[[y]]),
                                                          function(z) {
                                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                              row.names(temp[[y]][[z]])[index]
                                                          }))
                                        })
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        ncg.selected
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
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 0 ) top.rank.uv.genes.nb = 1
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x) {
                                temp <- get(x)
                                if (length(names(temp)) != 0) {
                                    ranks.data <- lapply(
                                        names(temp),
                                        function(y) {
                                            unlist(lapply(names(temp[[y]]),
                                                          function(z) {
                                                              index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                              row.names(temp[[y]][[z]])[index]
                                                          }))
                                        })
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        ncg.selected
                    }
                }
                ##### check selection ####
                if(length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 0)
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
        }
    }

    printColoredMessage(
        message = paste0('A set of ', sum(ncg.selected), ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose
    )

    ## Plotting ####
    printColoredMessage(
        message = '- Generating a heatmap plot of all the ranks of all the NCGs across all the variables',
        color = 'magenta',
        verbose = verbose
    )

    all.uv.bio.tests <- lapply(
        c(all.bio.tests, all.uv.tests),
        function(x){
            if (!is.null(x)){
                print(x)
                temp.data <- get(x)
                temp.data <- lapply(
                    names(temp.data),
                    function(y) {
                        tests.data <- sapply(
                            names(temp.data[[y]]),
                            function(z){
                                tm <- temp.data[[y]][[z]][ , 'ranked.genes']
                            })
                        tests.data <- rowMeans(tests.data)
                        tests.data <- data.frame(
                            ranks = tests.data,
                            group = rep(y, length(tests.data))
                            )
                        tests.data
                    })
                temp.data <- do.call(cbind , temp.data)
                temp.data
            }
        })
    all.uv.bio.tests <- do.call(cbind, all.uv.bio.tests)
    temp.data <- lapply(
        seq(1, ncol(all.uv.bio.tests), 2),
        function(x){
            temp.data <- all.uv.bio.tests[ , x, drop = FALSE]
            colnames(temp.data) <- all.uv.bio.tests[ , x+1][1]
            temp.data
        })
    temp.data <- do.call(cbind, temp.data)
    temp.data$ncg <- ncg.selected
    ha <- ComplexHeatmap::rowAnnotation(
        NCG = temp.data$ncg,
        col = list(NCG = c('TRUE' = 'gray10', 'FALSE' = 'gray'))
    )
    ncg.plot <- ComplexHeatmap::Heatmap(
        temp.data[ , seq_len(ncol(temp.data) - 1)],
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        right_annotation = ha,
        column_names_rot = 45,
        col = viridis::magma(n = 20),
        heatmap_legend_param = list(
            title = 'Ranks',
            title_gp = grid::gpar(fontsize = 14),
            by_row = TRUE,
            ncol = 1)
    )
    if (isTRUE(plot.output)) print(ncg.plot)

    # Assessment of selected set of NCG  ####
    ## pca ####
    if (isTRUE(assess.ncg)) {
        printColoredMessage(
            message = '-- Assessing the performance of selected NCG set:',
            color = 'magenta',
            verbose = verbose
            )
        if (is.null(variables.to.assess.ncg)) {
            all.variables <- c(bio.variables, uv.variables)
        } else all.variables <- variables.to.assess.ncg
        printColoredMessage(
            message = '- Performing PCA using only the selected genes as NCGs.',
            color = 'blue',
            verbose = verbose
            )
        if (isTRUE(apply.log)) {
            temp.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else temp.data <- assay(x = se.obj, i = assay.name)

        pca.data <- BiocSingular::runSVD(
            x = t(temp.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM =  bsparam(),
            center = center,
            scale = scale)$u

        ## regression and vector correlations ####
        printColoredMessage(
            message = paste0(
                '- Exploring the association of the first ',
                nb.pcs,
                '  PCs with the ',
                paste0(variables.to.assess.ncg, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose
            )
        all.corr <- lapply(
            all.variables,
            function(x) {
                if (class(se.obj[[x]]) %in% c('numeric', 'integer')) {
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared)
                } else if (class(se.obj[[x]]) %in% c('factor', 'character')) {
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <-catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y) {
                            cca <- cancor(x = pca.data[, 1:y, drop = FALSE], y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                }
            })
        pcs <- Groups <- NULL
        names(all.corr) <- all.variables
        pca.ncg <- as.data.frame(do.call(cbind, all.corr))
        pca.ncg['pcs'] <- c(1:nb.pcs)
        pca.ncg <- tidyr::pivot_longer(
            data = pca.ncg, -pcs,
            names_to = 'Groups',
            values_to = 'ls')
        p.assess.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), size = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab (expression("Correlations")) +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs))) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0, 1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(
                    size = 10,
                    angle = 25,
                    hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
            )
        if(isTRUE(plot.output )) print(p.assess.ncg)
    }
    # Save results ####
    ## add results to the SummarizedExperiment object ####
    if(is.null(output.name)){
        output.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|PbPbio:',
            ncg.selection.method,
            '|',
            assay.name)
    }
    if(is.null(ncg.group)){
        ncg.group <- paste0('ncg|supervised')
    }

    if (isTRUE(save.se.obj)) {
        printColoredMessage(message = '-- Saving a selected set of NCG to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if(!'supervised' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['supervised']] <- list()
        }
        if(!ncg.group %in% names(se.obj@metadata[['NCG']][['supervised']])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]] <- list()
        }
        if(!'ncg.set' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']] )){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']][[output.name]] <- ncg.selected

        if(isTRUE(assess.ncg)){
            if(!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['assessment.plot']] <- list()
            }
            if(!output.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['assessment.plot']] )){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['assessment.plot']][[output.name]] <- list()
            }
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['assessment.plot']][[output.name]] <- p.assess.ncg
        }
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The findNcgPerBiologyPerBatch function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## add results to the SummarizedExperiment object ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The findNcgPerBiologyPerBatch function finished.',
                            color = 'white',
                            verbose = verbose)
        return(ncg.selected)
    }
}


