#' Finds a set of negative control genes using two-way ANOVA.

#' @author Ramyar Molania

#' @description
#' This function utilizes two-way ANOVA to identify a set of genes suitable as negative control genes (NCG) for RUVIIIPRPS
#' normalization. Prior knowledge of biological and unwanted variation sources is necessary and should be specified. The
#' function begins by creating all possible sample groups based on biological and unwanted variation separately.
#' Subsequently, these groups are used as factors in two-way ANOVA to identify genes highly influenced by biological and
#' unwanted variation. Finally, the function selects genes with the possible highest F-statistics for unwanted variation and
#' lowest F-statistics for biological variation. Various approaches are employed for the final gene selection; please refer
#' to the details for more information.

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
#' * The product, sum or average of ranks of F-statistics is calculated. Then, the function selects 'nb.ncg' 'numbers of
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
#' Assess the performance of NCGS:
#' * The function can assess the initial performance of selected NCGs. This analysis involves principal component analysis
#' on only the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs' first principal components
#' and with the specified variables. Ideal NCGS, should show high and low R^2 or vector correlation for unwanted and
#' biological variation respectively.

#' @references
#' * Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' * Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A string indicating the name of the assay in the SummarizedExperiment object. This assay
#' should be the one that will be used as input data for the RUV-III normalization. We recommend using raw data,
#' if available.
#' @param bio.variables Character. Indicating the column names that contain known biological variable(s) in the
#' SummarizedExperiment object. These biological variables can be categorical or continuous. Continuous variables will be
#' divided into 'nb.bio.clusters' groups based on a clustering method selected in the 'bio.clustering.method' argument.
#' This argument cannot be empty.
#' @param uv.variables Character. Indicating the column names that contain unwanted variable(s) in the SummarizedExperiment
#' object. These unwanted variables can be categorical or continuous. Continuous variables will be divided into
#' 'nb.uv.clusters' groups based on a clustering method selected in the 'uv.clustering.method' argument. This argument
#' cannot be empty.
#' @param ncg.selection.method Character. Indicating how to summarize F-statistics obtained from two-way ANOVA and select a
#' set of genes as negative control genes. The options are 'prod', 'average', 'sum', 'non.overlap', 'auto' and 'quantile'.
#' The default is set to 'non.overlap'. For more information, refer to the details of the function.
#' @param nb.ncg Numeric. Indicating how many genes should be selected as NCG. The value represents the proportion of the
#' total genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param top.rank.bio.genes Numeric. Indicating the top-ranked genes that are highly affected by the biological variation.
#' This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'. The default
#' is set to 0.5.
#' @param top.rank.uv.genes Numeric. Indicating the top-ranked genes that are highly affected by the unwanted variation.
#' This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'. The default
#' is set to 0.5.
#' @param bio.percentile Numeric. The percentile cut-off of F-statistics to select genes that are highly affected by
#' the biological variation. The default is set to 0.8.
#' @param uv.percentile Numeric. The percentile cut-off of F-statistics to select genes that are highly affected by
#' the unwanted variation. The default is set to 0.8.
#' @param grid.group Character. Indicating whether the grid search should be performed on biological ('top.rank.bio.genes'),
#' unwanted ('top.rank.uv.genes') or both factors when the 'ncg.selection.method' is set to 'auto'. The options are
#' 'bio', 'uv', or 'both'. If set to 'both', the grid search will be performed on both biological and unwanted factors.
#' If set to 'bio' or 'uv', the grid search will be performed only on biological or unwanted factors. The default is
#' set to 'uv'.
#' @param grid.direction Character. Indicating whether the grid search should be performed in decreasing or increasing
#' order. The options are 'increase' or 'decrease'. The default is set to 'decrease'.
#' @param grid.nb Numeric. Indicating the number of genes for grid search when the 'ncg.selection.method' is set to 'auto'.
#' In the 'auto' approach, the grid search increases or decreases the initial values of 'top.rank.bio.genes' or
#' 'top.rank.uv.genes' or 'both' to find ~'nb.ncg' of genes as NCGs. The default is set to 20.
#' @param bio.clustering.method Character. Indicating which clustering methods should be used to group continuous sources
#' of biological variation. Refer to the 'createHomogeneousBioGroups' function for more details. The default is set to
#' 'kmeans' clustering.
#' @param nb.bio.clusters Numeric. Indicating the number of clusters for each continuous source of biological variation.
#' By default, it is set to 2. This means individual continuous sources of biological variation will be divided into two
#' groups.
#' @param uv.clustering.method Character. Indicates which clustering methods should be used to group continuous sources
#' of unwanted variation. Refer to the 'createHomogeneousUvGroups' function for more details. The default is set to
#' 'kmeans' clustering.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous source of unwanted variation.
#' By default, it is set to 2. This means individual continuous sources of unwanted variation will be divided into two
#' groups.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default is set to 'TRUE'.
#' @param pseudo.count Numeric. A value to be added as a pseudo count to all measurements before log transformation. The
#' default is 1.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not. This analysis involves
#' principal component analysis on only the selected NCG and then exploring the R² or vector correlation between the 'nb.pcs'
#' first principal components and the specified variables.
#' @param variables.to.assess.ncg Character. Indicates the column names of the SummarizedExperiment object that contain
#' variables whose association with the selected genes as NCG needs to be evaluated. The default is 'NULL'. This means all
#' the variables specified in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicating the number of the first principal components of selected NCG to be used to assess
#' the performance of NCGs. The default is 5.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis. Refer to the
#' 'computePCA' function for more details. The default is set to 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data before applying principal component analysis. Refer to the
#' 'computePCA' function for more details. The default is 'FALSE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeObj' will be applied. The default is set to 'TRUE'.
#' @param remove.na Character. Indicates whether to remove NA or missing values from either the 'assays', 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contain NA or missing values will be excluded. If
#' 'sample.annotation' is selected, the samples that contain NA or missing values for any 'bio.variables' and
#' 'uv.variables' will be excluded. The default is set to 'none'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result as a logical vector. The default is set to 'TRUE'. The file will be saved in
#' 'se.obj@metadata$NCG$supervised$output.name".
#' @param output.name Character. A string representing the output’s name. If set to 'NULL', the function will choose a name
#' automatically. In this case, the file name will be constructed as paste0(sum(ncg.selected), '|', paste0(bio.variables,
#' collapse = '&'), '|', paste0(uv.variables, collapse = '&'), '|TWAnova:', ncg.selection.method, '|', assay.name).
#' @param ncg.group Character. A string indicating the name of the group of NCG.
#' @param plot.output Logical. Indicates whether to plot the output or not.
#' @param save.imf Logical. Indicates whether to save the intermediate file in the SummarizedExperiment object. If set to 'TRUE',
#' the function saves the results of the two-way ANOVA. Subsequently, if users wish to adjust parameters such as 'nb.ncg',
#' 'ncg.selection.method', 'top.rank.bio.genes', and 'top.rank.uv.genes', the two-way ANOVA will not be recalculated. This
#' accelerates parameter tuning for NCG selection. The default value is 'FALSE'.
#' @param imf.name Character. Indicates the name to use when saving the intermediate file. If set to 'NULL', the function
#' will create a name. In this case, the file name will be constructed as
#' paste0(assay.name, '|TwoWayAnova|', ncg.selection.method). The default is 'NULL'.
#' @param use.imf Logical. Indicates whether to use the intermediate file
#' @param verbose Logical. If 'TRUE', shows messages of different steps of the function.

#' @importFrom dplyr mutate progress_estimated
#' @importFrom BiocSingular runSVD bsparam
#' @importFrom SummarizedExperiment assay
#' @importFrom fastDummies dummy_cols
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @importFrom stats aov
#' @import ggplot2
#' @export

findNcgByTwoWayAnova <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.7,
        top.rank.uv.genes = 0.7,
        bio.percentile = 0.5,
        uv.percentile = 0.5,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
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
    printColoredMessage(
        message = '------------The findNcgByTwoWayAnova function starts:',
        color = 'white',
        verbose = verbose
        )

    # Checking inputs ####
    if(!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name) | assay.name == 'all'){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty or "NULL".')
    }
    if (is.null(uv.variables)){
        stop('The "uv.variables" cannot be empty or "NULL".')
    }
    if(!is.vector(bio.variables) | !is.vector(uv.variables) ){
        stop('The "uv.variables" and "bio.variables" must be a vector of variables name(s) in the SummarizedExperiment object.')
    }
    if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('Individual specified variable must be either in the "bio.variables" or "uv.variables".')
    }
    if(!is.numeric(nb.ncg)){
        stop('The "nb.ncg" must be a positve numeric value 0 < nb.ncg < 1.')
    }
    if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" must be a positve value 0 < nb.ncg < 1.')
    }
    if (length(ncg.selection.method) > 1 | is.logical(ncg.selection.method)){
        stop('The "ncg.selection.method" muat be one of "prod", "sum", "average", "auto", "non.overlap" or "quantile".')
    }
    if (!ncg.selection.method %in% c('prod', 'sum', 'average', 'auto', 'non.overlap', 'quantile')){
        stop('The "ncg.selection.method" muat be one of "prod", "sum", "average", "auto", "non.overlap" or "quantile".')
    }
    if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE".')
    }
    if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE".')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    }
    if (isFALSE(is.logical(save.se.obj))) {
        stop('The "save.se.obj" must be "TRUE" or "FALSE".')
    }
    if (isFALSE(is.logical(use.imf))) {
        stop('The "use.imf" must be "TRUE" or "FALSE".')
    }
    if (isFALSE(is.logical(save.imf))) {
        stop('The "save.imf" must be "TRUE" or "FALSE".')
    }

    if (ncg.selection.method %in% c('non.overlap' , 'auto')){
        if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
            stop('The "top.rank.bio.genes" must be a positve value  0 < top.rank.bio.genes =< 1.')
        } else if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
            stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes =< 1.')
        }
    }

    if (isTRUE(apply.log)){
        if(length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
            stop('The "pseudo.count" must be 0 or a postive integer value.')
    }

    if (isTRUE(ncg.selection.method == 'quantile')){
        if(is.null(bio.percentile) | is.null(uv.percentile))
            stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
        if(bio.percentile > 1 | bio.percentile < 0)
            stop('The "bio.percentile" must be a postive value between 0 and 1.')
        if(uv.percentile > 1 | uv.percentile < 0)
            stop('The "uv.percentile" must be a postive value between 0 and 1.')
    }

    if (isFALSE(assess.se.obj)) {
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

    # Checking the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose
            )
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

    # Perform gene level two_way ANOVA ####
    if (isFALSE(use.imf)){
        ## Data transformation ####
        printColoredMessage(
            message = '-- Applying data transformation:',
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
                verbose = verbose
                )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0(
                    '- The ',
                    assay.name,
                    ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose
                )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        ## Creating all possible homogeneous sample groups ####
        ### biological groups ####
        printColoredMessage(
            message = '-- Creating all possible groups with respect to the specified sources of biological variation:',
            color = 'magenta',
            verbose = verbose
        )
        all.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = remove.na,
            verbose = verbose
            )

        ### unwanted groups ####
        printColoredMessage(
            message = '-- Creating all possible groups with respect to the specified sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose
            )
        all.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            nb.clusters = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )

        ## Two_way ANOVA ####
        printColoredMessage(
            message = '-- Two_way ANOVA:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Perform Two_way ANOVA:',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- This is between all individual gene-level expression',
                'and considering both biological and unwanted variables created above as factors.'),
            color = 'blue',
            verbose = verbose
            )
        ### two_way ANOVA ####
        all.aov <- aov(t(expr.data) ~ all.bio.groups + all.uv.groups)
        all.aov <- summary(all.aov)
        all.aov <- as.data.frame(t(sapply(
            c(1:nrow(se.obj)),
            function(x) all.aov[[x]]$`F value`[1:2]
            )))
        all.aov <- round(x = all.aov, digits = 1)
        colnames(all.aov) <- c('Biology', 'UV')
        row.names(all.aov) <- row.names(se.obj)
        ### rank the F-statistics ####
        bio.rank <- uv.rank <- Biology <- UV <- NULL
        set.seed(2190)
        all.aov$bio.rank <- rank(x = all.aov$Biology, ties.method = 'random')
        set.seed(2190)
        all.aov$uv.rank <- rank(x = -all.aov$UV, ties.method = 'random')
        all.aov$uv.rank.plot <- rank(x = all.aov$UV, ties.method = 'random')
    }
    # Plotting the F-statistics ####
    if (isTRUE(plot.output)){
        p.fvalues <- ggplot(all.aov, aes(x = log2(Biology), y = log2(UV)) ) +
            geom_point(alpha = .1) +
            xlab('Biology (log2 of F-values)') +
            ylab('Unwanted variation (log2 of F-values)') +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                plot.title = element_text(size = 16))
        p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot)) +
            geom_point(alpha = .1) +
            xlab('Biology (rank of F-values)') +
            ylab('Unwanted variation (rank of F-values)') +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                plot.title = element_text(size = 16))
        all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
        print(all.plots)
        rm(p.fvalues, p.fvalues.rank)
    }

    # Read the intermediate file ####
    if (isTRUE(use.imf)){
        printColoredMessage(
            message = '- Retrieving the results of two-way ANOVA from the the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (is.null(imf.name)){
            imf.name <- paste0(assay.name, '|TwoWayAnova|', ncg.selection.method)
        }
        if (is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.aov <- se.obj@metadata$IMF$NCG[[imf.name]]
    }
    # Save the intermediate file ####
    if (isTRUE(save.imf)){
        printColoredMessage(
            message = '-- Saving a intermediate file:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- The results of two-way ANOVA is saved in the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (length(se.obj@metadata$IMF) == 0 ) {
            se.obj@metadata[['IMF']] <- list()
        }
        if (!'NCG' %in% names(se.obj@metadata[['IMF']])){
            se.obj@metadata[['IMF']][['NCG']] <- list()
        }
        if (is.null(imf.name)){
            imf.name <- paste0(assay.name, '|TwoWayAnova|', ncg.selection.method)
        }
        if (!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
        }
        se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- all.aov
    }

    # Selection of NCG ####
    printColoredMessage(
        message = '-- Selecting a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose
        )
    ## Product, average and sum of ranks ####
    if (ncg.selection.method %in% c('prod', 'average', 'sum')){
        ### product of ranks ####
        if(isTRUE(ncg.selection.method == 'prod')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- all.aov$bio.rank * all.aov$uv.rank
            if(sum(is.infinite(all.aov$all.rank)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## average of ranks ####
        if (isTRUE(ncg.selection.method == 'average')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- rowMeans(all.aov[ , c('bio.rank', 'uv.rank')])
        }
        ## sum of ranks ####
        if(isTRUE(ncg.selection.method == 'sum')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- rowSums(all.aov[ , c('bio.rank', 'uv.rank')])
        }
        ### select top genes as NCGS ####
        nb.ncg <- round(x = nb.ncg * nrow(se.obj), digits = 0)
        printColoredMessage(
            message = paste0('- Select ', nb.ncg , ' genes as NCGS.'),
            color = 'blue',
            verbose = verbose)
        all.aov <- all.aov[order(all.aov$all.rank, decreasing = FALSE), ]
        ncg.selected <- row.names(all.aov)[1:nb.ncg]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        ncg <- NULL
    }
    ## Quantile approach ####
    if (ncg.selection.method == 'quantile'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "quantile" approach..',
            color = 'blue',
            verbose = verbose
            )
        ### find biological percentile ####
        bio.quan <- quantile(x = all.aov$Biology , probs = bio.percentile)
        top.bio.genes <- row.names(all.aov)[all.aov$Biology > bio.quan]

        ## find UV percentile ####
        uv.quan <- quantile(x = all.aov$UV , probs = uv.percentile)
        top.uv.genes <- row.names(all.aov)[all.aov$UV > uv.quan]
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                length(top.uv.genes),
                ' genes with the UV F-statistics higher than ',
                uv.quan,
                ' (' ,
                uv.percentile* 100,
                '% percentile), and exclude any genes presents in ',
                length(top.bio.genes),
                ' genes with the biological F-statistics higher than ',
                bio.quan,
                ' (' ,
                bio.percentile* 100,
                '% percentile).'),
            color = 'blue',
            verbose = verbose)

        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
        ## plot ####
        if(isTRUE(plot.output)){
            all.aov$ncg <- row.names(all.aov) %in%  row.names(se.obj)[ncg.selected]
            p.fvalues <- ggplot(all.aov, aes(x = log2(Biology), y = log2(UV), color = ncg)) +
                geom_point(alpha = .1) +
                xlab('Biology (log2 of F-values)') +
                ylab('Unwanted variation (log2 of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot, color = ncg)) +
                geom_point(alpha = .1) +
                xlab('Biology (rank of F-values)') +
                ylab('Unwanted variation (rank of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
        }
    }
    # Non.overlap approach ####
    if (ncg.selection.method == 'non.overlap'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### select genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1-top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
        ## select genes affected by unwanted variation ####
        top.rank.uv.genes <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes]
        ## select of NCGS ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        if(isTRUE(sum(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        ## plot ####
        if(isTRUE(plot.output)){
            all.aov$ncg <- row.names(all.aov) %in%  row.names(se.obj)[ncg.selected]
            p.fvalues <- ggplot(all.aov, aes(x = log2(Biology), y = log2(UV), color = ncg)) +
                geom_point(alpha = .1) +
                xlab('Biology (log2 of F-values)') +
                ylab('Unwanted variation (log2 of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot, color = ncg)) +
                geom_point(alpha = .1) +
                xlab('Biology (rank of F-values)') +
                ylab('Unwanted variation (rank of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
        }
    }
    ## Auto approach ####
    if (ncg.selection.method == 'auto'){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Select top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude all genes in top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### select genes affected by biological variation ####
        nb.ncg <- round(nb.ncg * nrow(se.obj), digits = 0)
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank < top.rank.uv.genes.nb]
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        printColoredMessage(
            message = paste0(
                '- ',
                length(ncg.selected),
                ' genes are found.'),
            color = 'blue',
            verbose = verbose
            )

        ncg.ranges <- round(x = 0.01 *nb.ncg, digits = 0)
        if (length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges) ){
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
                    verbose = verbose)

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
            #### grid search ####
            ##### grid group: both bio and uv variable ####
            if (grid.group == 'both'){
                printColoredMessage(
                    message = '- The grid search will be applied on both biological and unwanted factor. ',
                    color = 'blue',
                    verbose = verbose
                    )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                        )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if(top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                # genes selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                ## bio
                top.rank.bio.genes <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes/nrow(se.obj) * 100, digits = 2)
                if(top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                ## uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if(top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then excluding any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ##### grid group: bio ####
            if (grid.group == 'bio'){
                printColoredMessage(
                    message = '- The grid search will be applied on biological factor. ',
                    color = 'blue',
                    verbose = verbose
                    )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if(top.rank.bio.genes.nb < 1 ) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
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
                        if(top.rank.bio.genes.nb > nrow(se.obj) ) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.aov)[ all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 2)
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
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- nrow(se.obj) - top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if(top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if(top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if(top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes * 100,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ## plotting ####
            if(isTRUE(plot.output)){
                all.aov$ncg <- row.names(all.aov) %in%  row.names(se.obj)[ncg.selected]
                p.fvalues <- ggplot(all.aov, aes(x = log2(Biology), y = log2(UV), color = ncg)) +
                    geom_point(alpha = .1) +
                    xlab('Biology (log2 of F-values)') +
                    ylab('Unwanted variation (log2 of F-values)') +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        plot.title = element_text(size = 16))
                p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot, color = ncg)) +
                    geom_point(alpha = .1) +
                    xlab('Biology (rank of F-values)') +
                    ylab('Unwanted variation (rank of F-values)') +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        plot.title = element_text(size = 16))
                all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
                print(all.plots)
                }
            } else {
            printColoredMessage(
                message = paste0(length(ncg.selected), ' genes are selected as NCGs.'),
                color = 'blue',
                verbose = verbose)
        }
    }
    printColoredMessage(
        message = paste0(
            '- ',
            sum(ncg.selected),
            ' genes are selected as negative control genes.'),
        color = 'blue',
        verbose = verbose
        )

    # Performance assessment of selected NCG  ####
    if (isTRUE(assess.ncg)){
        printColoredMessage(
            message = '-- Assessing the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Performing PCA on only selected genes as NCG.',
            color = 'blue',
            verbose = verbose
            )
        if (is.null(variables.to.assess.ncg))
            variables.to.assess.ncg <- c(bio.variables, uv.variables)
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
        ### pca ####
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u

        # regression and vector correlations ####
        all.corr <- lapply(
            variables.to.assess.ncg,
            function(x){
                if (class(se.obj[[x]]) %in% c('numeric', 'integer')){
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared
                        )
                } else if (class(se.obj[[x]]) %in% c('factor', 'character')){
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y){ cca <- cancor(
                            x = pca.data[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)
                        })
                }
            })
        names(all.corr) <- variables.to.assess.ncg
        pca.ncg <- as.data.frame(do.call(cbind, all.corr))
        pcs <- Groups <- NULL
        pca.ncg['pcs'] <- c(1:nb.pcs)
        pca.ncg <- tidyr::pivot_longer(
            data = pca.ncg,
            -pcs,
            names_to = 'Groups',
            values_to = 'ls')
        # plot ####
        p.assess.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), linewidth = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab (expression("Correlations")) +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                plot.title = element_text(size = 16)
            )
        if (isTRUE(plot.output)) print(p.assess.ncg)
    }

    # Save results ####
    ### add results to the SummarizedExperiment object ####
    if (is.null(ncg.group)){
        ncg.group <- paste0('ncg|supervised')
    }
    if (is.null(output.name)){
        output.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|TWAnova:',
            ncg.selection.method,
            '|',
            assay.name)}

    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose
            )
        ## Check if metadata NCG already exists
        if (length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if (!'supervised' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['supervised']] <- list()
        }
        if (!ncg.group %in% names(se.obj@metadata[['NCG']][['supervised']])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]] <- list()
        }
        if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']] )){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['supervised']][[ncg.group]][['ncg.set']][[output.name]] <- ncg.selected

        if (isTRUE(assess.ncg)){
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group]])){
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
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ### export results as logical vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The set of NCGs is outpputed as a logical vector.',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose
            )
        return(ncg.selected)
    }
}




