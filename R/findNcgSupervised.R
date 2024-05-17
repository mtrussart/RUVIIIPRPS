#' Find a set of negative control genes using supervised approaches.

#' @author Ramyar Molania

#' @description
#' This function contains three different functions including 'findNcgByTwoWayAnova', 'findNcgAcrossSamples' and
#' 'findNcgPerBiologyPerBatch' to find a set of genes as negative control genes (NCG) for RUV-III-PRPS normalization. We
#' refer to each function for more details.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates the name of the assay in the SummarizedExperiment object. This assay should
#' be the one that will be used for RUV-III-PRPS normalization. We recommend to use raw data.
#' @param bio.variables Symbols. Indicating the column names that contain known biological variable(s) in the
#' SummarizedExperiment object. These biological variables can be categorical or continuous. Continuous variables will be
#' divided into 'nb.bio.clusters' groups based on a clustering method selected in the 'bio.clustering.method' argument.
#' This argument cannot be empty.
#' @param uv.variables Symbols. Indicating the column names that contain unwanted variable(s) in the SummarizedExperiment
#' object. These unwanted variables can be categorical or continuous. Continuous variables will be divided into
#' 'nb.uv.clusters' groups based on a clustering method selected in the 'uv.clustering.method' argument. This argument
#' cannot be empty.
#' @param approach Symbols. Indicats which NCGs selection method should be used. The options are 'AnovaCorr.PerBatchPerBio',
#' 'AnovaCorr.AcrossAllSamples' and'TwoWayAnova'. The default is set to 'TwoWayAnova'. Refer to details for more information.
#' @param ncg.selection.method Symbol. Indicating how to summarize F-statistics obtained from two-way ANOVA and select a
#' set genes as negative control genes. The options are 'prod', 'average', 'sum', 'non.overlap', 'auto' and 'quantile'.
#' The default is set to non.overlap'. For more information, refer to the details of the function.
#' @param nb.ncg Numeric. Indicating how many genes should be selected as NCG. The value represents the proportion of the
#' total genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param top.rank.bio.genes Numeric. Indicating the top ranked genes that are highly affected by the biological variation.
#' This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'. The default
#' is set to 0.5.
#' @param top.rank.uv.genes Numeric. Indicating the top ranked genes that are highly affected by the unwanted variation.
#' This is required to be specified when the 'ncg.selection.method' is set to either 'non.overlap' or 'auto'.The default
#' is set to 0.5.
#' @param bio.percentile Numeric. The percentile cut-off of F-statistics to select genes that are highly affected by
#' the biological variation. The default is set to 0.8.
#' @param uv.percentile Numeric. The percentile cut-off of F-statistics to select genes that are highly affected with
#' the unwanted variation. The default is set to 0.8.
#' @param grid.group Symbol. Indicating whether the grid search should be performed on biological ('top.rank.bio.genes'),
#' unwanted ('top.rank.uv.genes') or both factors, when the ncg.selection.method' is set to 'auto'.. The options are
#' 'bio', 'uv' or 'both'. If is set to 'both', the grid search will be performed on both biological and unwanted factors.
#' If is set to 'bio' or 'uv', the grid search will be performed only on biological or unwanted factors. The default is
#' set to 'uv'.
#' @param grid.direction Symbol. Indicating the grid search should be performed in decreasing or increasing order. The
#' option are 'increase' or 'decrease'. The default is set to 'decrease'.
#' @param grid.nb Numeric. Indicating the number of genes for grid search when the 'ncg.selection.method' is set to 'auto'.
#' In the 'auto' approach, the grid search increases or decreases the initial values of 'top.rank.bio.genes' or
#' 'top.rank.uv.genes' or 'both' to find ~'nb.ncg' of genes as NCGs. The default is set to 20.
#' @param bio.groups Symbol. A symbol or a vector of symbols indicating the columns names that contains biological variables
#' in the SummarizedExperiment object. If is specified, the 'bio.groups' will be used for grouping samples into different
#' homogeneous biological groups. If is 'NULL', the 'bio.variables' will be used for grouping samples into different
#' homogeneous biological groups.
#' @param bio.clustering.method Symbols. Indicating which clustering methods should be used to group continuous sources
#' of biological variation. Refer to the 'createHomogeneousBioGroups' function for more details. The default is set to
#' 'kmeans' clustering .
#' @param nb.bio.clusters Numeric. Indicating the number of clusters for each continuous sources of biological variation.
#' The by default it is set to 2. This means individual continuous sources of biological variation will be divided into two
#' groups.
#' @param uv.groups Symbol. A symbol or a vector of symbols indicating the columns names that contains biological variables
#' in the SummarizedExperiment object. If is specified, the 'uv.groups' will be used for grouping samples into possible
#' homogeneous sample groups with respect to unwanted variables. If is 'NULL', the 'uv.variables' will be used for grouping
#' samples.
#' @param uv.clustering.method Symbols. Indicates which clustering methods should be used to group continuous sources
#' of unwanted variation. Refer to the 'createHomogeneousUvGroups' function for more details. The default is set to
#' 'kmeans' cluster.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous sources of unwanted variation.
#' The by default it is set to 2. This means individual continuous sources of biological variation will be divided into two
#' groups.
#' @param normalization Symbol. Indicates which normalization method should be applied to the data before finding genes
#' that are affected by biological variation. The default is set to 'CPM'. If is 'NULL', no normalization will be applied.
#' Refer to the 'applyOtherNormalizations' function for more details.
#' @param regress.out.bio.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param regress.out.uv.variables Symbols. Indicates the columns names that contain unwanted variation variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by biological variation. The default is NULL, indicates the regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default it is set to 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param anova.method Indicates which anova method should be used to compute association between gene-level
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to be present in each group before applying
#' the ANOVA. The default is 3.
#' @param corr.method Symbol. Indicates which correlation method should be used to compute association between gene-level
#' expression and a continuous variable. The default is 'spearman'.
#' @param a The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to be considered in each group before
#' applying the correlation analysis. The default is 10.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not. This analysis involves
#' principal component analysis on only the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs'
#' first principal components and with the specified variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that contain
#' variables whose association with the selected genes as NCG needs to be evaluated. The default is 'NULL'. This means all
#' the variables specified in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicateing the number of the first principal components on selected NCG to be used to assess
#' the performance of NCGs. The default is 5.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis or not.
#' Refer to the 'computePCA' function for more details. The default is set to 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data before applying principal component analysis.
#' Refer to the 'computePCA' function for more details. The default is 'FALSE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeObj' will be applied. The default is set to 'TRUE'.
#' @param assess.variables Logical. Indicates whether to assess the correlation between biological and unwanted variation
#' variables separately when creating homogeneous sample groups. Refer to the function 'assessVariablesAssociation' for
#' more details. The default is set to 'FLASE'.
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
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If
#' 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and
#' 'uv.variables' will be excluded. The default is set to 'none'.
#' @param save.se.obj Logical. Indicating whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result as a logical vector. The default it is set to 'TRUE'. The file will be saved in
#' 'se.obj@metadata$NCG$supervised$output.name".
#' @param output.name Symbol. A representation for the output s name. If set to 'NULL', the function will choose a name
#' automatically. In this case, the file name will be constructed as paste0(sum(ncg.selected),'|', paste0(bio.variables,
#' collapse = '&'), '|', paste0(uv.variables, collapse = '&'),'|TWAnova:', ncg.selection.method, '|', assay.name).
#' @param use.imf Logical. Indicating whether to use the intermediate file or not. The default is set to 'FALSE'.
#' @param save.imf Logical. Indicating whether to save the intermediate file in the SummarizedExperiment object or not.
#' If set to 'TRUE', the function saves the results of the two-way ANOVA. Subsequently, if users wish to adjust parameters
#' such as 'nb.ncg', 'ncg.selection.method', 'top.rank.bio.genes', and 'top.rank.uv.genes', the two-way ANOVA will not
#' be recalculated. This accelerates parameter tuning for NCG selection. The default value is 'FALSE'.
#' @param imf.name Symbol. Indicating the name to use when saving the intermediate file. If set to 'NULL', the function
#' will create a name. In this case, the file name will be constructed as
#' paste0(assay.name, '|TwoWayAnova|', ncg.selection.method).The default is 'NULL'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return Either the SummarizedExperiment object containing the a set of negative control genes in the metadata  or a
#' logical vector of the selected negative control genes.

#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @export

findNcgSupervised <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        approach = 'TwoWayAnova',
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        normalization = 'CPM',
        regress.out.uv.variables = NULL,
        regress.out.bio.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        anova.method = 'aov',
        min.sample.for.aov = 3,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        min.sample.for.correlation = 10,
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.95, 0.95),
        cont.cor.coef = c(0.95, 0.95),
        remove.na = 'none',
        save.se.obj = TRUE,
        output.name = NULL,
        use.imf = FALSE,
        save.imf = FALSE,
        imf.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgSupervised function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(!approach %in% c('AnovaCorr.PerBatchPerBio', 'AnovaCorr.AcrossAllSamples', 'TwoWayAnova')){
        stop('The approach must be one of the "AnovaCorr.PerBatchPerBiology", "AnovaCorr.AcrossAllSamples" or "TwoWayAnova".')
    }

    # find NCGs ####
    ## find NCGs using AnovaCorr.PerBatchPerBio approach ####
    if(approach == 'AnovaCorr.PerBatchPerBiology'){
        se.obj <- findNcgPerBiologyPerBatch(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.bio.variables = regress.out.bio.variables,
            bio.groups = bio.groups,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            regress.out.uv.variables = regress.out.uv.variables,
            uv.groups = uv.groups,
            uv.clustering.method = uv.clustering.method,
            nb.uv.clusters = nb.uv.clusters,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            corr.method = corr.method,
            a = a,
            rho = rho,
            anova.method = anova.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    ## find NCGs using AnovaCorr.AcrossAllSamples approach ####
    if(approach == 'AnovaCorr.AcrossAllSamples'){
        se.obj <- findNcgAcrossSamples(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.bio.variables = regress.out.bio.variables,
            regress.out.uv.variables = regress.out.uv.variables,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            corr.method = corr.method,
            a = a,
            rho = rho,
            anova.method = anova.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    ## find NCGs using TwoWayAnova approach ####
    if(approach == 'TwoWayAnova'){
        se.obj <- findNcgByTwoWayAnova(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            uv.clustering.method = uv.clustering.method,
            nb.uv.clusters = nb.uv.clusters,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    printColoredMessage(message = '------------The findNcgSupervised function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
