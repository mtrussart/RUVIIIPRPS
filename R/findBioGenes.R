#' Find biological genes.

#' @author Ramyar Molania

#' @description
#' This function uses different approaches to find genes that highly affected by biological variables or highly
#' variable genes. Refer to the details for

#' @details
#' Additional details...
#' - TwoWayAnov: In this process, all the biological and unwanted variables will be grouped into two categorical variables
#' separately using the createHomogeneousBioGroups() and createHomogeneousUVGroups() functions. Then, a two-way ANOVA is
#' applied to the expression levels of individual genes, considering the summarized biological and unwanted variables as
#' two factors. A set of genes with the highest F-statistic the biological variables will be selected as biological genes.
#'
#'- AnovaCorr.AcrossAllSamples: computes gene-level correlation and ANOVA to find genes that are highly affected by continuous
#' and categorical sources biological variables across all samples, respectively.
#' Then, a set of genes with highest absolute correlation coefficients and F-statistics will be selected as biological.
#'
#' - AnovaCorr.PerBatchPerBio: computes gene-level correlation and ANOVA within groups of samples that are homogeneous
#' with respect to biological or unwanted variation. This is useful in situations where the biological and unwanted variation
#' are highly correlated.
#'
#' - mad.unsupervised: If sources of unwanted variation are unknown, the identifyUnknownUV must be applied to
#' estimate them. Subsequently, the functions conducts a Median Absolute Deviation (MAD) analysis on each gene within
#' sample groups homogeneous with respect to unwanted variables, to pinpoint genes highly
#' variable due to biological factors. The higher the MAD, the more likely the genes are to be biologically variable.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol that indicates the name of the assay in the SummarizedExperiment object.
#' @param approach Symbol. A symbol that indicates which biological gene selection methods should be used. The options
#' are 'AnovaCorr.PerBatchPerBio', 'AnovaCorr.AcrossAllSamples', 'TwoWayAnova' and 'mad.unsupervised'. The default is set
#' to 'TwoWayAnova'. Refer to the function details for more information.
#' @param nb.bio.genes Numeric. A numeric value indicates how many genes should be selected as biological genes. The value
#' represents the proportion of the total genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param bio.variables Symbol. A symbol or a vector of symbols indicating the column names that contain known biological
#' variable(s) in the SummarizedExperiment object. These biological variables can be categorical or continuous. Continuous
#' variables will be divided into 'nb.bio.clusters' groups based on a clustering method selected in the 'bio.clustering.method'
#' argument.
#' @param uv.variables Symbols. A symbol or a vector of symbols indicating the column names that contain unwanted variable(s)
#' in the SummarizedExperiment object. These unwanted variables can be categorical or continuous. Continuous variables
#' will be divided into nb.uv.clusters' groups based on a clustering method selected in the 'uv.clustering.method' argument.
#' @param bio.clustering.method Symbol. A symbol that indicates which clustering methods should be used to group continuous
#' sources of biological variation.  The default is set to kmeans' clustering. Refer to the 'createHomogeneousBioGroups'
#' function for more details.
#' @param nb.bio.clusters Numeric. A numeric value that indicates the number of clusters for each continuous sources of
#' biological variation. The default is set to 3. This means individual continuous sources of biological variation will
#' be divided into two groups based on the 'bio.clustering.method' method.
#' @param uv.groups Symbol. A symbol or a vector of symbols indicating the column names that contain known unwanted
#' variable(s) in the SummarizedExperiment object. These biological variables can be categorical or continuous. Continuous
#' variables will be divided into 'nb.bio.clusters' groups based on a clustering method selected in the 'bio.clustering.method'
#' argument.
#' @param uv.clustering.method Symbol. A symbol that indicates which clustering methods should be used to group continuous
#' sources of unwanted variation. The default is set to kmeans' cluster. Refer to the 'createHomogeneousUvGroups' function
#' for more details.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous sources of unwanted variation.
#' The by default it is set to 2. This means individual continuous sources of biological variation will be divided into two
#' groups.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default it is set to 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param normalization Symbol. Indicates which normalization method should be applied to the data before finding genes
#' that are affected by biological variation. The default is set to 'CPM'. If is 'NULL', no normalization will be applied.
#' Refer to the 'applyOtherNormalizations' function for more details.
#' @param regress.out.uv.variables Symbols. Indicates the columns names that contain unwanted variation variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by biological variation. The default is NULL, indicates the regression will not be applied.
#' @param assess.bio.genes Logical. Indicates whether to assess the performance of selected NCG or not. This analysis involves
#' principal component analysis on only the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs'
#' first principal components and with the specified variables.
#' @param variables.to.assess.bio.genes Symbols. Indicates the column names of the SummarizedExperiment object that contain
#' variables whose association with the selected genes as NCG needs to be evaluated. The default is 'NULL'. This means all
#' the variables specified in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicateing the number of the first principal components on selected NCG to be used to assess
#' the performance of NCGs. The default is 5.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis or not.
#' Refer to the 'computePCA' function for more details. The default is set to 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data before applying principal component analysis.
#' Refer to the 'computePCA' function for more details. The default is 'FALSE'.
#' @param anova.method Indicates which anova method should be used to compute association between gene-level
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to be present in each group before applying
#' the ANOVA. The default is 3.
#' @param corr.method Symbol. Indicates which correlation method should be used to compute association between gene-level
#' expression and a continuous variable. The default is 'spearman'.
#' @param a The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to be considered in each group before
#' applying the correlation analysis. The default is 10.
#' @param min.sample.for.mad Numeric. Indicates the minimum number of samples to perform Median Absolute Deviation on
#' each gene within homogeneous sample groups, considering unwanted variables. The default is 3.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeObj' will be applied. The default is set to 'TRUE'.
#' @param assess.variables Logical. Indicates whether to assess the correlation between biological and unwanted variation
#' variables separately when creating homogeneous sample groups. Refer to the function 'assessVariablesAssociation' for
#' more details. The default is set to 'FLASE'
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
#' collapse = '&'), '|', paste0(uv.variables, collapse = '&'),'|TWAnova:', '|', assay.name).
#' @param save.imf Logical. Indicating whether to save the intermediate file in the SummarizedExperiment object or not.
#' If set to 'TRUE', the function saves the results of the two-way ANOVA. Subsequently, if users wish to adjust parameters
#' such as 'nb.bio.genes', 'ncg.selection.method', 'top.rank.bio.genes', and 'top.rank.uv.genes', the two-way ANOVA will not
#' be recalculated. This accelerates parameter tuning for NCG selection. The default value is 'FALSE'.
#' @param imf.name Symbol. Indicating the name to use when saving the intermediate file. If set to 'NULL', the function
#' will create a name. In this case, the file name will be constructed as
#' paste0(assay.name, '|TwoWayAnova|'). The default is 'NULL'.
#' @param use.imf Logical. Indicating whether to use the intermediate file or not. The default is set to 'FALSE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return Either the SummarizedExperiment object containing the a set of biological genes in the metadata  or a
#' logical vector of the selected biological genes

#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @export


findBioGenes <- function(
        se.obj,
        assay.name,
        approach = 'TwoWayAnova',
        nb.bio.genes = 0.1,
        bio.variables,
        uv.variables,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        normalization = 'CPM',
        regress.out.uv.variables = NULL,
        assess.bio.genes = TRUE,
        variables.to.assess.bio.genes = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        anova.method = 'aov',
        min.sample.for.aov = 3,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        min.sample.for.correlation = 10,
        min.sample.for.mad = 3,
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        remove.na = 'none',
        save.se.obj = TRUE,
        output.name = NULL,
        save.imf = FALSE,
        imf.name = NULL,
        use.imf = FALSE,
        verbose = TRUE
        ){
    # Check inputs ####
    if(length(assay.name) > 1 | is.logical(assay.name) | assay.name == 'all'){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    } else if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty or "NULL".')
    } else if (is.null(uv.variables)){
        stop('The "uv.variables" cannot be empty or "NULL".')
    } else if(!is.vector(bio.variables) | !is.vector(uv.variables) ){
        stop('The "uv.variables" and "bio.variables" must be a vector of variables name(s) in the SummarizedExperiment object.')
    } else if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('Individual specified variable must be either in the "bio.variables" or "uv.variables".')
    } else if(!is.numeric(nb.bio.genes)){
        stop('The "nb.bio.genes" must be a positve numeric value 0 < nb.bio.genes < 1.')
    } else if (nb.bio.genes >= 1 | nb.bio.genes <= 0){
        stop('The "nb.bio.genes" must be a positve value 0 < nb.bio.genes < 1.')
    } else if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    } else if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(save.se.obj))) {
        stop('The "save.se.obj" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(use.imf))) {
        stop('The "use.imf" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(save.imf))) {
        stop('The "save.imf" must be "TRUE" or "FALSE".')
    }

    if(isTRUE(apply.log)){
        if(length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
            stop('The "pseudo.count" must be 0 or a postive integer value.')
    }

    if (is.null(assess.se.obj)) {
        if (isTRUE(sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables))) {
            stop('All or some of "bio.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (!is.null(variables.to.assess.bio.genes)) {
            if (isTRUE(sum(variables.to.assess.bio.genes %in% colnames(colData(se.obj))) != length(variables.to.assess.bio.genes)))
                stop('All or some of "variables.to.assess.bio.genes" cannot be found in the SummarizedExperiment object.')
        }
    }

    # check inputs ####
    if(!approach %in% c('AnovaCorr.PerBatchPerBio', 'AnovaCorr.AcrossAllSamples', 'TwoWayAnova', 'mad.unsupervised')){
        stop('The approach must be one of the "AnovaCorr.PerBatchPerBio", "AnovaCorr.AcrossAllSamples", "TwoWayAnova" or "mad.unsupervised".')
    }
    # Check the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, variables.to.assess.bio.genes)),
            remove.na = remove.na,
            verbose = verbose)
    }

    # Two-way anova approach ####
    if(approach == 'TwoWayAnova'){
        if (isFALSE(use.imf)){
            ## Data transformation ####
            printColoredMessage(
                message = '-- Data transformation:',
                color = 'magenta',
                verbose = verbose)
            ### apply log ####
            if (isTRUE(apply.log) & !is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Applying log2 + ', pseudo.count, ' (pseudo.count) on the ', assay.name, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Applying log2 on the ', assay.name, ' data.'),
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

            ## Create all possible homogeneous sample groups ####
            ### biological groups ####
            printColoredMessage(
                message = '-- Create all possible major homogeneous biological groups:',
                color = 'magenta',
                verbose = verbose)
            all.bio.groups <- createHomogeneousBioGroups(
                se.obj = se.obj,
                bio.variables = bio.variables,
                nb.clusters = nb.bio.clusters,
                clustering.method = bio.clustering.method,
                assess.variables = assess.variables,
                cat.cor.coef = cat.cor.coef,
                cont.cor.coef = cont.cor.coef,
                assess.se.obj = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)

            ### unwanted groups ####
            printColoredMessage(
                message = '-- Create all possible major groups with respect to sources of unwanted variation:',
                color = 'magenta',
                verbose = verbose)
            all.uv.groups <- createHomogeneousUVGroups(
                se.obj = se.obj,
                uv.variables = uv.variables,
                nb.clusters = nb.uv.clusters,
                clustering.method = uv.clustering.method,
                assess.se.obj = FALSE,
                assess.variables = assess.variables,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)

            ## Two_way ANOVA ####
            printColoredMessage(
                message = '-- Two_way ANOVA:',
                color = 'magenta',
                verbose = verbose)
            ## Perform gene level two_way ANOVA ####
                printColoredMessage(
                    message = '- Perform Two_way ANOVA:',
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = paste0('- This is between all individual gene-level expression',
                                     'and considering both biological and unwanted variables created above as factors.'),
                    color = 'blue',
                    verbose = verbose)
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
                ## rank the F-statistics ####
                bio.rank <- uv.rank <- Biology <- UV <- NULL
                set.seed(2190)
                all.aov$bio.rank <- rank(x = -all.aov$Biology, ties.method = 'random')
        }
        ## Read the intermediate file ####
        if (isTRUE(use.imf)){
            printColoredMessage(
                message = paste0('- Retrieve the results of two-way ANOVA from the the SummarizedExperiment object .'),
                color = 'blue',
                verbose = verbose)
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|TwoWayAnova|')
            }
            if(is.null(se.obj@metadata$IMF$BioGenes[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.aov <- se.obj@metadata$IMF$NCG[[imf.name]]
        }
        ## Save the intermediate file ####
        if(isTRUE(save.imf)){
            printColoredMessage(
                message = '-- Save a intermediate file:',
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = '- The results of two-way ANOVA is saved in the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose)
            if(length(se.obj@metadata$IMF) == 0 ) {
                se.obj@metadata[['IMF']] <- list()
            }
            if(!'BioGenes' %in% names(se.obj@metadata[['IMF']])){
                se.obj@metadata[['IMF']][['BioGenes']] <- list()
            }
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|TwoWayAnova|')
            }
            if(!imf.name %in% names(se.obj@metadata[['IMF']][['BioGenes']])){
                se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- all.aov
        }
        selected.bio.genes <- round(x = nb.bio.genes * nrow(se.obj) , digits = 0)
        selected.bio.genes <- row.names(all.aov)[all.aov$bio.rank < selected.bio.genes]
        selected.bio.genes <- row.names(se.obj) %in% selected.bio.genes
    }

    # Per biology and per batch ####
    if(approach == 'AnovaCorr.PerBatchPerBio'){
        if(isFALSE(use.imf)){
            ## Data transformation and normalization ####
            printColoredMessage(message = '-- Data transformation and normalization:',
                                color = 'magenta',
                                verbose = verbose)
            ### apply log ####
            if (isTRUE(apply.log) & !is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- apply log2 + ', pseudo.count, ' (pseudo.count) on the ', assay.name,' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- Apply log2 on the ', assay.name,' data.'),
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

            ### normalization ####
            if (!is.null(normalization)) {
                printColoredMessage(
                    message = '-- Data normalization:',
                    color = 'magenta',
                    verbose = verbose)
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
            ## Regress out variables ####
            ## regress out unwanted variables ####
            if(!is.null(regress.out.uv.variables)){
                printColoredMessage(
                    message = '- Regress out unwanted variables:',
                    color = 'blue',
                    verbose = verbose)
                if(!is.null(normalization)){
                    expr.data.reg.uv <- expr.data.nor
                } else expr.data.reg.uv <- expr.data
                printColoredMessage(
                    message = paste0('The ', paste0(regress.out.uv.variables, collapse = ' & '),
                                     ' will be regressed out from the data,',
                                     ' please make sure your data is log transformed.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = paste0(
                        'We do not recommend regressing out ',
                        paste0(regress.out.uv.variables, collapse = ' & '),
                        ' if they are largely associated with the ',
                        paste0(bio.variables, collapse = ' & '),
                        ' variables.'),
                    color = 'red',
                    verbose = verbose)
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
            ## Statistical analyses ####
            printColoredMessage(
                message = '-- Find biological genes:',
                color = 'magenta',
                verbose = verbose)

            ## select genes that are highly affected by biological variation  ####
            printColoredMessage(
                message = '-- Select genes that are highly affected by each source(s) of biological variation:',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(
                message = '- This step with be performed within each possible homogeneous unwanted groups.',
                color = 'blue',
                verbose = verbose)
            ### create all possible homogeneous uv groups ####
            printColoredMessage(
                message = '- reate all possible major groups with respect to sources of unwanted variation:',
                color = 'blue',
                verbose = verbose)
            if(is.null(uv.groups)){
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
                    assess.variables = assess.variables,
                    save.se.obj = FALSE,
                    verbose = verbose)
            } else if(!is.null(uv.groups)){
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(uv.groups, collapse = ' & '),
                        ' variables will be used as a major sources of unwanted variation',
                        ' to find all possible groups.'),
                    color = 'blue',
                    verbose = verbose)
                all.uv.groups <- createHomogeneousUVGroups(
                    se.obj = se.obj,
                    uv.variables = uv.groups,
                    nb.clusters = nb.uv.clusters,
                    clustering.method = uv.clustering.method,
                    assess.se.obj = FALSE,
                    assess.variables = assess.variables,
                    save.se.obj = FALSE,
                    verbose = verbose)
            }
            ### correlation between gene expression and all continuous source of biological variation with each uv groups ####
            bio.var.class <- unlist(lapply(
                bio.variables,
                function(x) class(colData(se.obj)[[x]])))
            continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
            if (length(continuous.bio) > 0) {
                printColoredMessage(message = '-- Correlation analyses:',
                                    color = 'magenta',
                                    verbose = verbose)
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
                                    corr.genes <-
                                        as.data.frame(correls(
                                            y = se.obj@colData[, x][selected.samples],
                                            x = t(data.to.use[, selected.samples]),
                                            type = corr.method,
                                            a = a ,
                                            rho = rho
                                        ))
                                    corr.genes <- cbind(round(x = corr.genes[, 1:4], digits = 3),
                                                        corr.genes[, 5, drop = FALSE])
                                    set.seed(2233)
                                    corr.genes$ranked.genes <- rank(-abs(corr.genes[, 'correlation']), ties.method = 'random')
                                    row.names(corr.genes) <- row.names(data.to.use)
                                    corr.genes
                                })
                            names(all.corr) <- selected.uv.groups
                            all.corr
                        })
                    names(corr.genes.bio) <- continuous.bio
                } else if (length(selected.uv.groups) == 0) {
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
                printColoredMessage(message = '-- ANOV analyses:',
                                    color = 'magenta',
                                    verbose = verbose)
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
                                message = paste0('There are not homogeneous unwanted groups that have at least ',
                                                 min.sample.for.aov , ' (min.sample.for.aov) samples within each batches of the ',
                                                 x,' variable. This may result in unsatisfactory NCG selection.'),
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
                                anova.gene.bio <- as.data.frame(
                                    row_oneway_equalvar(
                                        x = data.to.use[, selected.samples],
                                        g = se.obj@colData[, x][selected.samples]))
                                set.seed(2233)
                                anova.gene.bio$ranked.genes <- rank(-anova.gene.bio[, 'statistic'], ties.method = 'random')
                                anova.gene.bio
                            })
                        names(all.anova) <- selected.uv.groups
                        all.anova
                    })
                names(anova.genes.bio) <- categorical.bio
            } else anova.genes.bio <- NULL
        }

        ### Intermediate file ####
        ## read intermediate file ####
        if (isTRUE(use.imf)){
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|PerBiologyPerBatch|')
            }
            if(is.null(se.obj@metadata$IMF$BioGenes[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.tests <- se.obj@metadata$IMF$BioGenes[[imf.name]]
            anova.genes.bio <- all.tests$anova.genes.bio
            corr.genes.bio <- all.tests$corr.genes.bio
        }

        ## save intermediate file ####
        if(isTRUE(save.imf)){
            if(length(se.obj@metadata$IMF) == 0 ) {
                se.obj@metadata[['IMF']] <- list()
            }
            if(!'BioGenes' %in% names(se.obj@metadata[['IMF']])){
                se.obj@metadata[['IMF']][['BioGenes']] <- list()
            }
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|PerBiologyPerBatch|')
            }
            if(!imf.name %in% names(se.obj@metadata[['IMF']][['BioGenes']])){
                se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list(
                anova.genes.bio = anova.genes.bio,
                corr.genes.bio = corr.genes.bio)
        }
        ### select genes affected by biological variation ####
        selected.bio.genes <- round(x = nb.bio.genes * nrow(se.obj) , digits = 0)
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        selected.bio.genes <- unique(unlist(lapply(
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
                            row.names(se.obj)[all.ranks < selected.bio.genes]
                        })))
                }
            })))
        selected.bio.genes <- row.names(se.obj) %in% selected.bio.genes

    }
    # Across all samples ####
    if( approach == 'AnovaCorr.AcrossAllSamples'){
        ## Data transformation and normalization ####
        if(isFALSE(use.imf)){
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

            ## Gene-level ANOVA and correlation analyses  ####
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
                        anova.genes$ranked.genes <- rank(-anova.genes[, 'statistic'], ties.method = 'random')
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
                            x = -abs(corr.genes.var[, 'statistic']),
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
                imf.name <- paste0(assay.name, '|AcrossSamples|')
            }
            if(is.null(se.obj@metadata$IMF$BioGenes[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.tests <- se.obj@metadata$IMF$BioGenes[[imf.name]]
            anova.genes.bio <- all.tests$anova.genes.bio
            corr.genes.bio <- all.tests$corr.genes.bio
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
            if(!'BioGenes' %in% names(se.obj@metadata[['IMF']])){
                se.obj@metadata[['IMF']][['BioGenes']] <- list()
            }
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|AcrossSamples|')
            }
            if(!imf.name %in% names(se.obj@metadata[['IMF']][['BioGenes']])){
                se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list(
                anova.genes.bio = anova.genes.bio,
                corr.genes.bio = corr.genes.bio)
        }
        ### select genes affected by biological variation ####
        selected.bio.genes <- round(x = nb.bio.genes * nrow(se.obj) , digits = 0)
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        selected.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < selected.bio.genes
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        selected.bio.genes <- row.names(se.obj) %in% selected.bio.genes
    }
    # Unsupervised approach ####
    if(approach == 'mad.unsupervised'){
        ## Finding negative control genes ####
        if(isFALSE(use.imf)){
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
            printColoredMessage(
                message = '-- Find negative control genes',
                color = 'magenta',
                verbose = verbose)
            ## find genes that are not highly affected by biology ####
            if (!is.null(normalization)) {
                data.to.use <- expr.data.nor
            } else data.to.use <- expr.data

            #### regress out variables ####
            if(!is.null(regress.out.uv.variables)){
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(regress.out.uv.variables, collapse = ' & '),
                        ' variables will be regressed out from the data,',
                        ' please make sure your data is log transformed.'),
                    color = 'blue',
                    verbose = verbose)
                data.to.use <- t(data.to.use)
                lm.formula <- paste('se.obj', regress.out.uv.variables, sep = '$')
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
                clustering.method =  uv.clustering.method,
                nb.clusters = nb.uv.clusters,
                assess.variables = assess.variables,
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
                bio.genes$bio.ranks <- rank(x = -bio.genes$bio.mad, ties.method = 'random')
            } else if (isTRUE(length(groups) == 0))
                stop('There is no any homogenous sample groups with enough samples to perform MAD.')
        }
        ### Intermediate file ####
        ## read intermediate file ####
        if (isTRUE(use.imf)){
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|un.supervised|')
            }
            if(is.null(se.obj@metadata$IMF$BioGenes[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.tests <- se.obj@metadata$IMF$BioGenes[[imf.name]]
            bio.genes <- all.tests$bio.genes
        }

        ## save intermediate file ####
        if(isTRUE(save.imf)){
            if(length(se.obj@metadata$IMF) == 0 ) {
                se.obj@metadata[['IMF']] <- list()
            }
            if(!'BioGenes' %in% names(se.obj@metadata[['IMF']])){
                se.obj@metadata[['IMF']][['BioGenes']] <- list()
            }
            if(is.null(imf.name)){
                imf.name <- paste0(assay.name, '|un.supervised|')
            }
            if(!imf.name %in% names(se.obj@metadata[['IMF']][['BioGenes']])){
                se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['BioGenes']][[imf.name]] <- list(
                bio.genes = bio.genes)
        }
        ### select genes affected by biological variation ####
        selected.bio.genes <- round(x = nb.bio.genes * nrow(se.obj) , digits = 0)
        selected.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks < selected.bio.genes]
    }
    printColoredMessage(
        message = paste0('- ', sum(selected.bio.genes), ' genes are selected as negative control genes.'),
        color = 'blue',
        verbose = verbose)

    # Performance assessment of selected NCG  ####
    if(isTRUE(assess.bio.genes)){
        printColoredMessage(
            message = '-- Assess the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = paste0('- Perform PCA on only selected genes as NCG.'),
            color = 'blue',
            verbose = verbose)
        if(is.null(variables.to.assess.bio.genes))
            variables.to.assess.bio.genes <- c(bio.variables, uv.variables)
        printColoredMessage(
            message = paste0(
                '- Explore the association of the first ',
                nb.pcs,
                '  PCs with the ',
                paste0(variables.to.assess.bio.genes, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        ### pca ####
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[selected.bio.genes, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u

        ## regression and vector correlations ####
        all.corr <- lapply(
            variables.to.assess.bio.genes,
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
                        1 - prod(1 - cca$cor^2)
                        })
                }
            })
        names(all.corr) <- variables.to.assess.bio.genes
        pca.ncg <- as.data.frame(do.call(cbind, all.corr))
        pcs <- Groups <- NULL
        pca.ncg['pcs'] <- c(1:nb.pcs)
        pca.ncg <- tidyr::pivot_longer(
            data = pca.ncg,
            -pcs,
            names_to = 'Groups',
            values_to = 'ls')
        pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
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
        print(pca.ncg)
    }
    # Save results ####
    ### add results to the SummarizedExperiment object ####
    if(is.null(output.name)){
        output.name <- paste0(
            sum(selected.bio.genes),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|TWAnova:',
            '|',
            assay.name)}

    if(isTRUE(save.se.obj)){
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$BioGenes) == 0 ) {
            se.obj@metadata[['BioGenes']] <- list()
        }
        if(!'supervised' %in% names(se.obj@metadata[['BioGenes']])){
            se.obj@metadata[['BioGenes']][['supervised']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['BioGenes']][['supervised']])){
            se.obj@metadata[['BioGenes']][['supervised']][[output.name]] <- list()
        }
        se.obj@metadata[['BioGenes']][['supervised']][[output.name]] <- selected.bio.genes
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    }
    ### export results as logical vector ####
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The set of NCGs is outpputed as a logical vector.',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose)
        return(selected.bio.genes)
    }

}

