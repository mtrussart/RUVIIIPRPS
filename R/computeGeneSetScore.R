#' Compute sample wise score of a gene set.

#' @author Ramyar Molania

#' @description
#' The function uses the singscore function from the R/Bioconductor singscore to calculate sample wise score of
#' a gene set in for all the assay(s) in the SummarizedExperiment object. Refer to the details for more information.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol or a vector of symbols specifying the name(s) of the assay(s) in the SummarizedExperiment
#' object to calculate sample-wise scores. The default is set to "all," indicating all assays in the SummarizedExperiment
#' object will be selected.
#' @param upset.genes Vector. A character vector of gene names/ids that are up-regulated in the a gen set.
#' @param downset.genes Vector. A character vector of gene names/ids that are down-regulated in the a gen set.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param normalization Symbol. A symbol indicating the name of normalization to be used before computing the scores. The
#' default is set to 'NULL'. Refer to the 'applyOtherNormalizations' function for more details.
#' @param regress.out.variables Symbol or vector of symbols indicating the name(s) of the column(s) in the SummarizedExperiment
#' object to be regressed out from the data before computing the scores. The default is set to 'NULL'.
#' @param assess.score Logical. If 'TRUE', the association between the computed scores and specified variables will be
#' assessed. The default is set to 'NULL'. See the details for more information.
#' @param variables.to.assess Symbol. A symbol or vector of symbols indicating the name(s) of the column(s) in the SummarizedExperiment
#' object to be used for the assessment of the computed scores. It can comprise of continuous and categorical variables.
#' The default is set to 'NULL'.
#' @param corr.method Symbol. A symbol indicating which correlation method should be used for the correlation analysis of
#' the computed scores with the specified continuous variable(s). The options are "pearson", "kendall", "spearman". The
#' default is set to 'spearman'.
#' @param gene.set.name Symbol. A symbol indicating the name to be used to save the score in the SummarizedExperiment object.
#' if 'NULL', the function will select a name as follow :
#' gene.set.name <- paste0('singscore|',length(c(upset.genes, downset.genes)),'genes')
#' @param plot.output Logical. If 'TRUE', the assessment plot will be printed while running the function. The default is
#' set to 'TRUE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to 'TRUE'. Refer to the 'checkSeObj' function for more details.
#' @param save.se.obj Logical. Indicates whether to save the score results in the metadata of the SummarizedExperiment object
#' or to output the result as list. The default is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom singscore rankGenes simpleScore
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cor.test cancor
#' @import ggplot2
#' @export

computeGeneSetScore <- function(
        se.obj,
        assay.names = 'all',
        upset.genes,
        downset.genes = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        normalization = NULL,
        regress.out.variables = NULL,
        assess.score = FALSE,
        variables.to.assess = NULL,
        corr.method = 'spearman',
        gene.set.name = NULL,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeGeneSetScore function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    # Data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Data transformation and normalization:',
        color = 'magenta',
        verbose = verbose
        )
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            ## apply log ####
            if(is.null(normalization) & is.null(regress.out.variables)){
                do.log <- TRUE
            } else if (is.null(normalization) & !is.null(regress.out.variables)){
                do.log <- TRUE
            }
            if(isTRUE(do.log)){
                if (isTRUE(apply.log) & !is.null(pseudo.count)){
                    printColoredMessage(
                        message = paste0('- apply log2 + ', pseudo.count, ' (pseudo.count) on the ', x,' data.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- log2(assay(x = se.obj, i = x) + pseudo.count)
                } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                    printColoredMessage(
                        message = paste0('- Apply log2 on the "', x,'" data.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- log2(assay(x = se.obj, i = x))
                } else if (isFALSE(apply.log)) {
                    printColoredMessage(
                        message = paste0('The "', x, '" data will be used without any log transformation.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- assay(x = se.obj, i = x)
                }
            }
            ## normalization ####
            if (!is.null(normalization)) {
                printColoredMessage(
                    message = '-- Data normalization:',
                    color = 'magenta',
                    verbose = verbose)
                expr.data <- applyOtherNormalizations(
                    se.obj = se.obj,
                    assay.name = x,
                    method = normalization,
                    pseudo.count = pseudo.count,
                    apply.log = apply.log,
                    assess.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = verbose
                )
            }
            # Regress out variables ####
            if(!is.null(regress.out.variables)){
                printColoredMessage(
                    message = '-- Regress out unwanted or biological variables:',
                    color = 'magenta',
                    verbose = verbose
                )
                printColoredMessage(
                    message = paste0('The ', paste0(regress.out.variables, collapse = ' & '),
                                     ' will be regressed out from the data,',
                                     ' please make sure your data is log transformed.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = paste0(
                        'We do not recommend regressing out ',
                        paste0(regress.out.variables, collapse = ' & '),
                        ' if they are largely associated with the ',
                        paste0(regress.out.variables, collapse = ' & '),
                        ' variables.'),
                    color = 'red',
                    verbose = verbose)
                expr.data <- t(expr.data)
                uv.variables.all <- paste('se.obj', regress.out.variables, sep = '$')
                expr.data <- lm(as.formula(paste(
                    'expr.data',
                    paste0(uv.variables.all, collapse = '+') ,
                    sep = '~')))
                expr.data <- t(expr.data$residuals)
                colnames(expr.data) <- colnames(se.obj)
                row.names(expr.data) <- row.names(se.obj)
            }
            return(expr.data)
            }
        )
    names(all.assays) <- levels(assay.names)

    # Gene set scoring ####
    all.scores <- lapply(
        levels(assay.names),
        function(x){
            rank.data <- singscore::rankGenes(expreMatrix = all.assays[[x]])
            if(is.null(downset.genes)){
                gene.set.score <- singscore::simpleScore(
                    rankData = rank.data,
                    upSet = upset.genes
                )
            } else {
                gene.set.score <- singscore::simpleScore(
                    rankData = rank.data,
                    upSet = upset.genes,
                    downSet = downset.genes
                )
            }
            gene.set.score <- gene.set.score$TotalScore
            return(gene.set.score)
        })
    rm(all.assays)
    names(all.scores) <- levels(assay.names)

    # Assess the score ####
    if(isTRUE(assess.score)){
        class.variables <- sapply(
            variables.to.assess,
            function(x) class(se.obj[[x]]))
        continuous.variables <- variables.to.assess[class.variables %in% c('numeric', 'integer')]
        categorical.variables <- variables.to.assess[class.variables %in% c('factor', 'charachter')]

        all.assessment.plots <- lapply(
            levels(assay.names),
            function(x){
                ## continuous ####
                if(length(continuous.variables) > 1){
                    corr.continuous <- sapply(
                        continuous.variables,
                        function(y){
                            cor.test(x = all.scores[[x]], y = se.obj[[y]], method = corr.method)[[4]][[1]]
                        })
                }
                ## categorical ####
                if(length(categorical.variables) > 1){
                    corr.categorical <- sapply(
                        categorical.variables,
                        function(y){
                            catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca <- cancor(x = all.scores[[x]], y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                }
                ## plot #####
                corr <- variable.name <- NULL
                all.corr <- data.frame(
                    corr = c(corr.continuous, corr.categorical),
                    variable.name = names(c(corr.continuous, corr.categorical))
                )
                p.assess.geneset <- ggplot(data = all.corr, aes(x = variable.name, y = corr)) +
                    geom_col() +
                    ggtitle('Assessment of the gene set scoring') +
                    xlab('Variables') +
                    ylab('Correlations') +
                    ylim(c(-1,1)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        plot.title = element_text(size = 18),
                        axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                        axis.text.y = element_text(size = 12))
                if(isTRUE(plot.output)) print(p.assess.geneset)
                return(p.assess.geneset)
            })
        names(all.assessment.plots) <- levels(assay.names)
    }
    # Save the results ####
    if(isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Save all the gene set enrichment score into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if(is.null(gene.set.name)){
            gene.set.name <- paste0(
                'singscore|',
                length(c(upset.genes, downset.genes)),
                'genes')
        }
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            slot = 'Metrics',
            assessment = 'GeneSetScore',
            assessment.type = 'global.level',
            method = 'singscore',
            file.name = 'score',
            variables = gene.set.name,
            results.data = all.scores
        )
        if(isTRUE(assess.score)){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                assay.names = assay.names,
                slot = 'Metrics',
                assessment = 'GeneSetScore',
                assessment.type = 'global.level',
                method = 'singscore',
                file.name = 'assessment.plot',
                variables = gene.set.name,
                results.data = all.assessment.plots
            )
        }
        printColoredMessage(message = '------------The computeGeneSetScore function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## save the data as list ####
    if(isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The computeGeneSetScore function finished.',
                            color = 'white',
                            verbose = verbose)
        if(isTRUE(assess.score)){
            return(list(all.scores = all.scores, all.assessment.plots = all.assessment.plots))
        } else return(list(all.scores = all.scores))
    }
}
