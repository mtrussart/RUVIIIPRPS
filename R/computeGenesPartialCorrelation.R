#' Compute gene-gene partial correlation.

#' @author Ramyar Molania

#' @description
#' This function computes all possible gene-gene pairwise ordinary and partial correlation of the assays in the
#' SummarizedExperiment object.

#' @details
#' Partial correlation is used to estimate correlation between two variables while controlling for third
#' variables.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartile ranges. The default is set to "all, which
#' indicates all the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol that indicates the column name of the SummarizedExperiment object that contains a
#' continuous variable such as library size, tumor purity, ... .
#' @param method Symbol. A symbol that indicates which correlation methods should be used. The options are 'pearson',
#' 'kendall', or "spearman". The default is set to 'spearman'.
#' @param genes Vector. A vector comprising genes upon which the correlation analysis will be performed. The vector can be
#' logical, numeric or gene names. The default is 'NULL'.
#' @param select.genes Logical. If 'TRUE', the function will compute correlation between individual genes and the specified
#' variable, then select a subset of gene based the "corr.coff.cutoff" for down-stream analysis. The default is set to
#' 'TRUE'.
#' @param reference.data Symbol. A symbol specifying the name of the assay to be used to select genes using the correlation
#' analysis. The default is NULL. This means all the specified assays will be used.
#' @param corr.coff.cutoff Numeric. A numeric value to be used as a cutoff for selecting the genes. The default is set to
#' 0.7.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param apply.round Logical. Indicates whether to round the correlations coefficients or not. The default is set to 'TRUE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param override.check Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
#' determine if the PPcorr has already been computed for the current parameters. If it has, the metric will not be recalculated.
#' The default is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object or
#' to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return The SummarizedExperiment object that contains the correlation coefficients or a list of the correlation coefficient
#' for individual assays.

#' @importFrom Rfast correls
#' @importFrom stats cor
#' @export

computeGenesPartialCorrelation <- function(
        se.obj ,
        assay.names ,
        variable,
        method = 'spearman',
        genes = NULL,
        select.genes = TRUE,
        reference.data = NULL,
        corr.coff.cutoff = 0.3,
        apply.log = TRUE,
        pseudo.count = 1,
        apply.round = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        override.check = FALSE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeGenesPartialCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if(isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'PPcorr',
            method = method,
            variable = variable,
            file.name = 'correlations'
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
        if (is.list(assay.names)) {
            stop('The "assay.names" cannot be a list.')
        } else if (is.null(variable)) {
            stop('The "variable" cannot be NULL or empty.')
        } else if (length(variable) > 1) {
            stop('The "variable" must be the name of a variable in the SummarizedExperiment object.')
        } else if (!variable %in% colnames(colData(se.obj))) {
            stop('The "variable" cannot be found in the SummarizedExperiment object.')
        }
        if(isTRUE(select.genes)){
            if(is.null(corr.coff.cutoff)){
                stop('The "corr.coff.cutoff" cannot be NULL when the "select.genes" is set to "TRUE".')
            } else if (length(corr.coff.cutoff) > 1){
                stop('The "corr.coff.cutoff" must be a postive nmeric value between 0 and 1.')
            } else if (corr.coff.cutoff <= 0 | corr.coff.cutoff >= 1 ){
                stop('The "corr.coff.cutoff" must be a postive nmeric value between 0 and 1.')
            }
        }
        if(!is.null(genes) & isTRUE(select.genes)){
            stop('- Both "genes" and "select.genes" are specified, Please specifiy only one of them.')
        }
        if(!is.null(genes)){
            if(is.logical(genes)){
                if(length(genes) > nrow(se.obj)){
                    stop('The length of "genes" is larger than the total of genes in the SummarizedExperiment object.')
                } else if(sum(genes) == 0){
                    stop('All "genes" are "FALSE"')
                }
            } else if (is.character(genes)){
                if(sum(row.names(se.obj) %in% genes) != length(genes)){
                    stop('All or some of "genes" cannot be found in the SummarizedExperiment object.')
                }
            } else if (is.numeric(genes)){
                if(max(genes) > nrow(se.obj) ){
                    stop('The maximum value in "genes" is larger than the total of genes in the SummarizedExperiment object')
                } else if (sum(genes < 1) !=0 ){
                    stop('There are negative values in "genes".')
                }
            }
        }

        # Check the assays ####
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
                verbose = verbose)
        }

        # Data log transformation ####
        printColoredMessage(
            message ='-- Data log transformation:',
            color = 'magenta',
            verbose = verbose)
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assessment = 'partail gene-gene correlation',
            verbose = verbose
        )
        # Select genes ####
        printColoredMessage(
            message ='-- Select genes based the current parameters:',
            color = 'magenta',
            verbose = verbose
            )

        if(is.null(genes) & isFALSE(select.genes)){
            printColoredMessage(
                message = '- Both "genes" and "select.genes" are not specified, so all genes will be selected.',
                color = 'blue',
                verbose = verbose)
            selected.genes <- row.names(se.obj)
        }
        if(isTRUE(select.genes) & is.null(genes)){
            printColoredMessage(
                message = paste0('- Select genes that are highly correlated with the "', variable, '" variable.'),
                color = 'orange',
                verbose = verbose
            )
            if(!is.null(reference.data)){
                printColoredMessage(
                    message = paste0(
                        '* Find genes that have absolute correlations more than ',
                        corr.coff.cutoff, ' in the "', reference.data, '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                corr.genes.var <- correls(
                    y = se.obj@colData[[variable]],
                    x = t(all.assays[[reference.data]]) ,
                    type = method)[ , 'correlation', drop = FALSE]
                selected.genes <- row.names(corr.genes.var)[abs(corr.genes.var[ ,'correlation']) > corr.coff.cutoff]
                if(length(selected.genes) == 0){
                    stop('Any genes cannot be found based on the current "corr.coff.cutoff".')
                }
                printColoredMessage(
                    message = paste0(
                        '- ', length(selected.genes), ' genes are selected for pairwise-partial correlation analysis.' ),
                    color = 'blue',
                    verbose = verbose
                    )
            }
            if(is.null(reference.data)){
                printColoredMessage(
                    message = paste0(
                        '* Find genes that have absolute correlations more than ', corr.coff.cutoff, ' in each datasets.'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.corr.genes.var <- lapply(
                    levels(assay.names),
                    function(x){
                        corr.genes.var <- correls(
                            y = se.obj@colData[[variable]],
                            x = t(all.assays[[x]]) ,
                            type = method)[ , 'correlation', drop = FALSE]
                        selected.genes <- row.names(corr.genes.var)[abs(corr.genes.var[,'correlation']) > corr.coff.cutoff]
                    })
                selected.genes <- unique(unlist(all.corr.genes.var))
                if(length(selected.genes) == 0){
                    stop('Any genes cannot be found based on the current "corr.coff.cutoff".')
                }
                printColoredMessage(
                    message = paste0(
                        '- ', length(selected.genes), ' genes are selected for pairwise-partial correlation analysis.' ),
                    color = 'blue',
                    verbose = verbose
                    )
            }
        }
        if(!is.null(genes) & isFALSE(select.genes)){
            selected.genes <- genes
        }

        ## update all assays ####
        all.assays <- lapply(
            levels(assay.names),
            function(x) all.assays[[x]][ selected.genes, ])
        names(all.assays) <- levels(assay.names)


        # Compute gene-gene correlation ####
        printColoredMessage(
            message ='-- Compute all possible pairwise gene-gene correlations:',
            color = 'magenta',
            verbose = verbose
            )
        gene.gene.correlation <- lapply(
            levels(assay.names),
            function(x){
                printColoredMessage(
                    message = paste0('- Compute all gene-gene correlations for the "', x , '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                cor.matrix <- stats::cor(x = t(all.assays[[x]]) , method = method, use = "everything")
                upper.tri <- upper.tri(cor.matrix)
                variable.pairs <- which(upper.tri, arr.ind = TRUE)
                correlation.coefficients <- cor.matrix[upper.tri]
                if(isTRUE(apply.round))
                    correlation.coefficients <- round(x = correlation.coefficients, digits = 2)
                correlation.df <- data.frame(
                    gene1 = rownames(cor.matrix)[variable.pairs[, 1]],
                    gene2 = rownames(cor.matrix)[variable.pairs[, 2]],
                    p.cor = correlation.coefficients
                    )
            })
        names(gene.gene.correlation) <- levels(assay.names)

        # Compute gene-variable correlation ####
        printColoredMessage(
            message = paste0('-- Compute correlation between each genes with the "', variable, '" variable:'),
            color = 'magenta',
            verbose = verbose
            )
        gene.variable.correlation <- lapply(
            levels(assay.names),
            function(x){
                printColoredMessage(
                    message = paste0('- Compute the correlation for the "', x,  '" data:'),
                    color = 'blue',
                    verbose = verbose)
                corr.genes.var <- correls(
                    y = se.obj@colData[[variable]],
                    x = t(all.assays[[x]]) ,
                    type = method)[ , 'correlation', drop = FALSE]
                if(isTRUE(apply.round))
                    corr.genes.var[,1] <- round(x = corr.genes.var[,1], digits = 2)
                corr.genes.var
            })
        names(gene.variable.correlation) <-  levels(assay.names)

        # Compute gene-gene partial correlation ####
        printColoredMessage(
            message = '-- Compute all possible partial pairwise gene-gene correlations:',
            color = 'magenta',
            verbose = verbose
            )
        gene.gene.partial.correlation <- lapply(
            levels(assay.names),
            function(x){
                printColoredMessage(
                    message = paste0('- Compute the correlation for the "', x,  '" data:'),
                    color = 'blue',
                    verbose = verbose
                    )
                correlation.df <- gene.gene.correlation[[x]]
                index.1 <- match(correlation.df$gene1 , row.names(gene.variable.correlation[[x]]))
                index.2 <- match(correlation.df$gene2 , row.names(gene.variable.correlation[[x]]))
                correlation.df$gene1.corr <- gene.variable.correlation[[x]][,1][index.1]
                correlation.df$gene2.corr <- gene.variable.correlation[[x]][,1][index.2]
                b <- correlation.df$p.cor - c(correlation.df$gene1.corr * correlation.df$gene2.corr)
                d <- c(1 - correlation.df$gene1.corr^2) * c(1 - correlation.df$gene2.corr^2)
                correlation.df$pp.cor <- NULL
                correlation.df$pp.cor <- b/sqrt(d)
                if(isTRUE(apply.round))
                    correlation.df$pp.cor <- round(x = correlation.df$pp.cor, digits = 2)
                return(correlation.df)
            })
        names(gene.gene.partial.correlation) <- levels(assay.names)

        # Save the data ####
        ## add results to the SummarizedExperiment object ####
        printColoredMessage(
            message = '-- Save all the correlation coefficients data:',
            color = 'magenta',
            verbose = verbose
            )
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- Save all the correlation data to the "metadata" of the SummarizedExperiment object.',
                color = 'orange',
                verbose = verbose
                )
            ### for all assays
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'correlations',
                results.data = gene.gene.partial.correlation
            )
            printColoredMessage(
                message = paste0('* The correlations data of individual assay(s) is saved to the metadata@metric in SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(message = '------------The computeGenesPartialCorrelation function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj = se.obj)
        }
        ## output results as a list ####
        if(isFALSE(save.se.obj)){
            printColoredMessage(
                message = paste0('- The correlations data is outputed as a list.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(message = '------------The computeGenesPartialCorrelation function finished.',
                                color = 'white',
                                verbose = verbose)
            return(gene.gene.correlation = gene.gene.partial.correlation)
        }
    } else {
        printColoredMessage(message = '------------The computeGenesPartialCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
}


RamyarPcorr <- function(x, y, z){
    rxy <- cor.test(x , y)[[4]]
    rxz <- cor.test(x , z)[[4]]
    ryz <- cor.test(y , z)[[4]]
    rxyz.a <- rxy - c(rxz *ryz)
    rxyz.b <- c(1-rxz^2)*c(1-ryz^2)
    ppcor <- rxyz.a/sqr(rxyz.b)
    return(ppcor)
}





