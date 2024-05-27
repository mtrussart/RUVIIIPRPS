#' Compute gene-gene partial correlation.

#' @author Ramyar Molania

#' @description
#' This function computes all gene-gene pairwise correlation of the assays in the SummarizedExperiment object.


#' @details
#' This function first calculates thte or...
#'

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartile ranges. The default is set to "all, which
#' indicates all the assays of the SummarizedExperiment object will be selected.
#' @param genes vector.
#' @param variable Symbol. A symbol that indicates the column name of the SummarizedExperiment object that contains a
#' continuous variable such as library size, tumor purity, ... .
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param method Symbol. A symbol that indicates which correlation methods should be used. The options are 'pearson',
#' 'kendall', or "spearman". The default is set to 'spearman'.
#' @param apply.round Logical. Indicates whether to round the correlations coefficients or not. The default is set to 'TRUE'.
#' @param plot.output Logical. Indicates to print the histogram of partial correlation coefficients or not. The default is
#' set to 'FLASE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param output.name = Symbol.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object or
#' to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom Rfast correls
#' @importFrom stats cor
#' @import ggplot2
#' @export

computeGenePartialCorrelation <- function(
        se.obj ,
        assay.names ,
        genes = NULL,
        variable,
        apply.log = TRUE,
        pseudo.count = 1,
        method = 'spearman',
        apply.round = TRUE,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        output.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # Data log transformation ####
    printColoredMessage( message ='-- Data log transformation:',
                         color = 'magenta',
                         verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (isTRUE(apply.log) & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('- Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- Apply log2 on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr.data <- log2(assay(x = se.obj, i = x))
            } else if (isFALSE(apply.log)){
                printColoredMessage(
                    message = paste0('- The ', x, ' data will be used without log transformation.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = '- Please note, the assay should be in log scale before computing RLE.',
                    color = 'red',
                    verbose = verbose)
                expr.data <- assay(x = se.obj, i = x)
            }
        })
    names(all.assays) <- levels(assay.names)

    # gene-gene correlation
    gene.gene.correlation <- lapply(
        levels(assay.names),
        function(x){
            cor.matrix <- stats::cor(x = t(all.assays[[x]]) , method = method, use = "everything")
            upper.tri <- upper.tri(cor.matrix)
            variable.pairs <- which(upper.tri, arr.ind = TRUE)
            correlation.coefficients <- cor.matrix[upper.tri]
            if(isTRUE(apply.round))
                correlation.coefficients <- round(x = correlation.coefficients, digits = 2)
            correlation.df <- data.frame(
                gene1 = rownames(cor.matrix)[variable.pairs[, 1]],
                gene2 = rownames(cor.matrix)[variable.pairs[, 2]],
                p.cor = correlation.coefficients)
        })
    names(gene.gene.correlation) <-  levels(assay.names)

    # gene-gene correlation
    gene.variable.correlation <- lapply(
        levels(assay.names),
        function(x){
            corr.genes.var <- correls(
                y = se.obj@colData[, variable],
                x = t(all.assays[[x]]) ,
                type = method)[ , 'correlation', drop = FALSE]
            if(isTRUE(apply.round))
                corr.genes.var[,1] <- round(x = corr.genes.var[,1], digits = 2)
            corr.genes.var
        })
    names(gene.variable.correlation) <-  levels(assay.names)

    # gene-gene partial correlation
    gene.gene.partial.correlation <- lapply(
        levels(assay.names),
        function(x){
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

    # plotting
    ## each assay
    hists.plots <- lapply(
        levels(assay.names),
        function(x){
            p <- ggplot(gene.gene.partial.correlation[[x]], aes(x = pp.cor)) +
                geom_histogram() +
                xlab('Partial correlations') +
                ggtitle(x) +
                theme(
                    panel.background = element_blank(),
                    plot.title = element_text(size = 10),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10),
                    axis.text.y = element_text(size = 9),
                    axis.ticks.x = element_blank())
        })
    names(hists.plots) <- levels(assay.names)
    ## all assays
    all.hists.plots <- ggarrange(plotlist = hists.plots)

    # Save the data
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save the correlation coefficients data:',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save all the RLE data to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ### for all assays
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'PPcorr' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['PPcorr']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'correlation' %in% names(se.obj@metadata[['metric']][[x]][['PPcorr']])) {
                se.obj@metadata[['metric']][[x]][['PPcorr']][['correlation']] <- list()
            }
            se.obj@metadata[['metric']][[x]][['PPcorr']][['correlation']] <- gene.gene.partial.correlation[[x]]
            ## check if metadata metric already exist for this assay and this metric
            if (!'plot' %in% names(se.obj@metadata[['metric']][[x]][['PPcorr']])) {
                se.obj@metadata[['metric']][[x]][['PPcorr']][['plot']] <- list()
            }
            se.obj@metadata[['metric']][[x]][['PPcorr']][['plot']] <- hists.plots[[x]]
        }
        printColoredMessage(
            message = paste0('The RLE data of individual assay is saved to the metadata@metric in SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(gene.gene.partial.correlation)
    }

}






