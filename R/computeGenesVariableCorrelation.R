#' Computes the correlation between individual gene expression and a continuous variable.

#' @author Ramyar Molania

#' @description
#' This function computes Spearman or Pearson correlations between individual gene-level expression of each assay and
#' a continuous variable in SummarizedExperiment object.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object for computing the correlation. By default is set to 'all', which indicates all
#' assays in the SummarizedExperiment object will be selected.
#' @param variable Character. A character string indicating a column name in the SummarizedExperiment object that contains
#' a continuous variable, such as library size, tumor purity, etc.
#' @param method Character. Specifies which correlation method should be used. The options are 'pearson', 'kendall', or
#' "spearman". The default is set to 'spearman'.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation. The default is 0.05.
#' Refer to the 'correls' function from the Rfast R package for more details.
#' @param rho Numeric. The hypothesized correlation value to be used in hypothesis testing. The default is 0.
#' Refer to the 'correls' function from the Rfast R package for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before computing the correlation.
#' By default, the log transformation is applied.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before log transformation.
#' The default is set to 1.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression for the top or bottom genes based on
#' correlation values. The default is set to 'FALSE'.
#' @param nb.top.genes Numeric. Defines the number of genes with the highest or lowest correlation to the variable to plot.
#' The default is 3.
#' @param apply.round Logical. Indicates whether to round the correlation coefficients. The default is TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object.
#' Refer to the 'checkSeObj' function for more details. The default is TRUE.
#' @param remove.na Character. Specifies whether to remove NA or missing values from the assays and variable.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the current SummarizedExperiment
#' object or to output the result. By default, it is set to TRUE.
#' @param override.check Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
#' determine if the correlation has already been computed with the current parameters. If so, the metric will not be recalculated.
#' The default is FALSE.
#' @param verbose Logical. If 'TRUE', displays the messages for different steps of the function.

#' @return Either a SummarizedExperiment object or a list containing the correlation coefficients for compuated for individual
#' genes.


#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom tidyr pivot_longer %>%
#' @importFrom Rfast correls
#' @importFrom stats var
#' @import ggplot2
#' @export

computeGenesVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        method = 'spearman',
        a = 0.05,
        rho = 0,
        apply.log = TRUE,
        pseudo.count = 1,
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        apply.round = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeGenesVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if (!is.logical(override.check)){
        stop('The "override.check" must be logical (TRUE or FALSE)')
    }
    if (isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'Correlation',
            method = method,
            variable = variable,
            file.name = 'correlations.pvalues'
        )
        if (is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if (isTRUE(compute.metric)){
        # Checking the inputs ####
        if (is.null(assay.names)) {
            stop('The "assay.names" cannot be empty.')
        }
        if (is.list(assay.names)){
            stop('The "assay.names" must be a vector of the assay names(s) or "assay.names = all".')
        }
        if (is.null(variable)) {
            stop('The "variable" cannot be empty.')
        }
        if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
            stop(paste0('The ', variable, 'should be a continuous variable.'))
        }
        if (!method %in% c('pearson', 'spearman')) {
            stop('The method must be one of the "pearson" or "spearman".')
        }
        if (var(se.obj[[variable]]) == 0) {
            stop(paste0('The ', variable, ' appears to have no variation.'))
        }

        # Checking the assays ####
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
        if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }
        if (!is.logical(plot.top.genes)){
            stop('The "plot.top.genes" must be logical (TRUE or FALSE)')
        }
        if (!is.logical(apply.round)){
            stop('The "apply.round" must be logical (TRUE or FALSE)')
        }
        if (!is.logical(save.se.obj)){
            stop('The "save.se.obj" must be logical (TRUE or FALSE)')
        }
        if (!is.logical(verbose)){
            stop('The "verbose" must be logical (TRUE or FALSE)')
        }

        # Checking SummarizedExperiment object ####
        if (isTRUE(assess.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = variable,
                remove.na = remove.na,
                verbose = verbose
                )
        }

        # Data transformation ####
        printColoredMessage(
            message = '-- Data log transformation:',
            color = 'magenta',
            verbose = verbose)
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assessment = 'correlation',
            verbose = verbose
        )

        # Compute correlation analyses ####
        printColoredMessage(
            message = paste0(
                'Performing ' ,
                method,
                ' correlation between individual genes expression of the assay(s) and the "',
                variable,
                '" variable.'),
            color = 'magenta',
            verbose = verbose
            )
        all.correlations <- lapply(
            levels(assay.names),
            function(x) {
                # correlation ####
                printColoredMessage(
                    message = paste0(
                        '- Performing the correlation on the "',
                        x,
                        '" data and obtain the coefficients.'),
                    color = 'blue',
                    verbose = verbose
                    )
                corr.genes.var <- correls(
                    y = se.obj@colData[, variable],
                    x = t(all.assays[[x]]),
                    type = method,
                    a = a ,
                    rho = rho
                )
                row.names(corr.genes.var) <- row.names(se.obj)
                # round the correlation obtained to 2 digits ####
                if (isTRUE(apply.round)) {
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 2),
                        corr.genes.var[, 5, drop = FALSE])
                }
                # plot highly affected genes ####
                if (isTRUE(plot.top.genes)) {
                    printColoredMessage(
                        message = paste0(
                            '- Plotting top ' ,
                            nb.top.genes,
                            ' highly correlated (positive and negative) genes with the ',
                            variable,'.'),
                        color = 'blue',
                        verbose = verbose)
                    temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'], decreasing = TRUE, na.last = TRUE) , ]
                    ### positive correlation ####
                    p.pos <- as.data.frame(t(all.assays[[x]][row.names(temp.corr)[c(1:nb.top.genes)],]))
                    p.pos[, 'variable'] <- se.obj@colData[[variable]]
                    p.pos <- p.pos %>%
                        pivot_longer(
                            -variable,
                            names_to = 'genes',
                            values_to = 'expr')
                    p.pos <- ggplot(p.pos, aes(x = variable, y = expr)) +
                        geom_point() +
                        ylab(expression(Log[2] ~ 'gene expression')) +
                        xlab(variable) +
                        geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                        facet_wrap(~ genes) +
                        ggtitle(paste0("Top highly positively correlated genes with ", variable,'\n in the assay ', x)) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', size = 1),
                              axis.title.x = element_text(size = 14),
                              axis.title.y = element_text(size = 14),
                              axis.text.x = element_text(size = 10),
                              axis.text.y = element_text(size = 12),
                              strip.text.x = element_text(size = 10),
                              plot.title = element_text(size = 14))
                    print(p.pos)

                    # negative correlation ####
                    temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                                      decreasing = FALSE,
                                                      na.last = TRUE) ,]
                    p.neg <- as.data.frame(t(all.assays[[x]][row.names(temp.corr)[c(1:nb.top.genes)],]))
                    p.neg[, 'variable'] <- se.obj@colData[, variable]
                    p.neg <- p.neg %>%
                        tidyr::pivot_longer(
                            -variable,
                            names_to = 'genes',
                            values_to = 'expr')
                    p.neg <- ggplot(p.neg, aes(x = variable, y = expr)) +
                        geom_point() +
                        ylab(expression(Log[2] ~ 'gene expression')) +
                        xlab(variable) +
                        geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                        facet_wrap(~ genes) +
                        ggtitle(paste0('Top highly negatively correlated genes with the variable ', variable, '\n in the assay ', x)) +
                        theme( panel.background = element_blank(),
                               axis.line = element_line(colour = 'black', size = 1),
                               axis.title.x = element_text(size = 14),
                               axis.title.y = element_text(size = 14),
                               axis.text.x = element_text(size = 10),
                               axis.text.y = element_text(size = 12),
                               strip.text.x = element_text(size = 10),
                               plot.title = element_text(size = 14))

                    print(p.neg)
                    rm(temp.corr)
                }
                return(corr.genes.var[ , c('p-value', 'correlation')])
            })
        names(all.correlations) <- levels(assay.names)

        # Save the results ####
        printColoredMessage(
            message = '-- Saving the correlation results:',
            color = 'magenta',
            verbose = verbose
            )
        ## add results to the SummarizedExperiment object ####
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- The correlation results for the indiviaul assay(s) are saved to the "metadata" of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
                )
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'gene.level',
                assessment = 'Correlation',
                method = method,
                variables = variable,
                file.name = 'correlations.pvalues',
                results.data = all.correlations
            )
            printColoredMessage(
                message = paste0(
                    '- The correlation results of all the assays are saved to the ',
                    ' "se.obj@metadata$metric$AssayName$Correlation$',
                    method,
                    ' in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)

            printColoredMessage(message = '------------The computeGenesVariableCorrelation function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj = se.obj)
        }
        ## output the results as list ####
        if (isFALSE(save.se.obj)) {
            printColoredMessage(
                message = '-The correlation results for indiviaul assay are saved as list.',
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(message = '------------The computeGenesVariableCorrelation function finished.',
                color = 'white',
                verbose = verbose
                )
            return(all.correlations)
        }
    } else {
        printColoredMessage(message = '------------The computeGenesVariableCorrelation function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj = se.obj)
    }
}
