#' Computes ANOVA between individual gene expression and a categorical variable.--Finalized

#' @author Ramyar Molania

#' @description
#' This function calculates the ANOVA between individual gene expression of the assay(s) in a SummarizedExperiment object
#' and a categorical variable as factor.

#' @details
#' ANOVA enables us to assess the effects of a given qualitative variable (which we call a factor) on gene expression
#' measurements across any set of groups (labeled by the levels of the factor) under study. We use ANOVA F-statistics
#' to summarize the effects of a qualitative source of unwanted variation (for example, batches) on the expression levels
#' of individual genes, where genes having large F-statistics are deemed to be affected by the unwanted variation.
#' We also use ANOVA tests (the aov() function in R) to assign P values to the association between tumor purity and
#' molecular subtypes.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object to compute the ANOVA. By default is set to 'all', which means all assays of the
#' SummarizedExperiment object will be selected.
#' @param variable Character. A character string indicating a column name in the sample annotation of the SummarizedExperiment
#' object that contains a categorical variable, such as experimental batches, etc.
#' @param method Character. A character string specifying the method to be used for the ANOVA analysis. The options are
#' "aov" and "welch". The default is set to "aov". The function applies a fast implementation of one-way ANOVA using the
#' row_oneway_equalvar() or the row_oneway_welch() functions from the matrixTests R package. The first function assumes
#' equal variance between the groups while the second function applied the Welch correction for the groups with unequal
#' variances.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of the "nb.top.genes" genes, as visual
#' examples, with the highest and lowest F-statistics across the specified variable. The default is set to 'FALSE'.
#' @param nb.top.genes Numeric. A numeric value specifying the number of genes to plot from the top or bottom of the ANOVA
#' F-statistics list. The default is set to 3.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing ANOVA. The
#' default is set to 'TRUE'.
#' @param pseudo.count Numeric. A numeric value representing a pseudo count to be added to all measurements before applying
#' the log transformation. The default is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is set to 'TRUE'.
#' This means the function will apply the 'checkSeObj' function.
#' @param remove.na Character. A character string specifying whether to eliminate missing values from 'assays', 'sample.annotation',
#' 'both', or 'none'. When 'assays' is chosen, genes with missing values will be omitted. If 'sample.annotation' is selected,
#' samples with NA or missing values for each 'variable' will be excluded. The default is 'both'.
#' @param apply.round Logical. Indicates whether to round the values of F-statistics. The default is set to 'TRUE'.
#' @param save.se.obj Logical. Indicates whether to save the results, ANOVA F-statistics, and p-values in the metadata
#' of the SummarizedExperiment object or to output these results as a list or vector. The default is set to 'TRUE'.
#' @param override.check Logical. When set to TRUE, the function checks whether ANOVA has already been computed for the
#' current parameters on the SummarizedExperiment object. If it has, the metric will not be recalculated. The default is
#' set to 'FALSE'.
#' @param verbose Logical. If TRUE, displays the messages of different steps of the function.

#' @return Either a SummarizedExperiment object containing the log2 F-statistics and p-values of ANOVA for the continuous
#' variable or a list of these results.


#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer %>%
#' @importFrom dplyr mutate
#' @import ggplot2
#' @export

computeGenesVariableAnova <- function(
        se.obj,
        assay.names = 'all',
        variable,
        method = 'aov',
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
    printColoredMessage(message = '------------The computeGenesVariableAnova function starts:',
                        color = 'white',
                        verbose = verbose)

    # Checking to override or not ####
    if (!is.logical(override.check)){
        stop('The "override.check" must be logical (TRUE or FALSE)')
    }
    if (isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'ANOVA',
            method = method,
            variable = variable,
            file.name = 'fstatistics.pvalues'
        )
        if (is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if (isTRUE(compute.metric)){
        # Check the inputs ####
        if (is.null(assay.names)) {
            stop('The "assay.names" cannot be null.')
        }
        if (is.null(variable)) {
            stop('The "variable" cannot be empty.')
        }
        if (class(se.obj@colData[[variable]]) %in% c('numeric', 'integer')) {
            stop('The "variable" must be a categorical variable.')
        }
        if (length(unique(se.obj@colData[, variable])) < 2) {
            stop('The "variable" contains only one level. ANOVA cannot be performed.')
        }
        if (!method %in% c('aov', 'welch')) {
            stop('The method must be one of the "aov" or "welch".')
        }
        if (!is.logical(apply.log)){
            stop('The "apply.log" must be logical (TRUE or FALSE)')
        }
        if (isTRUE(apply.log)){
            if (length(pseudo.count) > 1){
                stop('The "pseudo.count" must be a numeric postive value.')
            }
            if (pseudo.count < 0 )
                stop('The value of "pseudo.count" cannot be negative.')
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

        # Checking the assays ####
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
        if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
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

        # Data transformation ####
        printColoredMessage(
            message = '-- Applying log transformation on the data:',
            color = 'magenta',
            verbose = verbose
            )
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            pseudo.count = pseudo.count,
            assessment = 'ANOVA',
            verbose = verbose
        )

        # Applying ANOVA ####
        printColoredMessage(
            message = paste0(
                '-- Performing the ANOVA between individual genes expression of the assay(s) and the "',
                variable,
                '" variable.'),
            color = 'magenta',
            verbose = verbose
        )
        all.aov <- lapply(
            levels(assay.names),
            function(x) {
                # compute anova ####
                if (method == 'aov') {
                    printColoredMessage(
                        message = paste0(
                            '- Performing the ANOVA with equal variance for the "',
                            x,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    anova.genes.var <- row_oneway_equalvar(
                        x = all.assays[[x]],
                        g = se.obj@colData[[variable]]
                    )
                } else if (method == 'welch') {
                    printColoredMessage(
                        message = paste0(
                            '- Performing ANOVA with Welch correction for the "',
                            x,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    anova.genes.var <- row_oneway_welch(
                        x =  all.assays[[x]],
                        g = se.obj@colData[[variable]]
                    )
                }
                row.names(anova.genes.var) <- row.names(se.obj)

                # rounding the obtained anova statistic  to 2 digits ####
                if (isTRUE(apply.round)) {
                    printColoredMessage(
                        message = '-- Rounding the obtained anova statistic to 2 digits',
                        color = 'magenta',
                        verbose = verbose
                    )
                    anova.genes.var <- cbind(
                        round(anova.genes.var[, 1:9], digits = 3),
                        anova.genes.var[, 10, drop = FALSE])
                }
                pvalue <- NULL
                anova.genes.var <- anova.genes.var[ , c('pvalue', 'statistic')]
                # plot top genes ####
                if (isTRUE(plot.top.genes)) {
                    temp.anova <- anova.genes.var[order(anova.genes.var[, 'statistic'],decreasing = TRUE, na.last = TRUE) , ]
                    var <- NULL
                    p.high <- as.data.frame(t(all.assays[[x]][row.names(temp.anova)[c(1:nb.top.genes)],]))
                    p.high <- mutate(p.high , var = se.obj@colData[, variable])
                    p.high <- pivot_longer(
                        data = p.high,
                        cols =  -var,
                        names_to = 'genes',
                        values_to = 'expr')
                    p.high <- ggplot(p.high, aes(x = var, y = expr)) +
                        geom_boxplot() +
                        ylab(expression(Log[2] ~ 'gene expression')) +
                        xlab(variable) +
                        facet_wrap( ~ genes) +
                        ggtitle(paste0(nb.top.genes," Top affected genes by the variable ", variable, " for ", x)) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', linewidth = 1),
                              axis.title.x = element_text(size = 14),
                              axis.title.y = element_text(size = 14),
                              axis.text.x = element_text(
                                  size = 10,
                                  angle = 45,
                                  hjust = 1,
                                  vjust = 1),
                              axis.text.y = element_text(size = 12),
                              legend.text = element_text(size = 10),
                              legend.title = element_text(size = 14),
                              strip.text.x = element_text(size = 10),
                              plot.title = element_text(size = 12))
                    print(p.high)
                    rm(temp.anova)
                }
                return(anova.genes.var)
            })
        names(all.aov) <- levels(assay.names)

        # Save the results ####
        printColoredMessage(
            message = '-- Saving the ANOVA results :',
            color = 'magenta',
            verbose = verbose
            )
        ## add results to the SummarizedExperiment object ####
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- The ANOVA results for the indiviaul assay(s) are saved to the "metadata" of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
                )
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'gene.level',
                assessment = 'ANOVA',
                method = method,
                variables = variable,
                file.name = 'fstatistics.pvalues',
                results.data = all.aov
            )
            printColoredMessage(
                message = paste0(
                    '- The ANOVA results of all the assays are saved to the ',
                    ' "se.obj@metadata$metric$AssayName$ANOVA$',
                    method,
                    ' in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
            printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                                color = 'white',
                                verbose = verbose
            )
            return(se.obj = se.obj)
        }
        ## return the results as a list ####
        if (isFALSE(save.se.obj)) {
            printColoredMessage(
                message = '- The ANOVA results for indiviaul assay are saved as list.',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                                color = 'white',
                                verbose = verbose
                                )
            return(genes.var.anova = all.aov)
        }
    } else {
        printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose
        )
        return(se.obj = se.obj)
    }
}
