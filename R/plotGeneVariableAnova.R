#' Plots F-statistics obtained from ANOVA.

#' @author Ramyar Molania

#' @description
#' This functions plots the results of the ANOVA analysis computed by the 'computeGenesVariableANOVA' function.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or vector. A character string or a vector of strings for selecting the name(s) of the
#' assay(s) in the SummarizedExperiment object. By default, all assays in the SummarizedExperiment object will be selected,
#' If the specifed ANOVA has already been computed for them using 'computeGenesVariableANOVA' function.
#' The ANOVA results for each specified assay will be plotted.
#' @param variable Character. Indicates the column name in the SummarizedExperiment object that contains a categorical
#' variable such as sample types or batches. The association between individual gene expression and the specified levels
#' must have already been computed by the 'computeGenesVariableANOVA' function.
#' @param anova.method Logical. Specifies which ANOVA method results should be used for plotting.
#' @param plot.ncol Numeric. A numeric value specifying the number of columns in the plot grid. The default is set to 4.
#' @param plot.nrow Numeric. A numeric value specifying the number of rows in the plot grid. The default is set to 1.
#' @param plot.output Logical. Indicates whether to display the boxplot(s) of the ANOVA results. By default, this
#' is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object
#' or to output the result as a list. By default, this is set to TRUE.
#' @param verbose Logical. If TRUE, displays process messages during execution.

#' @return Either a SummarizedExperiment object or a list containing the plots of ANOVA results for the specified
#' categorical variable for each assay.

#' @importFrom SummarizedExperiment assays assay
#' @import ggplot2
#' @export

plotGenesVariableAnova <- function(
        se.obj,
        assay.names = 'all',
        variable,
        anova.method = 'aov',
        plot.ncol = 4,
        plot.nrow = 1,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGenesVariableAnova function finished.',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }
    if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one variable.'))
    }
    if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable", must be a categorical variable.')
    }
    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain ANOVA F-statistics and p-values  ####
    printColoredMessage(
        message = paste0(
            '-- Obtaining the computed ANOVA F-statistics between genes and the "',
            variable,
            ' variable from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
        )
    all.aov.fvals.pvalues <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = assay.names,
        assessment = 'ANOVA',
        assessment.type = 'gene.level',
        method = anova.method,
        variables = variable,
        file.name = 'fstatistics.pvalues',
        sub.file.name = NULL,
        required.function = 'computeGenesPartialCorrelation',
        message.to.print = 'ANOVA'
    )
    # Generate different plots of F-statistic and p values of  ANOVA ####
    ## generate boxplots of the ANOVA F-statistics fo each assay ####
    printColoredMessage(
        message = '-- Generating boxplots of the lo2 of the ANOVA F-statistics:',
        color = 'magenta',
        verbose = verbose
    )
    aov.fvals <- NULL
    all.aov.fvals.boxplots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Generating a boxplot for the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            aov.fvals <- log2(all.aov.fvals.pvalues[[x]]$statistic)
            p.corr.coeff <- ggplot() +
                geom_boxplot(aes(y = aov.fvals, x = 1)) +
                ggtitle(variable) +
                xlab(x) +
                ylab(expression(Log[2]~'F-statistic')) +
                geom_hline(yintercept = 0) +
                theme(panel.background = element_blank(),
                      axis.line = element_line(colour = 'black', linewidth = 1),
                      axis.title.x = element_text(size = 18),
                      axis.title.y = element_text(size = 18),
                      plot.title = element_text(size = 15),
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_text(size = 0),
                      axis.text.y = element_text(size = 12))
        })
    names(all.aov.fvals.boxplots) <- levels(assay.names)

    ## put all boxplots together ####
    everything <- datasets <- NULL
    if (length(assay.names) > 1){
        printColoredMessage(
            message = '-- Putting all the boxplots of ANOVA F-statistics together.',
            color = 'magenta',
            verbose = verbose
        )
        all.aov.fvals <- sapply(levels(assay.names), function(x) log2(all.aov.fvals.pvalues[[x]][,2])) %>%
            data.frame(.) %>%
            tidyr::pivot_longer(
                everything(),
                names_to = 'datasets',
                values_to = 'aov.fvals')
        all.aov.fvals$datasets <- factor(
            x = all.aov.fvals$datasets,
            levels = levels(assay.names)
        )
        overall.aov.fvals.boxplots <- ggplot(all.aov.fvals, aes(x = datasets, y = aov.fvals)) +
            geom_boxplot(outlier.colour = 'gray') +
            ggtitle(variable) +
            xlab('Datasets') +
            ylab(expression(Log[2]~'F-statistics')) +
            geom_hline(yintercept = 0) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                plot.title = element_text(size = 16),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12))
        overall.aov.fvals.boxplots <- annotate_figure(
            p = overall.aov.fvals.boxplots,
            top = text_grob(
                label = "Analysis of variance (ANOVA)",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ',
                    'ANOA analysis between individaul gen expression an the  ',
                    variable,
                    ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
            )
        printColoredMessage(
            message = '- The boxplots of ANOVA F-statistics of each assays are combined into one plot.',
            color = 'blue',
            verbose = verbose
        )
        if (isTRUE(plot.output)) suppressMessages(print(overall.aov.fvals.boxplots))
    }

    # Generate p-value histograms of the ANOVA ####
    printColoredMessage(
        message = '-- Generating p-value histograms of the ANOVA F-statistics:',
        color = 'magenta',
        verbose = verbose
    )
    ## specified ylim ####
    breaks <- seq(from = 0, to = 1, by = .1)
    ylim.pvalue <- sapply(
        levels(assay.names),
        function(x) {
            binned <- cut(
                x = all.aov.fvals.pvalues[[x]][, 1],
                breaks = breaks,
                include.lowest = TRUE)
            frequency <- table(binned)[1]
        })
    ylim.pvalue <- ceiling(x = max(ylim.pvalue))

    ## generate histograms ####
    aov.fvals <- NULL
    all.aov.pvalues.histograms <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Generating boxplot for the "', x, '" data.'),
                color = 'blue',
                verbose = verbose
            )
            aov.all <- all.aov.fvals.pvalues[[x]]
            p.values.hist <- ggplot(data = aov.all, aes(x = pvalue)) +
                geom_histogram(binwidth = 0.1) +
                ggtitle(variable) +
                xlab('p-values') +
                ylab(expression('Frequency'~10^3)) +
                scale_y_continuous(labels = function(x) format(x /1000, scientific = F), limits = c(0, ylim.pvalue)) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    plot.title = element_text(size = 15),
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10),
                    axis.text.x = element_text(size = 8),
                    axis.text.y = element_text(size = 8))
        })
    names(all.aov.pvalues.histograms) <- levels(assay.names)

    ## put all histograms together ####
    if (length(assay.names) > 1){
        printColoredMessage(
            message = '-- Putting all the p-values histograms of the ANOVA together.',
            color = 'magenta',
            verbose = verbose
        )
        overall.aov.pvalues.histograms <- ggarrange(
            plotlist = all.aov.pvalues.histograms,
            ncol = plot.ncol,
            nrow = plot.nrow
            )
        if (class(overall.aov.pvalues.histograms)[[1]] == 'list'){
            plot.list <- lapply(
                seq(length(overall.aov.pvalues.histograms)),
                function(x){
                    annotate_figure(
                        p = overall.aov.pvalues.histograms[[x]],
                        top = text_grob(
                            label = "Analysis of variance (ANOVA)",
                            color = "orange",
                            face = "bold",
                            size = 18),
                        bottom = text_grob(
                            label = paste0(
                                'Analysis: ',
                                'ANOVA analysis between individaul gen expression an the  ',
                                variable,
                                ' variable.'),
                            color = "black",
                            hjust = 1,
                            x = 1,
                            size = 10))
                })
            overall.aov.pvalues.histograms <- ggarrange(
                plotlist = plot.list,
                ncol = 1,
                nrow = 1
                )
        } else {
            overall.aov.pvalues.histograms <- annotate_figure(
                p = overall.aov.pvalues.histograms,
                top = text_grob(
                    label = "Analysis of variance (ANOVA)",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'ANOA analysis between individaul gen expression an the  ',
                        variable,
                        ' variable.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10))
        }
        printColoredMessage(
            message = '- The p-values histograms of each assays are combined into one plot.',
            color = 'blue',
            verbose = verbose
        )
        if(isTRUE(plot.output)) suppressMessages(print(overall.aov.pvalues.histograms))
    }

    # Save the results ####
    printColoredMessage(
        message = '-- Saving the all the boxplots of the ANOVA F-statistics:',
        color = 'magenta',
        verbose = verbose
        )
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the boxplots of the ANOVA F-statistics of each assay(s) in the "metadata" in the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'ANOVA',
            method = anova.method,
            variables = variable,
            file.name = 'boxplot',
            results.data = all.aov.fvals.boxplots
        )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'ANOVA',
            method = anova.method,
            variables = variable,
            file.name = 'histogram',
            results.data = all.aov.pvalues.histograms
        )
        printColoredMessage(
            message = paste0(
                '- The boxplot of the ANOVA F-statistics for the indiviaul assay(s) is  saved to',
                '  "se.obj@metadata$metric$AssayName$ANOVA$',
                anova.method,
                ' in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)

        if (length(assay.names) > 1) {
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'ANOVA',
                method = anova.method,
                variables = variable,
                file.name = 'boxplot',
                plot.data = overall.aov.fvals.boxplots
            )
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'ANOVA',
                method = anova.method,
                variables = variable,
                file.name = 'histogram',
                plot.data = overall.aov.pvalues.histograms
            )
            printColoredMessage(
                message = paste0('- The combined boxplots of the ANOVA F-statistics and p-values of all the  indiviaul assay(s) is saved to',
                                 '  "se.obj@metadata$plot$ANOVA$', anova.method, ' in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(message = '------------The plotGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

    }
    if (isFALSE(save.se.obj)) {
        # return only the correlation result ####
        printColoredMessage(
            message = '- All the plots are saved as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(list(
                all.aov.fvals.boxplots = all.aov.fvals.boxplots,
                all.aov.pvalues.histograms = all.aov.pvalues.histograms)
                )
        } else {
            return(list(
                all.aov.fvals.boxplots = all.aov.fvals.boxplots,
                all.aov.pvalues.histograms = all.aov.pvalues.histograms,
                overall.aov.fvals.boxplots = overall.aov.fvals.boxplots,
                overall.aov.pvalues.histograms = overall.aov.pvalues.histograms)
            )
        }
    }
}
