#' Plots Spearman or Pearson correlations coefficients obtained from gene-level analysis.

#' @author Ramyar Molania

#' @description
#' This function generates boxplots of computed Spearman or Pearson correlation coefficients of individual assays in a
#' SummarizedExperiment object. The correlation coefficients must be computed by the 'computeGeneVariableCorrelation' function.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or vector. A character string or vector of strings for selecting the name(s) of the assay(s)
#' in the SummarizedExperiment for which gene-level correlation with the speciferd  for computing the correlation. By default, all assays in the SummarizedExperiment object will be selected.
#' @param variable Character. The name of the column in the sample annotation of the SummarizedExperiment object that contains a
#' continuous variable such as library size, tumor purity, etc.
#' @param correlation.method Character. Specifies which computed correlation coefficient should be used for plotting. The
#' default is 'gene.spearman.corr'. Refer to the 'computeGenesVariableCorrelation' function for more details.
#' @param plot.ncol Numeric. A numeric value specifying the number of columns in the plot grid. The default is set to 4.
#' @param plot.nrow Numeric. A numeric value specifying the number of rows in the plot grid. The default is set to 1.
#' @param plot.output Logical. Indicates whether to display the plot. By default, it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object
#' or to output the result. By default, this is set to TRUE.
#' @param verbose Logical. If TRUE, process messages will be displayed during the function's execution.


#' @return A SummarizedExperiment object or a list that containing the boxplots of the Spearman or Pearson correlations
#' coefficients for individual assays.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr everything
#' @import ggplot2
#' @export

plotGenesVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        correlation.method = 'spearman',
        plot.ncol = 4,
        plot.nrow = 1,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGenesVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    }
    if (is.null(variable)) {
        stop('Please provide a variable.')
    }
    if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one variable.'))
    }
    if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0(
            'The ',
            variable,
            ', is a numeric, but this should a categorical variable'
        ))
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain correlations coefficients ####
    printColoredMessage(
        message = paste0(
            '-- Obtaining the computed correlation coefficient between genes and the "',
            variable,
            '" variable from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose)

    all.corr.coeff.pvalues <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = levels(assay.names),
        assessment = 'Correlation',
        assessment.type = 'gene.level',
        method = correlation.method,
        variables = variable,
        file.name = 'correlations.pvalues',
        sub.file.name = NULL,
        required.function = 'computeGenesVariableCorrelation',
        message.to.print = 'correlations and pvalues'
        )

    # Generate boxplots of the correlation coefficients ####
    printColoredMessage(
        message = '-- Generating boxplots of the correlation coefficients:',
        color = 'magenta',
        verbose = verbose
        )
    all.corr.coeff.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Generating boxplot for the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            corr.coeff <- all.corr.coeff.pvalues[[x]][,'correlation']
            p.corr.coeff <- ggplot() +
                geom_boxplot(aes(y = corr.coeff, x = 1), outlier.color = 'gray') +
                ylab('Spearman correlation coefficient') +
                xlab(x) +
                geom_hline(yintercept = 0) +
                ggtitle(variable) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    plot.title = element_text(size = 16),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_text(size = 0),
                    axis.text.y = element_text(size = 12))
            if(isTRUE(plot.output) & length(assay.names) == 1) print(p.corr.coeff)
            return(p.corr.coeff)
        })
    names(all.corr.coeff.plots) <- levels(assay.names)

    ## put all boxplots  together ####
    everything <- datasets <- corr.coff <- NULL
    if (length(assay.names) > 1){
        printColoredMessage(
            message = '-- Putting all the boxplots of correlation coefficients together.',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.coeff <- sapply(levels(assay.names), function(x) all.corr.coeff.pvalues[[x]][,2]) %>%
            data.frame(.) %>%
            tidyr::pivot_longer(
                dplyr::everything(),
                names_to = 'datasets',
                values_to = 'corr.coff'
                )
        all.corr.coeff$datasets <- factor(
            x = all.corr.coeff$datasets,
            levels = levels(assay.names)
            )
        overall.corr.coeff.plot <- ggplot(all.corr.coeff, aes(x = datasets, y = corr.coff)) +
            geom_boxplot(outlier.color = 'gray') +
            ylab('Spearman correlation coefficients') +
            xlab('Datasets') +
            geom_hline(yintercept = 0) +
            ggtitle(variable) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                plot.title = element_text(size = 16),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12))
        overall.corr.coeff.plot <- annotate_figure(
            p = overall.corr.coeff.plot,
            top = text_grob(
                label = "Correlation analysis",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ',
                    'Correlation analysis between individaul gen expression an the  ',
                    variable,
                    ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
        )
        printColoredMessage(
            message = '- The boxplots of correlation coefficients of each assays are combined into one plot.',
            color = 'blue',
            verbose = verbose
            )
        if (isTRUE(plot.output)) suppressMessages(print(overall.corr.coeff.plot))
    }

    # Generate p-value histograms of the correlation coefficients ####
    printColoredMessage(
        message = '-- Generating p-value histograms of the correlations:',
        color = 'magenta',
        verbose = verbose
    )
    ## specified ylim ####
    breaks <- seq(from = 0, to = 1, by = .1)
    ylim.pvalue <- sapply(
        levels(assay.names),
        function(x) {
            binned <- cut(
                x = all.corr.coeff.pvalues[[x]][, 1],
                breaks = breaks,
                include.lowest = TRUE)
            frequency <- table(binned)[1]
        })
    ylim.pvalue <- ceiling(x = max(ylim.pvalue))
    ## generate histograms ####
    all.pvalue.hist.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Generating p-value histograms for the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            pvalues <- all.corr.coeff.pvalues[[x]]
            p.pvalues <- ggplot(pvalues, aes(x = pvalues[,1])) +
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
            if (isTRUE(plot.output) & length(assay.names) == 1) print(p.pvalues)
            return(p.pvalues)
        })
    names(all.pvalue.hist.plots) <- levels(assay.names)

    ## put all histograms together ####
    everything <- datasets <- corr.coff <- NULL
    if (length(assay.names) > 1){
        printColoredMessage(
            message = '-- Putting all the p-values histograms of the correlations together.',
            color = 'magenta',
            verbose = verbose
        )
        overall.pvalue.hist.plot <- ggarrange(
            plotlist = all.pvalue.hist.plots,
            ncol = plot.ncol,
            nrow = plot.nrow)
        if (class(overall.pvalue.hist.plot)[[1]] == 'list'){
            plot.list <- lapply(
                seq(length(overall.pvalue.hist.plot)),
                function(x){
                    annotate_figure(
                        p = overall.pvalue.hist.plot[[x]],
                        top = text_grob(
                            label = "Correlation analysis",
                            color = "orange",
                            face = "bold",
                            size = 18),
                        bottom = text_grob(
                            label = paste0(
                                'Analysis: ',
                                'Correlation analysis between individaul gen expression an the  ',
                                variable,
                                ' variable.'),
                            color = "black",
                            hjust = 1,
                            x = 1,
                            size = 10))
                })
            overall.pvalue.hist.plot <- ggarrange(
                plotlist = plot.list,
                ncol = 1,
                nrow = 1)
        } else {
            overall.pvalue.hist.plot <- annotate_figure(
                p = overall.pvalue.hist.plot,
                top = text_grob(
                    label = "Correlation analysis",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ', 'Correlation analysis between individaul gen expression an the  ', variable, ' variable.'),
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
        if(isTRUE(plot.output)) suppressMessages(print(overall.pvalue.hist.plot))
    }

    # Save the results ####
    printColoredMessage(
        message = '-- Save the all the boxplots of the correlation coefficients:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        ### save plot per assay ####
        printColoredMessage(
            message = paste0(
                '- Saving all the boxplots of the correlation coefficients and p-values histograms',
                ' of each assay(s) in the "metadata" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'Correlation',
            method = correlation.method,
            variables = variable,
            file.name = 'cor.coef.boxplot',
            results.data = all.corr.coeff.plots
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'Correlation',
            method = correlation.method,
            variables = variable,
            file.name = 'cor.pvalues.histogram',
            results.data = all.pvalue.hist.plots
            )
        printColoredMessage(
            message = paste0(
                '- The boxplot of the correlation coefficients and p-values histograms for the indiviaul assay(s) is  saved to',
                '  "se.obj@metadata$metric$AssayName$Correlation$',
                correlation.method,
                ' in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )

        ## save combined plots ####
        if (length(assay.names) > 1) {
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'Correlation',
                method = correlation.method,
                variables = variable,
                file.name = 'cor.coef.boxplot',
                plot.data = overall.corr.coeff.plot
                )
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'Correlation',
                method = correlation.method,
                variables = variable,
                file.name = 'cor.pvalues.histogram',
                plot.data = overall.pvalue.hist.plot)
            printColoredMessage(
                message = paste0(
                    '- The combined boxplots of the correlation coefficients of all the  indiviaul assay(s) is saved to',
                    '  "se.obj@metadata$plot$Correlation$',
                    correlation.method,
                    ' in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)

        }
        printColoredMessage(message = '------------The plotGenesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

    }
    ## out put the results as list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- All the plots are saved as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The genesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        if (length(assay.names) == 1){
            return(all.corr.coeff.plots = all.corr.coeff.plots)
        } else{
            return(gene.var.corr.plot = list(
                all.corr.coeff.plots = all.corr.coeff.plots,
                overall.corr.coeff.plot = overall.corr.coeff.plot
                ))
        }
    }
}
