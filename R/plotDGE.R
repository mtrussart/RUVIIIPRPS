#' Plot p-values histograms of DGE analysis.

#' @author Ramyar Molania

#' @description
#' This function plots the p-values histograms of the DGE analysis.

#' @details
#' DE analyses is performed using the Wilcoxon signed-rank test with log-transformed data e.g. raw counts, normalized data, ....
#' To evaluate the effects of the different sources of unwanted variation on the data, DE analyses is performed across
#' batches. In the absence of any batch effects, the histogram of the resulting un-adjusted P values should be uniformly
#' distributed

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute DGE. The default is set to 'all' indicating all the the assays of the
#' SummarizedExperiment object will be selected.
#' @param variable Symbol. Specifies the column name in the SummarizedExperiment object containing a categorical variable
#' , such as sample types or batches
#' @param plot.ncol Numeric. A numeric value indicating the number of columns in the plot grid. This setting applies when
#' more than one data is provided. The default is set to 1.
#' @param plot.nrow Numeric. A numeric value indicating the number of rows in the plot grid. This setting applies when
#' more than one data is provided. The default is set to 1.
#' @param plot.output Logical. Indicates whether to plot the p-values histograms when the functions is running, by default
#' it is set to 'FALSE'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export


plotDGE <- function(
        se.obj,
        assay.names = 'all',
        variable,
        plot.ncol = 1,
        plot.nrow = 2,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotDGE function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain computed p-values for each DE contrasts ####
    printColoredMessage(
        message = paste0('-- Obtain computed p-values of all contrasts of the DGE analysis for the "', variable, '" variable.') ,
        color = 'magenta',
        verbose = verbose
    )
    all.de.tests <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = assay.names,
        assessment = 'DGE',
        assessment.type = 'gene.level',
        method = 'Wilcoxon',
        variables = variable,
        file.name = 'p.values',
        sub.file.name = NULL,
        required.function = 'computeDGE',
        message.to.print = 'DGE'
    )
    p.vals <- everything <- NULL

    # Generate p-values histograms for each contrast ####
    ## specified ylim for histograms ####
    breaks <- seq(from = 0, to = 1, by = .1)
    ylim.pvalue <- sapply(
        levels(assay.names),
        function(x) {
            sapply(
                names(all.de.tests[[x]]),
                function(y){
                    binned <- cut(
                        x = all.de.tests[[x]][[y]][, 3],
                        breaks = breaks,
                        include.lowest = TRUE)
                    frequency <- table(binned)[1]
                })
        })
    ylim.pvalue <- ceiling(x = max(ylim.pvalue))

    ## generate p-values histograms for each assay  ####
    all.pval.histograms <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Generate p-values histograms for the ', x, ' data.'),
                color = 'blue',
                verbose = verbose
            )
            if(length(unique(colData(se.obj)[[variable]])) == 2 ){
                pval.data <- do.call(rbind, all.de.tests[[x]])
                pval.data$contrasts <- rep(names(all.de.tests[[x]]), each = nrow(se.obj))
                pval.plot <- ggplot(pval.data, aes(x = pvalue)) +
                    geom_histogram(binwidth = 0.1) +
                    scale_y_continuous(labels = function(x) format(x / 1000, scientific = F), limits = c(0, ylim.pvalue)) +
                    ggtitle(x) +
                    xlab('p-values') +
                    ylab(expression('Frequency'~10^3)) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8)
                        )
            } else if(length(unique(colData(se.obj)[[variable]])) > 2 ) {
                pval.data <- do.call(rbind, all.de.tests[[x]])
                pval.data$contrasts <- rep(names(all.de.tests[[x]]), each = nrow(se.obj))
                pval.plot <- ggplot(pval.data, aes(x = pvalue)) +
                    geom_histogram(binwidth = 0.1) +
                    scale_y_continuous(labels = function(x) format(x /1000, scientific = F), limits = c(0, ylim.pvalue)) +
                    ggtitle(x) +
                    xlab('p-values') +
                    ylab(expression('Frequency'~10^3)) +
                    facet_wrap(~contrasts) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        strip.text = element_text(size = 8),
                        strip.text.x.top = element_text(size = 8),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8))
            }
            if(isTRUE(plot.output) & length(assay.names) == 1) print(pval.plot)
            return(pval.plot)
        })
    names(all.pval.histograms) <- levels(assay.names)

    ## put all the p-values histograms into one  ####
    printColoredMessage(
        message = '- Put all the p-values histograms together:' ,
        color = 'magenta',
        verbose = verbose
        )
    overall.pval.histograms <- ggpubr::ggarrange(
        plotlist = all.pval.histograms,
        ncol = plot.ncol,
        nrow = plot.nrow
        )
    if(class(overall.pval.histograms)[[1]] == 'list'){
        plot.list <- lapply(
            seq(length(overall.pval.histograms)),
            function(x){
                annotate_figure(
                    p = overall.pval.histograms[[x]],
                    top = text_grob(
                        label = "Differential gene expression analysis.",
                        color = "orange",
                        face = "bold",
                        size = 18),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ', 'differential gene expression analysis using Wilcoxon for all possible contrasts\n',
                            "Variable: ", variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            })
        overall.pval.histograms <- ggarrange(
            plotlist = plot.list,
            ncol = 1,
            nrow = 1
            )
    } else {
        overall.pval.histograms <- annotate_figure(
            p = overall.pval.histograms,
            top = text_grob(
                label = "Differential gene expression analysis.",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ', 'differential gene expression analysis using Wilcoxon for all possible contrasts\n',
                    "Variable: ", variable),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10))
    }
    if(isTRUE(plot.output)) suppressMessages(print(overall.pval.histograms))

    # Save the plots ####
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'DGE',
            method = 'Wilcoxon',
            variables = variable,
            file.name = 'plot',
            results.data = all.pval.histograms
            )
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) > 1){
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'DGE',
                method = 'Wilcoxon',
                variables = variable,
                file.name = 'histogram',
                plot.data = overall.pval.histograms
            )
            printColoredMessage(
                message = paste0('The p-value histograms of all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(
            message = '------------The plotDGE function finished.',
            color = 'white',
            verbose = verbose
        )
        return(se.obj = se.obj)
    }
    ## return the results as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = 'All the p-value histograms are outputed as list.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotDGE function finished.',
            color = 'white',
            verbose = verbose
            )
        if(length(assay.names) == 1){
            return(all.pval.histograms = all.pval.histograms)
        } else {
            return(list(
                all.pval.histograms = all.pval.histograms,
                overall.pval.histograms = overall.pval.histograms)
            )
        }
    }
}




