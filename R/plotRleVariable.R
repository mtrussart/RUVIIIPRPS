#' Plots a variable against the medians and IQRs of relative log expression (RLE) data.

#' @author Ramyar Molania

#' @description
#' This function plots a variable against the medians and IQR of a relative log expression (RLE) data. Because of the
#' sensitivity of the RLE medians and IQR to unwanted variation, we  examine the relationships between RLE medians and
#' IQR with potential sources of unwanted variation. In the absence of any influence of unwanted variation in the data,
#' we should see no such associations.

#' @details
#' If the variable is categorical, boxplots of the RLE medians and IQR against the variable will be generated. If
#' the variable is continuous, scatter plots will be created.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or a vector of character strings representing the names of the assay(s)
#' in the SummarizedExperiment object to plot. The default is "all", which indicates that all assays in the SummarizedExperiment
#' object will be selected.
#' @param variable Character. The name of the column(s) in the sample annotation of the SummarizedExperiment object.
#' The variable can be either categorical or continuous. If a categorical variable is provided, boxplots of the RLE medians
#' and IQRs will be generated across the variable separately. If a continuous variable is provided, scatter plots of the RLE
#' medians and IQRs will be produced against the variable separately.
#' @param rle.data.type Character. Indicates which RLE data should be used for plotting. The options are 'rle.medians',
#' 'rle.iqr', or 'both'. If 'rle.medians' is selected, the RLE medians will be plotted against the variable. If 'rle.iqr'
#' is selected, the RLE IQRs will be plotted against the variable, and if 'both', both RLE medians and IQRs will be plotted
#' against the variable. The default is 'both'.
#' @param ylim.rle.med.plot Numeric. Specifies the ylim of the boxplot or scatter plots when the RLE medians are used. If
#' NULL, the function will automatically determine an appropriate ylim for the plots.
#' @param ylim.rle.iqr.plot Numeric. Specifies the ylim of the boxplot or scatter plots when the RLE IQRs are used. If
#' NULL, the function will automatically determine an appropriate ylim for the plots. The default is NULL.
#' @param points.size Numeric. Specifies the point size of the scatter plots. The default is 1.
#' @param plot.ncol Numeric. Specifies the number of columns in the plot grid.
#' @param plot.nrow Numeric. Specifies the number of rows in the plot grid.
#' @param plot.output Logical. If TRUE, individual RLE plot(s) will be printed while the function is running.
#' @param save.se.obj Logical. Specifies whether to save the plots in the metadata of the SummarizedExperiment object or
#' to output them as a list. The default is TRUE.
#' @param verbose Logical. If TRUE, messages for different steps of the function will be shown.
#'
#' @return A SummarizedExperiment object that contains all the plot(s) in the metadata or a list containing all the plot(s).


#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays
#' @importFrom ggpubr ggarrange stat_cor stat_compare_means
#' @import ggplot2
#' @export

plotRleVariable <- function(
        se.obj,
        assay.names = "all",
        variable,
        rle.data.type = 'both',
        ylim.rle.med.plot = NULL,
        ylim.rle.iqr.plot = NULL,
        points.size = 1,
        plot.ncol = 3,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotRleVariable function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of assay names(s) or "assay.names = all".')
    }
    if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }
    if(length(rle.data.type) > 1){
        stop('The "rle.data.type" must be one of the "rle.medians", "rle.iqrs" or "both".')
    }
    if (!rle.data.type %in% c('rle.medians', 'rle.iqrs', 'both')) {
        stop('The "rle.data.type" must be one of the "rle.medians", "rle.iqrs" or "both".')
    }
    if (!is.null(ylim.rle.med.plot)) {
        if (length(ylim.rle.med.plot) != 2)
            stop('Please specify the "ylim.rle.med.plot" argument.')
    }
    if (!is.null(ylim.rle.iqr.plot)) {
        if (length(ylim.rle.iqr.plot) != 2)
            stop('Please specify the "ylim.rle.iqr.plot" argument.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" variable contains NA. ',
                    'Run the checkSeObj function with "remove.na = both"',
                    ', then re-run the computeRLE function.'))
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain the rle data ####
    printColoredMessage(
        message = paste0('-- Obtaining the computed RLE medians and IQR from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    if (rle.data.type %in% c('rle.medians', 'both')){
        all.rle.medians <- getMetricFromSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment = 'RLE',
            assessment.type = 'global.level',
            method = 'gene.median.center',
            file.name = 'data',
            sub.file.name = 'rle.med',
            variables = 'general',
            required.function = 'computeRLE',
            message.to.print = 'medians of RLE'
        )
    }
    if (rle.data.type %in% c('rle.iqrs', 'both')){
        all.rle.iqrs <- getMetricFromSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment = 'RLE',
            assessment.type = 'global.level',
            method = 'gene.median.center',
            file.name = 'data',
            sub.file.name = 'rle.iqr',
            variables = 'general',
            required.function = 'computeRLE',
            message.to.print = 'IQR of RLE'
        )
    }
    # Generate the plots ####
    ## RLE medians ####
    if (rle.data.type %in% c('rle.medians', 'both')){
        ### generate plots between the RLE medians and the variable ####
        printColoredMessage(
            message = paste0('-- Generating plots between the RLE medians and the "variable":'),
            color = 'magenta',
            verbose = verbose
            )
        ## specify ylim ####
        if(is.null(ylim.rle.med.plot)){
            ylim.rle.med <- unlist(lapply(
                levels(assay.names),
                function(x){
                    c(min(all.rle.medians[[x]]),
                      max(all.rle.medians[[x]]))
                }))
            ylim.rle.med.plot <- c(min(ylim.rle.med), max(ylim.rle.med))
        }
        samples <- rle <- everything <- sample <- rle.medians <- r.label <- NULL
        all.rle.med.var.plots <- lapply(
            levels(assay.names),
            function(x) {
                rle.med.data <- data.frame(
                    rle.medians = all.rle.medians[[x]],
                    var = colData(se.obj)[[variable]])
                ## scatter plot ####
                if (class(colData(se.obj)[[variable]]) %in% c('numeric', 'integr')){
                    printColoredMessage(
                        message = paste0(
                            '- Generate a scatter plot between the RLE medians of the "',
                            x,
                            '" data and the "',
                            variable,
                            '" variable.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.medians)) +
                        geom_point(size = points.size, color = 'gray40', stroke = .2, pch = 21, alpha = .6) +
                        ggtitle(x) +
                        xlab(variable) +
                        ylab('RLE medians') +
                        geom_smooth(formula = y ~ x, method = 'lm', colour = "darkgreen") +
                        ggpubr::stat_cor(
                            aes(label = r.label),
                            color = "navy",
                            label.y = max(ylim.rle.med)) +
                        coord_cartesian(ylim = ylim.rle.med.plot) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', linewidth = 1),
                              axis.title.x = element_text(size = 12),
                              axis.title.y = element_text(size = 12),
                              axis.text.x = element_text(size = 9),
                              axis.text.y = element_text(size = 9),
                              legend.position = 'bottom')
                }
                ## boxplot plot ####
                if (class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                    printColoredMessage(
                        message = paste0(
                            '- Generating a boxplot between the RLE medians of the "',
                            x,
                            '" data and the "',
                            variable,
                            '" variable.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.medians)) +
                        geom_boxplot(outlier.color = 'gray') +
                        stat_compare_means(method = 'anova', label.x  = 1, label.y = ylim.rle.med.plot[2]-.1, color = 'navy') +
                        ggtitle(x) +
                        xlab(variable) +
                        ylab('RLE medians') +
                        coord_cartesian(ylim = ylim.rle.med.plot) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', linewidth = 1),
                              axis.title.x = element_text(size = 12),
                              axis.title.y = element_text(size = 12),
                              axis.text.x = element_text(size = 8, angle = 35, vjust = 1, hjust = 1),
                              axis.text.y = element_text(size = 9),
                              legend.position = 'bottom')
                }
                if (isTRUE(plot.output) & length(assay.names) == 1) print(p.rle)
                p.rle
            })
        names(all.rle.med.var.plots) <- levels(assay.names)
        ### put all the RLE plots of all assays ####
        if(length(assay.names) > 1){
            printColoredMessage(
                message = '-- Putting all the plots together:',
                color = 'magenta',
                verbose = verbose
                )
            overall.rle.med.var.plots <- ggpubr::ggarrange(
                plotlist = all.rle.med.var.plots,
                ncol = plot.ncol,
                nrow = plot.nrow,
                common.legend = TRUE
                )
            if (class(overall.rle.med.var.plots)[[1]] == 'list'){
                plot.list <- lapply(
                    seq(length(overall.rle.med.var.plots)),
                    function(x){
                        annotate_figure(
                            p = overall.rle.med.var.plots[[x]],
                            top = text_grob(
                                label = "The medians of RLE against variable",
                                color = "orange",
                                face = "bold",
                                size = 18),
                            bottom = text_grob(
                                label = paste0(
                                    'Analysis: ',
                                    'the medians of the RLE against the variable. \n',
                                    'Variable ',
                                    variable),
                                color = "black",
                                hjust = 1,
                                x = 1,
                                size = 10))
                    })
                overall.rle.med.var.plots <- ggarrange(
                    plotlist = plot.list,
                    ncol = 1,
                    nrow = 1)
            } else {
                overall.rle.med.var.plots <- annotate_figure(
                    p = overall.rle.med.var.plots,
                    top = text_grob(
                        label = "The medians of the RLE against variable",
                        color = "orange",
                        face = "bold",
                        size = 18),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ', 'the medians of the RLE against the variable. \n',
                            'Variable ', variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            }
            printColoredMessage(
                message = '- The plots of individual assays are combined into one.',
                color = 'blue',
                verbose = verbose)
            if (isTRUE(plot.output)) print(overall.rle.med.var.plots)
        }
    }

    ## RLE IQRs ####
    if(rle.data.type %in% c('rle.iqrs', 'both')){
        ### generate plots between the RLE medians and the variable ####
        printColoredMessage(
            message = paste0('-- Generating plots between the "variable" and the RLE IQRs:'),
            color = 'magenta',
            verbose = verbose
        )
        ## specify ylim ####
        if(is.null(ylim.rle.iqr.plot)){
            ylim.rle.iqr <- unlist(lapply(
                levels(assay.names),
                function(x){
                    c(min(all.rle.iqrs[[x]]),
                      max(all.rle.iqrs[[x]]))
                }))
            ylim.rle.iqr.plot <- c(min(ylim.rle.iqr), max(ylim.rle.iqr))
        }
        samples <- rle <- everything <- sample <- rle.iqr <- r.label <- NULL
        all.rle.iqr.var.plots <- lapply(
            levels(assay.names),
            function(x) {
                rle.iqr.data <- data.frame(
                    rle.iqr = all.rle.iqrs[[x]],
                    var = colData(se.obj)[[variable]])
                ## scatter plot ####
                if (class(colData(se.obj)[[variable]]) %in% c('numeric', 'integr')){
                    printColoredMessage(
                        message = paste0(
                            '- Generating a scatter plot between the RLE IQR of the ',
                            x,
                            ' data and the ',
                            variable,
                            ' variable.'),
                        color = 'blue',
                        verbose = verbose)
                    p.rle <- ggplot(rle.iqr.data, aes(x = var, y = rle.iqr)) +
                        geom_point(size = points.size, color = 'gray40', stroke = .2, pch = 21, alpha = .6) +
                        ggtitle(x) +
                        xlab(variable) +
                        ylab('RLE IQRs') +
                        geom_smooth(formula = y ~ x, method = 'lm', colour = "darkgreen") +
                        ggpubr::stat_cor(aes(label = r.label), color = "navy", label.y = max(ylim.rle.iqr)) +
                        coord_cartesian(ylim = ylim.rle.iqr.plot) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', linewidth = 1),
                              axis.title.x = element_text(size = 12),
                              axis.title.y = element_text(size = 12),
                              axis.text.x = element_text(size = 9),
                              axis.text.y = element_text(size = 9),
                              legend.position = 'bottom')
                }
                ## boxplot ####
                if (class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                    printColoredMessage(
                        message = paste0(
                            '- Generate a boxplot between the RLE IQR of the ', x,
                            ' data and the ', variable, ' variable.'),
                        color = 'blue',
                        verbose = verbose)
                    p.rle <- ggplot(rle.iqr.data, aes(x = var, y = rle.iqr)) +
                        geom_boxplot(outlier.color = 'gray') +
                        stat_compare_means(method = 'anova', label.x  = 1, label.y = ylim.rle.iqr.plot[2]-.1, color = 'navy') +
                        ggtitle(x) +
                        xlab(variable) +
                        ylab('RLE IQRs') +
                        coord_cartesian(ylim = ylim.rle.iqr.plot) +
                        theme(panel.background = element_blank(),
                              axis.line = element_line(colour = 'black', linewidth = 1),
                              axis.title.x = element_text(size = 12),
                              legend.position = 'bottom',
                              axis.title.y = element_text(size = 12),
                              axis.text.x = element_text(size = 8, angle = 35, vjust = 1, hjust = 1),
                              axis.text.y = element_text(size = 9))
                }
                if (isTRUE(plot.output) & length(assay.names) == 1) print(p.rle)
                p.rle
            })
        names(all.rle.iqr.var.plots) <- levels(assay.names)

        ### put all the RLE plots of all assays ####
        if (length(assay.names) > 1){
            printColoredMessage(
                message = '-- Putting all the plots together:',
                color = 'magenta',
                verbose = verbose
                )
            overall.rle.iqr.var.plots <- ggpubr::ggarrange(
                plotlist = all.rle.iqr.var.plots,
                ncol = plot.ncol,
                nrow = plot.nrow,
                common.legend = TRUE
                )
            if (class(overall.rle.iqr.var.plots)[[1]] == 'list'){
                plot.list <- lapply(
                    seq(length(overall.rle.iqr.var.plots)),
                    function(x){
                        annotate_figure(
                            p = overall.rle.iqr.var.plots[[x]],
                            top = text_grob(
                                label = "The IQR of the RLE against variable",
                                color = "orange",
                                face = "bold",
                                size = 18),
                            bottom = text_grob(
                                label = paste0(
                                    'Analysis: ',
                                    ' the IQR of the RLE data against the variable. \n',
                                    'Variable: ',
                                    variable),
                                color = "black",
                                hjust = 1,
                                x = 1,
                                size = 10))
                    })
                overall.rle.iqr.var.plots <- ggarrange(
                    plotlist = plot.list,
                    ncol = 1,
                    nrow = 1)
            } else {
                overall.rle.iqr.var.plots <- annotate_figure(
                    p = overall.rle.iqr.var.plots,
                    top = text_grob(
                        label = "The IQR of the RLE against variable",
                        color = "orange",
                        face = "bold",
                        size = 18),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ',
                            ' the IQR of the RLE data against the variable. \n',
                            'Variable: ',
                            variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            }
            printColoredMessage(
                message = '- The individual assay plots are combined into one.',
                color = 'blue',
                verbose = verbose)
            if (plot.output) print(overall.rle.iqr.var.plots)
        }
    }

    # Save the plots ####
    printColoredMessage(
        message = '-- Saving all the RLE plots:',
        color = 'magenta',
        verbose = verbose
        )
    ## add plots to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '-- Save all the RLE plots to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (rle.data.type %in% c('rle.medians', 'both')){
            if (class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                file.name <- 'boxplot'
            } else  file.name <- 'scatter.plot'
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                assay.names = assay.names,
                slot = 'Metrics',
                assessment.type = 'global.level',
                assessment = 'RLE',
                method = 'corr.medians.variable',
                file.name = file.name,
                variables = variable,
                results.data = all.rle.med.var.plots
            )
        }
        if (rle.data.type %in% c('rle.iqrs', 'both')){
            if(class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                file.name <- 'boxplot'
            } else  file.name <- 'scatter.plot'
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                assay.names = assay.names,
                slot = 'Metrics',
                assessment.type = 'global.level',
                assessment = 'RLE',
                method = 'corr.iqr.variable',
                file.name = file.name,
                variables = variable,
                results.data = all.rle.iqr.var.plots
            )
        }
        printColoredMessage(
            message = paste0(
                '- The plot of individual assay(s) are saved to the',
                ' "se.obj@metadata$metric$AssayName$RLE$rle.plot$rle.var.plot" in the SummarizedExperiment object..'),
            color = 'blue',
            verbose = verbose)
        }

        ## add overall RLE plots of all assays ####
        if (length(assay.names) > 1){
            if(rle.data.type %in% c('rle.medians', 'both')){
                if(class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                    file.name <- 'boxplot'
                } else  file.name <- 'scatter.plot'
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'RLE',
                    method = 'corr.medians.variable',
                    variables = variable,
                    file.name = file.name,
                    plot.data = overall.rle.med.var.plots
                )
            }
            if (rle.data.type %in% c('rle.iqrs', 'both')){
                if(class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
                    file.name <- 'boxplot'
                } else  file.name <- 'scatter.plot'
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'RLE',
                    method = 'corr.iqrs.variable',
                    variables = variable,
                    file.name = file.name,
                    plot.data = overall.rle.iqr.var.plots
                )
            printColoredMessage(
                message = paste0('- The combined plots of all assays are saved to the',
                                 ' "se.obj@metadata$plot$RLE$rle.var.plot" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(
            message = '------------The plotRleVariable function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj = se.obj)

        } else {
            printColoredMessage(message = '------------The plotRleVariable function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj = se.obj)
        }
    ## save the plots as list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '-- All the plots are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotRleVariable function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            if(rle.data.type == 'both'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots)
            } else if (rle.data.type == 'rle.medians'){
                rle.var.plots <- list(all.rle.med.var.plots = all.rle.med.var.plots)
            } else if (rle.data.type == 'rle.iqrs'){
                rle.var.plots <- list(all.rle.iqr.var.plots = all.rle.iqr.var.plots)
            }
        } else if (length(assay.names) > 1){
            if(rle.data.type == 'both'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots,
                    overall.rle.med.var.plots = overall.rle.med.var.plots,
                    overall.rle.iqr.var.plots = overall.rle.iqr.var.plots)
            } else if (rle.data.type == 'rle.medians'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    overall.rle.med.var.plots = overall.rle.med.var.plots)
            } else if (rle.data.type == 'rle.iqrs'){
                rle.var.plots <- list(
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots,
                    overall.rle.iqr.var.plots = overall.rle.iqr.var.plots)
            }
        }
    }
}
