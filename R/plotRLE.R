#' Generates boxplot of the relative log expression (RLE) distributions of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function generates a boxplot for individual RLE data. The function does not plot the outliers in the boxplot in
#' order to enhance the visibility of the variation in the RLE medians and interquartile ranges. The RLE data can be
#' obtained using the 'computeRLE' function. We refer to the 'computeRLE' function for more detail about RLE.

#' @details
#' An ideal RLE plot should have its medians centered around zero, and its box widths and their interquartile ranges
#' (IQRs) should be similar in magnitude. We refer to Gandolfo L. C. & Speed, T. P., PLoS ONE, 2018 for more details and
#' assumption about the RLE plots.

#' @references
#' Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector of character strings for the selection of the name(s) of
#' the assay(s) in the SummarizedExperiment object to plot the RLE data. The default is set to "all", which indicates that
#' all the assays of the SummarizedExperiment object will be selected.
#' @param variable Character. Indicates the name of the column in the sample annotation of the SummarizedExperiment object.
#' The interquartile ranges of the RLE boxplots will be colored based on the specified variable. The variable must be a
#' categorical variable. The default is set to 'NULL'.
#' @param variable.colors Character. A vector of colors for the interquartile ranges of the RLE boxplots if the 'variable'
#' is specified. The default is set to 'NULL'. This means the function will select required colors.
#' @param ylim.rle.plot Numeric. A vector of two values to specify the ylim of the RLE plot(s). If 'NULL', the function
#' uses the minimum and maximum interquartile ranges of all the RLE data to specify the ylim. The default is 'NULL'. The
#' ylim of the RLE plots should be the same to be able to compare them against each other.
#' @param iqr.width Numeric. A numeric value indicating the width size of RLE interquartile ranges in the plot. The default
#' is set to 1.
#' @param median.points.size Numeric. A numeric value indicating the size of the points of the RLE medians in the boxplot(s).
#' The default is set to 1.
#' @param median.points.color Character. Specifies the color of the points representing RLE medians in the boxplot(s).
#' @param geom.hline.color Character. Indicates the color of the horizontal line (geom.hline) crossing 0 in the RLE boxplot(s).
#' This line aids in visualizing the deviation of the RLE medians within the RLE boxplot(s).
#' @param plot.ncol Numeric. A numeric value indicating the number of columns in the plot grid. When the number of selected
#' assays is more than 1, the function puts all the RLE boxplots in one grid. The default is set to 2.
#' @param plot.nrow Numeric. A numeric value indicating the number of rows in the plot grid. When the number of selected
#' assays is more than 1, the function puts all the RLE boxplots in one grid. The default is set to 3.
#' @param plot.output Logical. If 'TRUE', the individual RLE plot(s) will be printed while the function is running.
#' @param save.se.obj Logical. Indicates whether to save the RLE plots in the metadata of the SummarizedExperiment object
#'  or to output the result as a list. By default, it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object that contains all the RLE plot(s) in the "metadata" or a list that contains all
#' the RLE plot(s).


#' @importFrom matrixStats colQuantiles colMedians colIQRs
#' @importFrom SummarizedExperiment assays
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @importFrom scales hue_pal
#' @import ggplot2
#' @export

plotRLE <- function(
        se.obj,
        assay.names = "all",
        variable = NULL,
        variable.colors = NULL,
        ylim.rle.plot = NULL,
        iqr.width = 1,
        median.points.size = 1,
        median.points.color = 'black',
        geom.hline.color = 'cyan',
        plot.ncol = 2,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotRLE function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty or "NULL".')
    } else if(is.logical(assay.names) | is.list(assay.names)){
        stop('The "assay.names" must be either a vector of assay name(s) or set to "all".')
    }
    if (!is.null(variable)){
        if(is.logical(variable)){
            stop('The "variable" must contain only one variable name.')
        } else if (length(variable) > 1){
            stop('The "variable" must contain only one variable name.')
        } else if (!class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
            stop('The "variable" must be a categorical variable.')
        } else if (sum(is.na(se.obj@colData[[variable]])) > 0){
            stop(paste0('The "', variable, '" variable contains NA. ',
                        'Run the checkSeObj function with "remove.na = both"',
                        ', then re-run the "computeRLE" function.'))
        }
    }
    if (!is.null(ylim.rle.plot) & !length(ylim.rle.plot) != 2) {
        stop('Please specify the approprate "ylim.rle.plot" argument.')
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical(TRUE or FALSE).')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical(TRUE or FALSE).')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical(TRUE or FALSE).')
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Select colors ####
    if (!is.null(variable)){
        printColoredMessage(
            message = paste0('-- Checking the colors for the boxplots:'),
            color = 'magenta',
            verbose = verbose
        )
        if(is.null(variable.colors)){
            printColoredMessage(
                message = paste0('- Selecting required colors for the boxplots:'),
                color = 'blue',
                verbose = verbose
            )
            n.colors <- length(unique(colData(se.obj)[[variable]]))
            rle.plot.colors <- scales::hue_pal()(n.colors*2)
            rle.plot.colors <- rle.plot.colors[seq(1, n.colors*2, 2)]
        } else if (!is.null(variable.colors)){
            printColoredMessage(
                message = paste0('- Checking the provided colors for the boxplots:'),
                color = 'blue',
                verbose = verbose
            )
            if (length(variable.colors) != length(unique(colData(se.obj)[[variable]]))){
                stop('The provided colors do not match the levels of the provided variable.')
            } else rle.plot.colors <- variable.colors
        }
    }

    # Obtain rle data ####
    printColoredMessage(
        message = paste0('-- Obtaining the computed RLE data from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
        )
    all.rle.data <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = levels(assay.names),
        assessment = 'RLE',
        assessment.type = 'global.level',
        method = 'gene.median.center',
        file.name = 'data',
        sub.file.name = 'rle.data',
        variables = 'general',
        required.function = 'computeRLE',
        message.to.print = 'RLE data'
        )
    # Generate the RLE plots ####
    printColoredMessage(
        message = '-- Generating the boxplots of the RLE data  :',
        color = 'magenta',
        verbose = verbose
        )
    # Specify ylim for the RLE plots ####
    if (is.null(ylim.rle.plot)){
        printColoredMessage(
            message = '- ylim is not provided, then specifying the same ylim for all the RLE plots:',
            color = 'blue',
            verbose = verbose
            )
        ylim.rle.plot <- abs(unlist(lapply(
            levels(assay.names),
            function(x){
                samples.quantiles <- matrixStats::colQuantiles(
                    x = all.rle.data[[x]],
                    probs = c(0.2, 0.8))
                c(max(samples.quantiles), min(samples.quantiles))
            })))
        ylim.rle.plot <- c(-max(ylim.rle.plot)[1], max(ylim.rle.plot)[1])
    }
    rle <- sample <- NULL
    all.rle.plots <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('- Generating the RLE plot for the "', x , '" data.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- all.rle.data[[x]]
            samples.quantiles <- matrixStats::colQuantiles(
                x = rle.data,
                probs = seq(from = 0, to = 1, by = 0.25))
            samples.quantiles <- as.data.frame(samples.quantiles[, c(2:4)])
            colnames(samples.quantiles) <- c('per25', 'medians', 'per75')
            samples.quantiles$sample <- paste0('Sam', 1:ncol(rle.data))
            if(!is.null(variable)){
                # colored RLE plots ####
                samples.quantiles$variable <- colData(se.obj)[[variable]]
                samples.quantiles <- tidyr::pivot_longer(
                    data = samples.quantiles,
                    -c(sample, variable),
                    values_to = 'rle',
                    names_to = 'range')
                samples.quantiles$sample <- factor(samples.quantiles$sample, levels = paste0('Sam', 1:ncol(rle.data)))
                p.rle <- ggplot(samples.quantiles, aes(x = sample, y = rle, group = sample, color = variable)) +
                    geom_line(linewidth = iqr.width) +
                    geom_point(data = samples.quantiles[samples.quantiles$range == 'medians',],
                               aes(group = range),
                               size = median.points.size ,
                               colour = median.points.color) +
                    scale_color_manual(name = 'Groups:', values = rle.plot.colors) +
                    ylab('RLE') +
                    xlab('Samples') +
                    ggtitle(x) +
                    coord_cartesian(ylim = ylim.rle.plot) +
                    geom_hline(yintercept = 0, colour = geom.hline.color) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 9),
                        axis.ticks.x = element_blank(),
                        legend.position = 'bottom') +
                guides(color = guide_legend(override.aes = list(linewidth = iqr.width *3)))
            } else if (is.null(variable)){
                # general RLE plots ####
                samples.quantiles <- tidyr::pivot_longer(
                    data = samples.quantiles,
                    -sample,
                    values_to = 'rle',
                    names_to = 'range')
                samples.quantiles$sample <- factor(samples.quantiles$sample, levels = paste0('Sam', 1:ncol(rle.data)))
                p.rle <- ggplot(samples.quantiles, aes(x = sample, y = rle, group = sample)) +
                    geom_line(linewidth = iqr.width) +
                    geom_point(data = samples.quantiles[samples.quantiles$range == 'medians',],
                               aes(group = range),
                               size = median.points.size ,
                               colour = median.points.color) +
                    ggtitle(x) +
                    xlab('Samples') +
                    ylab('RLE') +
                    labs(caption = paste0(
                        'Analysis: ', 'centred gene medians and then boxplots of samples.\n',
                        'Highlighted ', median.points.color, ' dots are boxplots medians.')) +
                    coord_cartesian(ylim = ylim.rle.plot) +
                    geom_hline(yintercept = 0, colour = geom.hline.color) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        plot.caption = element_text(hjust = 0, vjust = 0),
                        plot.title = element_text(size = 12),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 9),
                        axis.ticks.x = element_blank())
            }
            if (isTRUE(plot.output) & length(assay.names) == 1) print(p.rle)
            p.rle
        })
    names(all.rle.plots) <- levels(assay.names)

    # Generate the overall RLE plots ####
    if(length(assay.names) > 1){
        for(i in levels(assay.names) ){
            all.rle.plots[[i]]$labels$caption <- NULL
        }
        printColoredMessage(
            message = '-- Put all the RLE plots of all assays together.',
            color = 'magenta',
            verbose = verbose
            )
        overall.rle.plot <- ggpubr::ggarrange(
            plotlist = all.rle.plots,
            ncol = plot.ncol,
            nrow = plot.nrow,
            common.legend = TRUE
            )
        if(class(overall.rle.plot)[[1]] == 'list'){
            plot.list <- lapply(
                seq(length(overall.rle.plot)),
                function(x){
                    annotate_figure(
                        p = overall.rle.plot[[x]],
                        top = text_grob(
                            label = "Relative log expression plots",
                            color = "orange",
                            face = "bold",
                            size = 18),
                        bottom = text_grob(
                            label = paste0(
                                'Analysis: ', 'centred gene medians and then boxplots of samples.\n',
                                'Highlighted ', median.points.color, ' dots are boxplots medians.'),
                            color = "black",
                            hjust = 1,
                            x = 1,
                            size = 10))
                })
            overall.rle.plot <- ggarrange(
                plotlist = plot.list,
                ncol = 1,
                nrow = 1)
        } else {
            overall.rle.plot <- annotate_figure(
                p = overall.rle.plot,
                top = text_grob(
                    label = "Relative log expression plots",
                    color = "black",
                    face = "bold",
                    size = 18
                    ),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'centred gene medians and then boxplots of samples.\n',
                        'Highlighted ',
                        median.points.color,
                        ' dots are boxplots medians.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10))
        }
        printColoredMessage(
            message = '- All the RLE plots of all assays together are combined into one plot.',
            color = 'blue',
            verbose = verbose
        )
        if(isTRUE(plot.output)) print(overall.rle.plot)
    }

    # Save the plots ####
    printColoredMessage(
        message = '-- Save all the RLE plots:',
        color = 'magenta',
        verbose = verbose
        )
    ## add all the plots to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the RLE plots to the "metadata" of the SummarizedExperiment object.',
            color = 'black',
            verbose = verbose
            )
        ## add RLE plots of individual assays ####
        if(is.null(variable)){
            file.name <- 'un.colored'
            variables <- 'general'
        } else {
            file.name <- 'colored'
            variables <- variable
        }
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            slot = 'Metrics',
            assessment = 'RLE',
            assessment.type = 'global.level',
            method = 'gene.median.center',
            file.name = file.name,
            variables = variables,
            results.data = all.rle.plots
        )
        printColoredMessage(
            message = paste0(
                '* The RLE plot of individual assay(s) are saved to the ',
                '"se.obj@metadata$metric$AssayName$RLE$rle.plot" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )

        ## add overall RLE plots of all assays ####
        if(length(assay.names) > 1){
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'global.level',
                assessment = 'RLE',
                method = 'gene.median.center',
                variables = variables,
                file.name = file.name,
                plot.data = overall.rle.plot
                )
            printColoredMessage(
                message = paste0(
                    '* The combined RLE plots of all assays are saved to the',
                    ' "se.obj@metadata$plot$RLE" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
                )
        }
        printColoredMessage(
            message = '------------The plotRLE function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj = se.obj)
    }

    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '-- All the RLE plots are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(rle.plots = list(all.rle.plots = all.rle.plots))
        } else{
            return(rle.plots = list(
                all.rle.plots = all.rle.plots,
                overall.rle.plot = overall.rle.plot))
        }
    }
}
