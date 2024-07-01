#' Generates scatter and boxplot plot of principal components.

#' @description
#' This functions generates scatter and boxplot the first principal components of the assay(s) of a SummarizedExperiment
#' object. The function can generate pairwise scatter plots of the first principal components colored by a categorical
#' variable or creates boxplots of each PC across the variable. Boxplots of principal components become beneficial when
#' the levels of a categorical variable are numerous, making visualization through colored scatter plots challenging. If
#' a continuous variable is provided, the function creates a scatter plot of each PC against the variable.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to PC. The default is "all, which indicates all the assays of the SummarizedExperiment
#' object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable can be either a categorical or continuous. If a continuous variable is provided, the function creates
#' a scatter plot of each PC against the variable. If a categorical variable is given, the function can generate pairwise
#' scatter plots of the first principal components colored by a categorical variable or creates boxplots of each PC across
#' the variable.
#' @param fast.pca Logical. Indicates whether to use the computed fast PC ro not. The default is 'TRUE'. We refe to the
#' computePCA function for more details.
#' @param nb.pcs Numeric. A numeric value indicating the number of first PCs to be used for plotting. The default is set to 3.It's important to note
#' that the variation of PCs for a portion of all selected PCs will be based solely on those selected PCs.
#' @param plot.type Symbol. A symbol specifying which plot type should be generated. Options are: 'scatter' and 'boxplot'.
#' The default is set to 'scatter'. Please note, the 'plot.type' cannot be set to 'scatter' of a categorical variable is
#' provided.
#' @param variable.colors Symbol. A vector of color names in scatter plots of PCs when a categorical variable is specified.
#' If is 'NULL', the function will use the default colors.
#' @param points.size Numeric. A numeric value specifying the size of points in the scatter PCA plots. The default is set
#' to 1.
#' @param stroke.color Symbol. A symbol specifying the color of the stroke of the points in the scatter PCA plots. The
#' default is set to 'gray'.
#' @param stroke.size Numeric. A numeric value indicating the size of the stroke of the points in the scatter PCA plots.
#' The default is set to 0.1.
#' @param points.alpha Numeric. A numeric value indicating the transparency of the points in the scatter PCA plots. The
#' default is set to 0.5.
#' @param densities.alpha Numeric.  A numeric value indicating the transparency of the densities in the scatter PCA plots.
#' The default is set to 0.5.
#' @param legend.position Symbol. A symbol indicating the location of the plot legend.  This setting will apply when a categorical
#' variable is provided and the 'plot.type' is set to 'scatter'. Options are: 'top', 'right', 'bottom' and 'left'. The default
#' is set to "bottom".
#' @param plot.ncol Numeric. A numeric value indicating the number of columns in the plot grid. When the number of selected
#' assay is more than 1, the function puts all the PCA plots in one grid.
#' @param plot.nrow Numeric. A numeric value indicating number of rows in the plot grid. When the number of selected assay is more than
#' 3, the function puts all the PCA plots in one grid.
#' @param plot.output Logical. If 'TRUE', the individual PCA plot(s) will be printed while functions is running.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object
#' or to output the results as list. By default it is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object that contains all the PCA plot(s) in the metadata or a list that contains all
#' the PCA plot(s).

#' @importFrom ggpubr ggarrange theme_pubr stat_cor
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom tidyr pivot_longer
#' @importFrom utils combn
#' @import ggplot2
#' @export

plotPCA <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        plot.type = 'scatter',
        variable.colors = NULL,
        points.size = 1,
        stroke.color = 'gray',
        stroke.size = 0.1,
        points.alpha = 0.5,
        densities.alpha = 0.5,
        legend.position = 'bottom',
        plot.ncol = 3,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotPCA function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)){
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(colData(se.obj))) {
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (is.null(nb.pcs)){
        stop('The "nb.pcs" cannot be empty.')
    } else if (!plot.type %in% c('scatter', 'boxplot')){
        stop('The "plot.type" must be one of "scatter" or "boxplot".')
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Select colors ####
    if(!is.null(variable)){
        if(is.null(variable.colors)){
            n.colors <- length(unique(colData(se.obj)[[variable]]))
            pca.plot.colors <- scales::hue_pal()(n.colors*2)
            pca.plot.colors <- pca.plot.colors[seq(1, n.colors*2, 2)]
        } else if (!is.null(variable.colors)){
            pca.plot.colors <- variable.colors
        }
    }

    # Obtain computed PCs from the SummarizedExperiment object ####
    printColoredMessage(
        message = paste0('-- Obtain the first ', nb.pcs, ' computed PCs from the SummarizedExperiment object.'),
        color = 'magenta',
        verbose = verbose
        )
    if(isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    all.pca.data <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = assay.names,
        slot = 'Metrics',
        assessment = 'PCA',
        assessment.type = 'global.level',
        method = method,
        file.name = 'data',
        variables = 'general',
        sub.file.name = 'svd',
        required.function = 'computePCA',
        message.to.print = 'PCs'
        )
    all.pca.percentage <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = assay.names,
        slot = 'Metrics',
        assessment = 'PCA',
        assessment.type = 'global.level',
        method = method,
        file.name = 'data',
        variables = 'general',
        sub.file.name = 'percentage.variation',
        required.function = 'computePCA',
        message.to.print = 'PC percentage variation'
    )

    # Plot different PCA plots ####
    ## categorical variable ####
    if(class(colData(se.obj)[[variable]]) %in% c('character','factor')){
        ### scatter plots for individual assays  ####
        if(plot.type == 'scatter'){
            printColoredMessage(
                message = paste0(
                    '-- Create scatter PCA plots colored by the "',
                    variable, '" variable for the individual assay(s):'),
                color = 'magenta',
                verbose = verbose)
            all.scat.pca.plots <- lapply(
                levels(assay.names),
                function(x) {
                    pca.data <- all.pca.data[[x]]$u[ , seq_len(nb.pcs)]
                    pair.pcs <- combn(x = ncol(pca.data), m = 2)
                    printColoredMessage(
                        message = paste0( '- Create all possible pairwise scatter plots using the first ',
                                          nb.pcs, ' PCs for the ', x, ' data'),
                        color = 'blue',
                        verbose = verbose)
                    plot.per.data <- lapply(
                        1:ncol(pair.pcs),
                        function(i) {
                            main.plot <- ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], y = pca.data[, pair.pcs[2, i]])) +
                                geom_point(
                                    aes(fill = se.obj@colData[, variable]),
                                    color = stroke.color,
                                    pch = 21,
                                    stroke = stroke.size,
                                    size = points.size,
                                    alpha = points.alpha) +
                                scale_fill_manual(values = pca.plot.colors, name = variable) +
                                scale_x_continuous(
                                    name = paste0('PC', pair.pcs[1, i], ' (', all.pca.percentage[[x]][pair.pcs[1, i]], '%)' ),
                                    breaks = scales::pretty_breaks(n = 5)) +
                                scale_y_continuous(
                                    name = paste0('PC', pair.pcs[2, i], ' (', all.pca.percentage[[x]][pair.pcs[2, i]], '%)'),
                                    breaks = scales::pretty_breaks(n = 5)) +
                                ggtitle(x) +
                                theme_pubr() +
                                theme(
                                    legend.background = element_blank(),
                                    plot.title = element_text(size = 8),
                                    legend.text = element_text(size = 8),
                                    legend.title = element_text(size = 10, face = "bold"),
                                    axis.text.x = element_text(size = 6),
                                    axis.text.y = element_text(size = 6),
                                    axis.title.x = element_text(size = 8),
                                    axis.title.y = element_text(size = 8),
                                    plot.margin = unit(c(0,0,0,0), 'lines'),
                                    legend.position = "none")
                            dense.x <- ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], fill = se.obj@colData[[variable]])) +
                                geom_density(alpha = densities.alpha) +
                                theme_void() +
                                theme(
                                    legend.position = "none",
                                    legend.text = element_text(size = 10),
                                    legend.title = element_text(size = 12, face = "bold")) +
                                guides(fill = guide_legend(override.aes = list(size = 6, shape = 21))) +
                                scale_fill_manual(values = pca.plot.colors, name = variable)

                            dense.y <- ggplot(mapping = aes(x = pca.data[, pair.pcs[2, i]], fill = se.obj@colData[[variable]])) +
                                geom_density(alpha = densities.alpha) +
                                theme_void() +
                                theme(
                                    legend.position = "none",
                                    legend.text = element_text(size = 10),
                                    legend.title = element_text(size = 12, face = "bold")) +
                                coord_flip() +
                                guides(fill = guide_legend(override.aes = list(size = 6, shape = 21))) +
                                scale_fill_manual(values = pca.plot.colors, name = variable)

                            dense.x + plot_spacer() + main.plot + dense.y +
                                plot_layout(
                                    ncol = 2,
                                    nrow = 2,
                                    widths = c(4, 1),
                                    heights = c(1, 4))
                        })
                    return(plot.per.data)
                })
            names(all.scat.pca.plots) <- levels(assay.names)

            ### put all the scatter PCA plot for individual assay(s) ####
            printColoredMessage(
                message = '- Put all the scatter PCA plots togather for each assay.',
                color = 'blue',
                verbose = verbose
                )
            all.scat.pca.plots.assays <- lapply(
                levels(assay.names),
                function(x) {
                    ggarrange(
                        plotlist = all.scat.pca.plots[[x]],
                        common.legend = TRUE,
                        legend = legend.position,
                        nrow = ceiling(ncol(combn(x = nb.pcs, m = 2))/3)
                        )
                })
            names(all.scat.pca.plots.assays) <- levels(assay.names)

            ### put all the scatter PCA of all the assays together ####
            if(length(assay.names) > 1){
                printColoredMessage(
                    message = '-- Combine all the scatter PCA plots of all the assyas together:',
                    color = 'magenta',
                    verbose = verbose)
                overall.scat.pca.plot <- c()
                for (p in levels(assay.names)){
                    overall.scat.pca.plot <- c(
                        overall.scat.pca.plot,
                        all.scat.pca.plots[[p]])
                }
                overall.scat.pca.plot <- ggarrange(
                    plotlist = overall.scat.pca.plot,
                    common.legend = TRUE,
                    legend = legend.position,
                    nrow = plot.nrow,
                    ncol = plot.ncol
                    )
                if(class(overall.scat.pca.plot)[[1]] == 'list'){
                    plot.list <- lapply(
                        seq(length(overall.scat.pca.plot)),
                        function(x){
                            annotate_figure(
                                p = overall.scat.pca.plot[[x]],
                                top = text_grob(
                                    label = "Principal component plots",
                                    color = "orange",
                                    face = "bold",
                                    size = 18),
                                bottom = text_grob(
                                    label = paste0(
                                        'Analysis: ', 'scatter plots of the principal components.\n',
                                        'Variable ', variable, '.'),
                                    color = "black",
                                    hjust = 1,
                                    x = 1,
                                    size = 10))
                        })
                    overall.scat.pca.plot <- ggarrange(
                        plotlist = plot.list,
                        ncol = 1,
                        nrow = 1)
                } else {
                    overall.scat.pca.plot <- annotate_figure(
                        p = overall.scat.pca.plot,
                        top = text_grob(
                            label = "Principal component plots",
                            color = "orange",
                            face = "bold",
                            size = 18),
                        bottom = text_grob(
                            label = paste0(
                                'Analysis: ', 'scatter plots of the principal components.\n',
                                'Variable ', variable, '.'),
                            color = "black",
                            hjust = 1,
                            x = 1,
                            size = 10))
                }
                printColoredMessage(
                    message = '- The individual assay scatter PCA plots are combined into one.',
                    color = 'blue',
                    verbose = verbose)
                if(isTRUE(plot.output)) suppressMessages(print(overall.scat.pca.plot))
            }
        }
        ## boxplot of PCs of individual assays ####
        if (plot.type == 'boxplot'){
            printColoredMessage(
                message = paste0('-- Create boxplot of PCs acorss the "', variable, '" for the individual assay(s):'),
                color = 'magenta',
                verbose = verbose)
            all.boxplot.pca.plots <- lapply(
                levels(assay.names),
                function(x) {
                    var <- PC <- NULL
                    pca.data <- as.data.frame(all.pca.data[[x]]$u[ , seq_len(nb.pcs)])
                    colnames(pca.data) <- paste0('PC', seq_len(nb.pcs), ' (', all.pca.percentage[[x]][seq_len(nb.pcs)], '%)')
                    pca.data$var <- colData(se.obj)[, variable]
                    printColoredMessage(
                        message = paste0(
                            '- Create all possible boxplots using the first ', nb.pcs, ' PCs for the ', x, ' data'),
                        color = 'blue',
                        verbose = verbose)
                    pca.data <- tidyr::pivot_longer(
                        data = pca.data,
                        -var,
                        names_to = 'PCs',
                        values_to = 'PC')
                    pca.data <- as.data.frame(pca.data)
                    plot.p <- ggplot(pca.data, aes(x = var, y =  PC, fill = var)) +
                        geom_boxplot() +
                        stat_compare_means(method = 'anova', color = 'navy') +
                        facet_wrap(~PCs, scales = 'free') +
                        xlab(variable) +
                        ylab('PC') +
                        ggtitle(x) +
                        theme(
                            panel.background = element_blank(),
                            axis.line = element_line(colour = 'black', linewidth = 1),
                            axis.title.x = element_text(size = 10),
                            axis.title.y = element_text(size = 10),
                            plot.title = element_text(size = 12),
                            axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1),
                            axis.text.y = element_text(size = 10),
                            strip.text.x = element_text(size = 12),
                            legend.position = 'none')
                    if(isTRUE(plot.output) & length(assay.names) == 1) print(plot.p)
                    return(plot.p)
                })
            names(all.boxplot.pca.plots) <- levels(assay.names)
            ### put all boxplot PCA of assays together ####
            if(length(assay.names) > 1){
                printColoredMessage(
                    message = '-- Put all the plots of each assyas together:',
                    color = 'magenta',
                    verbose = verbose
                    )
                overall.boxplot.pca.plot <- ggarrange(
                    plotlist = all.boxplot.pca.plots,
                    nrow = plot.nrow,
                    ncol = plot.ncol
                    )
                if(class(overall.boxplot.pca.plot)[[1]] == 'list'){
                    plot.list <- lapply(
                        seq(length(overall.boxplot.pca.plot)),
                        function(x){
                            annotate_figure(
                                p = overall.boxplot.pca.plot[[x]],
                                top = text_grob(
                                    label = "Principal component plots",
                                    color = "orange",
                                    face = "bold",
                                    size = 18),
                                bottom = text_grob(
                                    label = paste0(
                                        'Analysis: ', 'boxplots of the principal components.\n',
                                        'Variable ', variable, '.'),
                                    color = "black",
                                    hjust = 1,
                                    x = 1,
                                    size = 10))
                        })
                    overall.boxplot.pca.plot <- ggarrange(
                        plotlist = plot.list,
                        ncol = 1,
                        nrow = 1)
                } else {
                    overall.boxplot.pca.plot <- annotate_figure(
                        p = overall.boxplot.pca.plot,
                        top = text_grob(
                            label = "Principal component plots",
                            color = "orange",
                            face = "bold",
                            size = 18),
                        bottom = text_grob(
                            label = paste0(
                                'Analysis: ', 'boxplots of the principal components.\n',
                                'Variable ', variable, '.'),
                            color = "black",
                            hjust = 1,
                            x = 1,
                            size = 10))
                }
                printColoredMessage(
                    message = '- The individual assay boxplot of PCs are combined into one.',
                    color = 'blue',
                    verbose = verbose)
                if(isTRUE(plot.output))
                    suppressMessages(print(overall.boxplot.pca.plot))
            }
        }
    }
    ## continuous variable ####
    if(class(colData(se.obj)[[variable]]) %in% c('numeric','integer')){
        #### scatter plots for individual assays ####
        printColoredMessage(
            message = paste0(
                '-- Create scatter plots between each PCs and the "',
                variable, '" variable for the individual assay(s):'),
            color = 'magenta',
            verbose = verbose)
        all.scat.pca.plots.assays <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('- PCA plots of for ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                var <- PC <- r.label <- NULL
                pca.data <- as.data.frame(all.pca.data[[x]]$u[ , seq_len(nb.pcs)])
                colnames(pca.data) <- paste0('PC', seq_len(nb.pcs), ' (', all.pca.percentage[[x]][seq_len(nb.pcs)], '%)')
                pca.data$var <- colData(se.obj)[[variable]]
                printColoredMessage(
                    message = paste0('- Create all possible scatter plots of the first ', nb.pcs, ' PCs.'),
                    color = 'blue',
                    verbose = verbose)
                pca.data <- tidyr::pivot_longer(
                    data = pca.data,
                    -var,
                    names_to = 'PCs',
                    values_to = 'PC')
                pca.data <- as.data.frame(pca.data)
                plot.p <- ggplot(pca.data, aes(x = var, y = PC)) +
                    geom_point(
                        color = 'gray40',
                        pch = 21,
                        stroke = .2,
                        size = 1.5,
                        alpha = points.alpha) +
                    facet_wrap(~PCs, scales = 'free') +
                    xlab(variable) +
                    ylab('PCs') +
                    ggtitle(x) +
                    geom_smooth(formula = y ~ x, method = 'lm', colour = "darkgreen")  +
                    ggpubr::stat_cor(
                        aes(label = r.label),
                        color = "navy") +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        plot.title = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8),
                        strip.text.x = element_text(size = 10),
                        legend.position = 'none')
                if(isTRUE(plot.output) & length(assay.names) == 1) print(plot.p)
                return(plot.p)
            })
        names(all.scat.pca.plots.assays) <- levels(assay.names)
        if(length(assay.names) > 1){
            printColoredMessage(
                message = '-- Put all the plots of each assyas together:',
                color = 'magenta',
                verbose = verbose
                )
            overall.scat.pca.plot <- ggarrange(
                plotlist = all.scat.pca.plots.assays,
                nrow = plot.nrow,
                ncol = plot.ncol
                )
            if(class(overall.scat.pca.plot)[[1]] == 'list'){
                plot.list <- lapply(
                    seq(length(overall.scat.pca.plot)),
                    function(x){
                        annotate_figure(
                            p = overall.scat.pca.plot[[x]],
                            top = text_grob(
                                label = "Principal component plots",
                                color = "orange",
                                face = "bold",
                                size = 18),
                            bottom = text_grob(
                                label = paste0(
                                    'Analysis: ', 'boxplots of the principal components.\n',
                                    'Variable ', variable, '.'),
                                color = "black",
                                hjust = 1,
                                x = 1,
                                size = 10))
                    })
                overall.scat.pca.plot <- ggarrange(
                    plotlist = plot.list,
                    ncol = 1,
                    nrow = 1
                    )
            } else {
                overall.scat.pca.plots <- annotate_figure(
                    p = overall.scat.pca.plot,
                    top = text_grob(
                        label = "Principal component plots",
                        color = "orange",
                        face = "bold",
                        size = 18),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ', 'boxplots of the principal components.\n',
                            'Variable ', variable, '.'),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            }
            printColoredMessage(
                message = '- The individual scatter plots of PCs are combined into one.',
                color = 'blue',
                verbose = verbose
                )
            if(isTRUE(plot.output))
                suppressMessages(print(overall.scat.pca.plot))
        }
    }
    # Save the results ####
    printColoredMessage(
        message = '-- Save the PCA plots:',
        color = 'magenta',
        verbose = verbose
        )
    ## add the pca plots to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        if(isTRUE(fast.pca)){
            method <- 'fast.svd'
        } else method <- 'ordinary.svd'
        if(plot.type == 'scatter'){
            results.data <- all.scat.pca.plots.assays
        } else results.data <- all.boxplot.pca.plots
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            slot = 'Metrics',
            assessment.type = 'global.level',
            assessment = 'PCA',
            method = method,
            file.name = paste0(plot.type, '.plot'),
            variables = variable,
            results.data = results.data
        )
        printColoredMessage(
            message = paste0('- The PCA plots of individual assay (s) are saved to the',
                             ' "se.obj@metadata$metric$AssayName$PCA$', method, '$pca.plot" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        ### overall pca plots ####
        if(length(assay.names) > 1){
            if(plot.type == 'scatter'){
                results.data <- overall.scat.pca.plot
            } else results.data <- overall.boxplot.pca.plot
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment = 'PCA',
                assessment.type = 'global.level',
                method = method,
                variables = variable,
                file.name = paste0(plot.type, '.plot'),
                plot.data = results.data
                )
            printColoredMessage(message = '------------The plotPCA function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj)
        }
    }
    ## return a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- The PCA plots are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotPCA function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            if(plot.type == 'scatter'){
                return(all.scat.pca.plots = all.scat.pca.plots)
            } else return(all.boxplot.pca.plots = all.boxplot.pca.plots)
        } else {
            if(plot.type == 'scatter'){
                return(pca.plots = list(
                    all.scat.pca.plots = all.scat.pca.plots.assays,
                    overall.scat.pca.plots = overall.scat.pca.plots))
            } else return(pca.plots = list(
                all.boxplot.pca.plots = all.boxplot.pca.plots,
                overall.boxplot.pca.plot = overall.boxplot.pca.plot))
        }
    }
}

