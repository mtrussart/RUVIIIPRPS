#' Generate line-dot plot of the PC variable regression analysis.

#' @author Ramyar Molania

#' @description
#' This function generate a dot-line plot between the first cumulative PCs and the regression R squared obtained
#' from the regression analysis. An ideal normalization should results a low R squared with unwanted variation
#' variables and high R squared with known biology.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to generate a vector correlation. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol indicating the name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable must be a continuous variable such as library size, tumor purity, ... .
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. A numeric value specifying the number of first PCs to be used to plot the r squared.
#' The default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the correlation statistics, by default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object containing the line-dot plots for the continuous variable and if requested
#'  the associated plot.

#' @importFrom tidyr pivot_longer
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @export

plotPCVariableRegression <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotPCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty')
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Select colors ####
    if(length(levels(assay.names)) < 9 ){
        data.sets.colors <- RColorBrewer::brewer.pal(8, 'Dark2')[1:length(levels(assay.names))]
        names(data.sets.colors) <- levels(assay.names)
    } else {
        colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, 'Dark2'))
        data.sets.colors <- colfunc(n = length(levels(assay.names)))
        names(data.sets.colors) <- levels(assay.names)
    }

    # Obtain regression r squared ####
    printColoredMessage(
        message = '-- Obtain the computed R squared of linear regressiob for each assay(s) from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
    )
    if(isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    all.reg.rseq <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = levels(assay.names),
        assessment = 'LRA',
        assessment.type = 'global.level',
        method = method,
        variables = variable,
        file.name = 'r.squared',
        sub.file.name = NULL,
        required.function = 'computePCVariableRegression',
        message.to.print = 'computed R squared'
        )
    # Gene individual line-dot plots ####
    printColoredMessage(
        message = paste0('-- Generate the R squared vs. PCs line-dot plot for each assay(s):'),
        color = 'magenta',
        verbose = verbose)
    rseq <- NULL
    all.reg.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Generate the R squared vs. PCs line-dot plot for the "', x , '" data.'),
                color = 'blue',
                verbose = verbose)
            to.plot <- data.frame(
                rseq = all.reg.rseq[[x]],
                pcs = seq_len(nb.pcs))
            reg.plot <- ggplot(to.plot, aes(x = pcs, y = rseq, group = 1)) +
                geom_line(color = 'gray80', size = 1) +
                geom_point(color = 'gray40', size = 3) +
                ggtitle(variable) +
                xlab('Cumulative PCs') +
                ylab(expression(R^2)) +
                scale_x_continuous(breaks = seq_len(nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    plot.title = element_text(size = 16),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                    axis.text.y = element_text(size = 12))
            if(isTRUE(plot.output) & length(assay.names) == 1) print(reg.plot)
            return(reg.plot)
        })
    names(all.reg.plots) <- levels(assay.names)

    # overall plots ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = paste0('-- Put all the R squared vs. PCs line-dot plots togather:'),
            color = 'magenta',
            verbose = verbose)
        datasets <- pcs <- vec.corr <- NULL
        all.reg.rseq <- as.data.frame(do.call(cbind, all.reg.rseq))
        all.reg.rseq <- all.reg.rseq %>%
            mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'datasets',
                values_to = 'vec.corr') %>%
            mutate(datasets = factor(datasets, levels = levels(assay.names)))
        overall.reg.plot <- ggplot(all.reg.rseq, aes(x = pcs, y = vec.corr, group = datasets)) +
            geom_line(aes(color = datasets), size = 1) +
            geom_point(aes(color = datasets), size = 3) +
            xlab('Cumulative PCs') +
            ylab(expression(R^2)) +
            ggtitle(variable) +
            scale_color_manual(values = c(data.sets.colors), name = 'Datasets') +
            scale_x_continuous(breaks = seq_len(nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                plot.title = element_text(size = 16),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 12))
        overall.reg.plot <- annotate_figure(
            p = overall.reg.plot,
            top = text_grob(
                label = "Linear regression analysis",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ', 'linear regression analysis between the first ', nb.pcs, ' PCs and the ', variable, ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
        )
        printColoredMessage(
            message = '- The individual assay the R squared vs. PCs line-dot plots are combined into one.',
            color = 'blue',
            verbose = verbose)
        if(isTRUE(plot.output))
            suppressMessages(print(overall.reg.plot))
    }

    # Save the plots ####
    printColoredMessage(
        message = '-- Save all the R squared vs. PCs line-dot plots:',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save all the R squared vs. PCs line-dot plot(s) to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## add vector correlation plots of individual assays ####
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'LRA',
            method = method,
            variables = variable,
            file.name = 'plot',
            results.data = all.reg.plots
        )

        printColoredMessage(
            message = paste0('- The R squared vs. PCs line-dot plot for the individal assay(s) is saved to the ',
                             ' "se.obj@metadata$metric$AssayName$LRA$VariableName$lra.plot" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)

        ## add overall vector correlation plot of all assays ####
        if(length(assay.names) > 1){
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'global.level',
                assessment = 'LRA',
                method = method,
                variables = variable,
                file.name = 'line.dotplot',
                plot.data = overall.reg.plot
            )
            printColoredMessage(
                message = paste0('- The combined R squared vs. PCs line-dot plots of all the assays are saved to the',
                                 ' "se.obj@metadata$plot$LRA$VaribleName" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotPCVariableRegressions function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## Return only the correlation result
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '-- Save the plots as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotPCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.pc.var.reg.plots = list(all.reg.plots = all.reg.plots))
        } else{
            return(all.pc.var.reg.plots = list(
                all.reg.plots = all.reg.plots,
                overall.reg.plot = overall.reg.plot))
        }

    }
}

