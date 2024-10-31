#' Plot the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function generate a dot-line plot between the first cumulative PCs and the correlation coefficient (R2) obtained
#' from the vector correlation analysis. An ideal normalization should results a low correlation with unwanted variation
#' variables and high correlation with known biology.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to generate a vector correlation. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol that indicates the name of the column in the sample annotation of the SummarizedExperiment
#' object. The variable must be a categorical variable.
#' @param fast.pca Logical. Indicates whether to use the fast PCA or PCA results computed by the computePCA function. The
#' default is 'TRUE'.
#' @param nb.pcs Numeric. A numeric value indicating the number of first PCs to use to plot the vector correlation. The
#' default is set to 10.
#' @param plot.output Logical. If 'TRUE', the vector correlation line-dot plot of individual assay(s) will be printed
#' while the function is running.
#' @param save.se.obj Logical. Indicates whether to save the vector correlation plots to the meta data of the
#' SummarizedExperiment object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all the vector correlation plots for the individual
#' assay(s).

#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @import RColorBrewer
#' @import ggplot2
#' @export

plotPCVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotPCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if(length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical varible.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop('The "variable" must have at least two levels.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
        ' Run the checkSeObj function with "remove.na = both"',
        ', then "computePCA"-->"computePCVariableCorrelation"-->"plotPCVariableCorrelation".'))
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

    # Obtain vector correlations ####
    printColoredMessage(
        message = '-- Obtain the computed vector correlations for each assay(s) from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
        )
    if(isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    all.pcs.vect.corr <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = levels(assay.names),
        slot = 'Metrics',
        assessment.type = 'global.level',
        assessment = 'VCA',
        method = method,
        variables = variable,
        file.name = 'vector.correlations',
        sub.file.name = NULL,
        required.function = 'computePCVariableCorrelation',
        message.to.print = 'vector correlation'
        )

    # Gene individual line-dot plots ####
    printColoredMessage(
        message = paste0('-- Generate the vector correlations vs. PCs line-dot plot for each assay(s):'),
        color = 'magenta',
        verbose = verbose
        )
    all.vect.corr.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Generate the vector correlations vs. PCs line-dot plot for the "', x , '" data.'),
                color = 'blue',
                verbose = verbose)
            data.to.plot <- data.frame(
                vec.corr = all.pcs.vect.corr[[x]],
                pcs = seq_len(nb.pcs)
                )
            vect.corr.plot <- ggplot(data.to.plot, aes(x = pcs, y = vec.corr, group = 1)) +
                geom_line(color = 'gray80', linewidth = 1) +
                geom_point(color = 'gray40', size = 3) +
                xlab('Cumulative PCs') +
                ylab('Vector correlations') +
                ggtitle(variable) +
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
            if(isTRUE(plot.output) & length(assay.names) == 1) print(vect.corr.plot)
            return(vect.corr.plot)
        })
    names(all.vect.corr.plots) <- levels(assay.names)

    # Generate overall plot ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = paste0('-- Put all the vector correlation line-dot plots togather:'),
            color = 'magenta',
            verbose = verbose
            )
        datasets <- pcs <- vec.corr <- NULL
        all.pcs.vect.corr <- as.data.frame(do.call(cbind, all.pcs.vect.corr))
        all.pcs.vect.corr <- all.pcs.vect.corr %>%
            mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'datasets',
                values_to = 'vec.corr') %>%
            dplyr::mutate(datasets = factor(datasets, levels = levels(assay.names)))
        overall.vect.corr.plot <- ggplot(all.pcs.vect.corr, aes(x = pcs, y = vec.corr, group = datasets)) +
            geom_line(aes(color = datasets), size = 1) +
            geom_point(aes(color = datasets), size = 3) +
            xlab('Cumulative PCs') +
            ylab('Vector correlations') +
            ggtitle(variable) +
            scale_color_manual(values = c(data.sets.colors), name = 'Datasets') +
            scale_x_continuous(breaks = seq_len(nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                plot.title = element_text(size = 16),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                axis.text.x = element_text(size = 11, angle = 35, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14))
        overall.vect.corr.plot <- annotate_figure(
            p = overall.vect.corr.plot,
            top = text_grob(
                label = "Vector correlation analsis",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ', 'vector correlation analysis between the first ', nb.pcs, ' PCs and the ', variable, ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
            )
        printColoredMessage(
            message = '- The individual assay vector correlation line-dots plot are combined into one.',
            color = 'blue',
            verbose = verbose)
        if(isTRUE(plot.output)) suppressMessages(print(overall.vect.corr.plot))
    }

    # Save the plots ####
    printColoredMessage(
        message = '-- Save all the vector correlation line-dots plots:',
        color = 'magenta',
        verbose = verbose)
    if (save.se.obj == TRUE) {
        ## add results to the SummarizedExperiment object ####
        printColoredMessage(
            message = '- Save all the vector correlation line-dots plot(s) to the "metadata" of the SummarizedExperiment object.',
            color = 'orange',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'VCA',
            method = method,
            variables = variable,
            file.name = 'plot',
            results.data = all.vect.corr.plots
            )

        printColoredMessage(
            message = paste0('* The vector correlation line-dot plot for the individal assay(s) is saved to the ',
                             ' "se.obj@metadata$metric$AssayName$VCA$VariableName$vect.corr.plot" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)

        ## add overall vector correlation plot of all assays ####
        if(length(assay.names) > 1){
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'global.level',
                assessment = 'VCA',
                method = method,
                variables = variable,
                file.name = 'line.dotplot',
                plot.data = overall.vect.corr.plot
                )
            printColoredMessage(
                message = paste0('* The combined vector correlation line-dots plots of all the assays are saved to the',
                                 ' "se.obj@metadata$plot$VCA$VaribleName" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotPCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- The vector correlation plots are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotPCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.vec.corr = list(all.vect.corr.plots = all.vect.corr.plots))
        } else {
            return(all.vec.corr = list(
                all.vect.corr.plots = all.vect.corr.plots,
                overall.vect.corr.plot = overall.vect.corr.plot))
        }
    }
}

