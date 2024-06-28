#' Plot gene-gene partial correlation.

#' @author Ramyar Molania

#' @description
#' This function computes all possible gene-gene pairwise ordinary and partial correlation of the assays in the
#' SummarizedExperiment object.

#' @details
#' Partial correlation is used to estimate correlation between two variables while controlling for third
#' variables.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to obtain the computed correlations. The default is set to "all, which
#' indicates all the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. A symbol that indicates the column name of the SummarizedExperiment object that has been used
#' to compute partial correlation.
#' @param method Symbol. A symbol that indicates which correlation methods has been computed for the "variable". The options
#' are 'pearson' or "spearman". The default is set to 'spearman'.
#' @param plot.type Symbol. A symbol specifying how to plot the results of the partial pairwise gene correlation. The
#' options are:
#' @param filter.genes Logical. Indicates to filter genes before generating the plots or not. The default is set to 'TRUE'.
#' @param corr.dif.cutoff Numeric. A numeric value to be used as a cutoff for filtering genes. The default is set to
#' 0.4.
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid. When the number of selected assay is more than
#' 1, the function puts all the RLE boxplots in one grid.
#' @param plot.nrow Numeric. Indicates number of rows in the plot grid. When the number of selected assay is more than
#' 3, the function puts all the RLE boxplots in one grid.
#' @param plot.output Logical. Indicates to print the histogram of partial correlation coefficients or not. The default is
#' set to 'FLASE'.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object or
#' to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return The SummarizedExperiment object that contains the correlation plots or a list of the correlation plots
#' for individual assays.

#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggjoy geom_joy
#' @import ggplot2
#' @export

plotGenesPartialCorrelation <- function(
        se.obj,
        assay.names = "all",
        variable = NULL,
        method = 'spearman',
        plot.type = 'barplot',
        filter.genes = TRUE,
        corr.dif.cutoff = 0.4,
        plot.ncol = 2,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGenesPartialCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (!is.null(variable)){
        if (length(variable) > 1){
            stop('The "variable" must contain only one variable.')
        }
        if (!class(colData(se.obj)[[variable]]) %in% c('numeric', 'integer')){
            stop('The "variable" must be a continous variable.')
        }
        if (sum(is.na(se.obj@colData[[variable]])) > 0){
            stop(paste0('The "', variable, '" variable contains NA. ',
                        'Run the checkSeObj function with "remove.na = both"',
                        ', then re-run the computeRLE function.'))
        }
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain the correlations data ####
    printColoredMessage(
        message = paste0('-- Obtain the computed correlation data from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    all.corr.data <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = assay.names,
        assessment = 'PPcorr',
        assessment.type = 'gene.level',
        method = method,
        variables = variable,
        file.name = 'correlations',
        sub.file.name = NULL,
        required.function = 'computeGenesPartialCorrelation',
        message.to.print = 'partial correlation'
        )

    # Filter the correlation
    if(isTRUE(filter.genes)){
        all.corr.data <- lapply(
            levels(assay.names),
            function(x) {
                corr.dif <- all.corr.data[[x]]$p.cor - all.corr.data[[x]]$pp.cor
                select.genes <- abs(corr.dif) > corr.dif.cutoff
                if(sum(select.genes) == 0){
                    printColoredMessage(
                        message = paste0('- All gene-gene correlations are filtered for the "', x, '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    all.corr.data[[x]] <- NULL
                } else if (sum(select.genes) != 0){
                    printColoredMessage(
                        message = paste0('- ', nrow(all.corr.data[[x]]) - sum(select.genes)  , ' gene-gene correlation are filtered for "', x, '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    all.corr.data[[x]][select.genes , ]
                }
            })
        names(all.corr.data) <- levels(assay.names)
        all.corr.data <- Filter(Negate(is.null), all.corr.data)
        if(length(all.corr.data) == 0)
            stop('Any assay contain gene-gene correlations for plotting.')
        assay.names <- droplevels(assay.names[assay.names %in% names(all.corr.data)])
    }

    # Generate different plots of partial and ordinary correlations ####
    ## generate histograms for each assay  ####
    if(plot.type == 'histogram'){
        printColoredMessage(
            message = '-- Generate the histograms of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.hist <- lapply(
            levels(assay.names),
            function(x){
                corr.data <- as.data.frame(all.corr.data[[x]][ , c('p.cor', 'pp.cor')])
                colnames(corr.data) <- c('Correlation', 'Partial correlation')
                corr.data <- pivot_longer(corr.data, everything(), names_to = 'corr.type', values_to = 'corr')
                corr <- corr.type <- NULL
                p.joy <- ggplot(data = corr.data, aes(x = corr, y = corr.type)) +
                    ggjoy::geom_joy(color = 'gray60') +
                    xlim(-1,1) +
                    ggtitle(x) +
                    xlab('Correlation coefficients') +
                    ylab('') +
                    labs(caption = paste0(
                        'Analysis: ', 'histograms of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ", variable)) +
                    theme_minimal() +
                    theme(
                        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
                        plot.caption = element_text(hjust = 0, vjust = 0),
                        plot.title = element_text(size = 16),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14, angle = 65),
                    )
                if(isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.joy)
                return(p.joy)
            })
        names(all.corr.hist) <- levels(assay.names)
        ### overall histograms of all assays ####
        if(length(levels(assay.names)) > 1){
            for(i in levels(assay.names) ){
                all.corr.hist[[i]]$labels$caption <- NULL
            }
            overall.hist.plots <- ggarrange(
                plotlist = all.corr.hist ,
                ncol = plot.ncol,
                nrow = plot.nrow
            )
            overall.hist.plots <- annotate_figure(
                p = overall.hist.plots,
                top = text_grob(
                    label = "Partial pairwise correlation",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ', 'histograms of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ", variable),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10))
            if(isTRUE(plot.output)) print(overall.hist.plots)
        }
    }

    ## generate scatter plots for each assay  ####
    if(plot.type == 'scatter.plot'){
        printColoredMessage(
            message = '-- Generate the scatter plots of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.scatter <- lapply(
            levels(assay.names),
            function(x){
                corr.data <- as.data.frame(all.corr.data[[x]][ , c('p.cor', 'pp.cor')])
                Correlation <- Partial.correlation <- NULL
                colnames(corr.data) <- c('Correlation', 'Partial.correlation')
                p.hex <- ggplot(data = corr.data, aes(x = Correlation, y = Partial.correlation)) +
                    geom_hex(alpha = 0.6) +
                    xlim(-1,1) +
                    ylim(-1,1) +
                    ggtitle(x) +
                    geom_abline(intercept = 0, slope = 1, color = "red") +
                    xlab('Correlations') +
                    ylab('Partial correlations') +
                    labs(caption = paste0(
                        'Analysis: ', 'scatter plots of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ", variable)) +
                    theme_minimal() +
                    theme(
                        plot.caption = element_text(hjust = 0, vjust = 0),
                        plot.title = element_text(size = 16),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12))
                if(isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.hex)
                return(p.hex)
            })
        names(all.corr.scatter) <- levels(assay.names)

        ### overall scatter plots of all assays ####
        if(length(levels(assay.names)) > 1){
            for(i in levels(assay.names) ){
                all.corr.scatter[[i]]$labels$caption <- NULL
            }
            overall.scatter.plots <- ggarrange(
                plotlist = all.corr.scatter ,
                common.legend = TRUE,
                ncol = plot.ncol,
                nrow = plot.nrow
            )
            overall.scatter.plots <- annotate_figure(
                p = overall.scatter.plots,
                top = text_grob(
                    label = "Partial pairwise correlation",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ', 'scatter plots of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ", variable),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10))
            if(isTRUE(plot.output)) print(overall.scatter.plots)
        }
    }
    ## generate barplot for each assay ####
    if(plot.type == 'barplot'){

        printColoredMessage(
            message = '-- Generate the barplot of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.barplots <- lapply(
            levels(assay.names),
            function(x) {
                corr.dif <- all.corr.data[[x]]$p.cor - all.corr.data[[x]]$pp.cor
                select.genes <- abs(corr.dif) > corr.dif.cutoff
                p.barplot <- ggplot() +
                    geom_col(aes(x = 1, y = sum(select.genes))) +
                    ggtitle(x) +
                    xlab(x) +
                    ylab('Number of genes') +
                    labs(caption = paste0(
                        'Analysis: ', 'barplot of the number of genes that show difference between partial and ordinary correlations\n',
                        "Variable: ", variable)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12))
                if(isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.barplot)
                return(p.barplot)
            })
        names(all.corr.barplots) <- levels(assay.names)

        ### generate barplot for all assays ####
        datasets <- corr.dif.group <- dd <- NULL
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('-- Put all the barplots togather:'),
                color = 'magenta',
                verbose = verbose)
            all.corr.data <- lapply(
                levels(assay.names),
                function(x) {
                    temp.data <- all.corr.data[[x]][ , c('p.cor', 'pp.cor')]
                    temp.data$datasets <- rep(x, nrow(temp.data))
                    return(temp.data)
                })
            all.corr.data <- do.call(rbind, all.corr.data)
            all.corr.data$corr.dif <- all.corr.data$p.cor - all.corr.data$pp.cor
            all.corr.data$corr.dif.group <- ifelse(test = all.corr.data$corr.dif > corr.dif.cutoff, TRUE, FALSE)
            all.corr.data <- group_by(.data = all.corr.data, datasets)
            all.corr.data <- summarise(.data = all.corr.data, dd = sum(corr.dif.group))
            all.corr.data <- as.data.frame(all.corr.data)
            all.corr.data$datasets <- factor(
                x =  all.corr.data$datasets,
                levels = levels(assay.names))
            overall.barplots <- ggplot(all.corr.data, aes(x = datasets, y = dd)) +
                geom_col() +
                xlab('Datasets') +
                ylab('Number of genes') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(
                        size = 10,
                        angle = 25,
                        hjust = 1),
                    axis.text.y = element_text(size = 12))
            overall.barplots <- annotate_figure(
                p = overall.barplots,
                top = text_grob(
                    label = "Partial pairwise correlation",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ', 'number of paired genes that show at least ', corr.dif.cutoff,
                        ' difference between partial and ordinary correlation\n',
                        "Variable: ", variable),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
                )
            printColoredMessage(
                message = '- The individual assay ARI barplot are combined into one.',
                color = 'blue',
                verbose = verbose
                )
            if(isTRUE(plot.output))
                suppressMessages(print(overall.barplots))
        }
    }

    # Save the plots ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save all the plots :',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save all the plots to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        # for all assays

        if(plot.type == 'barplot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'barplot',
                results.data = all.corr.barplots
            )
        } else if(plot.type == 'histogram'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'histogram',
                results.data = all.corr.hist
            )
        } else if (plot.type == 'scatter.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'scatter.plot',
                results.data = all.corr.scatter
            )
        }
        printColoredMessage(
            message = paste0('- The plots of individual assay is saved to the metadata@metric in SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
        )
        # Overall plot
        if(length(levels(assay.names)) > 1){
            if(plot.type == 'barplot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'barplot',
                    plot.data = overall.barplots
                )
            } else if(plot.type == 'histogram'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'histogram',
                    plot.data = overall.hist.plots
                )
            } else if (plot.type == 'scatter.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'scatter.plot',
                    plot.data = overall.scatter.plots
                )
            }
            printColoredMessage(
                message = paste0('- The combined plots of all assays are saved to the',
                                 ' "se.obj@metadata$plot$RLE" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotGenesPartialCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }

    ## output the plots as alist ####
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = paste0('- The plots of individual assay is outputed as list.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The plotGenesPartialCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        if(plot.type == 'barplot'){
            if(length(levels(assay.names)) > 1){
                return(list(
                    all.corr.barplots = all.corr.barplots,
                    overall.barplots = overall.barplots)
                    )
            } else return(all.corr.barplots = all.corr.barplots)
        }
        if(plot.type == 'histogram'){
            if(length(levels(assay.names)) > 1){
                return(list(
                    all.corr.hist = all.corr.hist,
                    overall.hist.plots = overall.hist.plots)
                    )
            } else return(all.corr.hist = all.corr.hist)
        }
        if(plot.type == 'scatter.plot'){
            if(length(levels(assay.names)) > 1){
                return(list(
                    all.corr.scatter = all.corr.scatter,
                    overall.scatter.plots = overall.scatter.plots)
                    )
            } else return(all.corr.scatter = all.corr.scatter)
        }
    }
}


