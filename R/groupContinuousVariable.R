#' Group a continuous variable.

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param variable Symbol. A symbol specifying the column in the SummarizedExperiment object. This should be a continuous
#' variable.
#' @param nb.clusters Numeric. A numeric value indicating the number of groups for continuous sources of unwanted variation.
#' The default is 3. This implies that each continuous source will be split into 3 groups using the specified
#' 'clustering.method' method.
#' @param clustering.method Symbol. A symbol specifying the clustering method to be applied for grouping each continuous
#' source of unwanted variables. Options include 'kmeans', 'cut', and 'quantile'. The default is 'kmeans' clustering.
#' @param perfix Symbol. A symbol to add to each groups.
#' @param plot.output Logical. If is 'TRUE', the function plots the variable values colored by the groups.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

groupContinuousVariable <- function(
        se.obj,
        variable,
        nb.clusters = 3,
        clustering.method = 'kmeans',
        perfix = '_group',
        plot.output = TRUE,
        verbose = TRUE
        ){
    # checking the function input ####
    if(!is.numeric(se.obj[[variable]])){
        stop(paste0('The variable ', variable, 'must be a continuous variable.'))
    }
    if(var(se.obj[[variable]]) == 0){
        stop(paste0('The variance of the variable ', variable, ' is 0. This variable must contain variation.'))
    }

    # grouping the variable ####
    initial.values <- colData(se.obj)[[variable]]
    printColoredMessage(
        message = paste0(
            '- The "',
            variable,
            '" is a continouse variable, then this will be divided into ',
            nb.clusters,
            ' groups using the ',
            clustering.method,
            ' clustering.'),
        color = 'blue',
        verbose = verbose
    )
    # kmeans
    if (clustering.method == 'kmeans') {
        set.seed(3456)
        kmeans.clusters <- kmeans(
            x = colData(se.obj)[[variable]],
            centers = nb.clusters,
            iter.max = 1000
        )
        groups <- factor(x = paste0(variable, perfix, kmeans.clusters$cluster))
    }
    # cut
    if (clustering.method == 'cut') {
        cut.clusters <- as.numeric(cut(
            x = colData(se.obj)[[variable]],
            breaks = nb.clusters,
            include.lowest = TRUE
        ))
        groups <- factor(x = paste0(variable, perfix, cut.clusters))
    }
    # quantile
    if (clustering.method == 'quantile') {
        quantiles <- quantile(
            x = colData(se.obj)[[variable]],
            probs = seq(0, 1, 1 / nb.clusters)
        )
        quantiles.clusters <- as.numeric(cut(
            x = colData(se.obj)[[variable]],
            breaks = quantiles,
            include.lowest = TRUE
        )
        )
        groups <- factor(x = paste0(variable, perfix, quantiles.clusters))
    }

    # plot the output ###
    p.plot <- data.frame(
        variable = initial.values,
        samples.nb = 1:ncol(se.obj),
        Groups = groups) %>%
        ggplot(data = ., aes(x = samples.nb, y = variable, color = Groups )) +
        geom_point() +
        xlab('Samples') +
        ylab(variable) +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),

        )
    if(isTRUE(plot.output)) print(p.plot)
    return(groups)
}




