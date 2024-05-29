#' Generates study outline of a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function generates a heatmap showing the study outline.

#' @param se.obj A SummarizedExperiment object.
#' @param variables Symbol. A symbol or symbols representing the label of variable(s) within the SummarizedExperiment
#' object. This can comprise a vector containing either categorical, continuous, or a combination of both variables.
#' @param sort.variable Symbol. A symbol indicating the label of variable that is used to sort the sample information.
#' The default is 'NULL'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the 'checkSeObj()'
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param plot.output Logical. Determines whether to plot the study outline. The default is set to 'TRUE'.
#' @param save.se.obj Logical. Indicates whether to save the study outline plot in the metadata of the SummarizedExperiment
#' object or to output the result as a plot. The default is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap Heatmap draw ht_opt
#' @importFrom dplyr arrange pick
#' @import viridis
#' @export

plotStudyOutline <- function(
        se.obj,
        variables,
        sort.variable = NULL,
        assess.se.obj = TRUE,
        remove.na = 'none',
        plot.output = TRUE,
        legend.font.size = 14,
        legend.ncol = 4,
        legend.direction = 'horizontal',
        heatmap.legend.side = 'bottom',
        column.names.rot = 25,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Sample information ####
    sample.info <- as.data.frame(colData(x = se.obj)[ , variables, drop = FALSE])
    if(!is.null(sort.variable))
        sample.info <- sample.info %>%  arrange(dplyr::pick(sort.variable))

    var.class <- sapply(variables, function(x) class(sample.info[ , x]))
    cat.var <- var.class %in% c('factor', 'character')
    cont.var <- var.class %in% c('integer', 'numeric')

    # Select colors ####
    ## Categorical
    palette.colors <- brewer.pal.info
    palette.colors <- palette.colors[order(palette.colors$category, decreasing = FALSE) , ]
    names(variables)[cat.var] <- seq(sum(cat.var))

    ## Continuous
    viridis.colors <- c('plasma', 'cividis', 'inferno', 'magma', 'mako','rocket', 'turbo', 'viridis')
    viridis.colors <- rep(viridis.colors, c(ceiling(sum(cont.var) /8 )))
    names(variables)[cont.var] <- seq(sum(cont.var))

    # Generate heatmaps ####
    ht.list <- NULL
    for(i in 1:length(variables)){
        selected.variable <- variables[i]
        if(class(sample.info[ , selected.variable]) %in% c('factor', 'character') ){
            num.colors <- length(unique(sample.info[ , selected.variable]))
            colfunc <- grDevices::colorRampPalette(
                RColorBrewer::brewer.pal(
                    n = palette.colors$maxcolors[as.numeric(names(selected.variable))],
                    name = row.names(palette.colors)[as.numeric(names(selected.variable))]))
            color.plates <- colfunc(num.colors)
            color.plates <- color.plates[as.integer(as.factor(sample.info[ , selected.variable]))]
            names(color.plates) <- sample.info[ , selected.variable]
            ht.list <-  ht.list + ComplexHeatmap::Heatmap(
                sample.info[ , selected.variable],
                cluster_rows = FALSE,
                name = selected.variable,
                column_names_rot = column.names.rot,
                col = color.plates,
                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend.font.size), by_row = TRUE, ncol = legend.ncol))
        } else if (class(sample.info[ , selected.variable])  %in% c('integer', 'numeric')){
            ht.list <-  ht.list + ComplexHeatmap::Heatmap(
                sample.info[ , selected.variable],
                cluster_rows = FALSE,
                name = selected.variable,
                column_names_rot = column.names.rot,
                col = viridis(n = 10, option = viridis.colors[as.numeric(names(selected.variable))]),
                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend.font.size), legend_direction = legend.direction))
        }
    }
    if(isTRUE(plot.output))
        ht_opt$message = FALSE
        print(draw(
            object = ht.list,
            heatmap_legend_side = heatmap.legend.side,
            auto_adjust = TRUE))
    # Save the data ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save the study outline plot:',
        color = 'magenta',
        verbose = verbose)
    ## check if metadata metric already exist
    if(isTRUE(save.se.obj)){
        if (!'StudyOutline' %in% names(se.obj@metadata)) {
            se.obj@metadata[['StudyOutline']] <- list()
        }
        se.obj@metadata[['StudyOutline']] <- draw(
            object = ht.list,
            heatmap_legend_side = heatmap.legend.side,
            auto_adjust = TRUE)
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(ht.list)
    }
}


