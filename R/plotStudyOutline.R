#' Generates study outline of a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function generates a heatmap of sample-level features such as batches, biological populations, library size, and
#' tumor purity. This plot is helpful for exploring how these factors are distributed across samples and for examining
#' visible unwanted variation in the data. Further, it can visually reveal confounder factors in the data. It can be initially
#' generated based on the known variables and updated with any estimated biological and unwanted variables during the
#' normalization process. We highly recommend generating the plot before performing downstream analyses.

#' @param se.obj A SummarizedExperiment object.
#' @param variables Character or character vector. The label(s) of variable(s) within the SummarizedExperiment
#' object. This can include a vector containing categorical, continuous, or a combination of both types of variables. This
#' cannot be empty or NULL.
#' @param variable.to.sort Character. The label of the variable used to sort the sample information. The default is 'NULL'.
#' This means the samples will be plotted as they are.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. See the 'checkSeObj()'
#' function for more details. The default is set to TRUE.
#' @param remove.na Character. Specifies whether to remove NA or missing values from the assays. The options are 'assays'
#' , 'sample.annotation', 'both' or 'none'. The default is set to "none". See the 'checkSeObj' function for more details.
#' @param plot.output Logical. Determines whether to plot the study outline. The default is set to 'TRUE'.
#' @param legend.font.size Numeric. A numeric value indicating the size of the font of the legend in the plot. The default
#' is set to 14.
#' @param legend.ncol Numeric. A numeric value indicating the number of columns in the legend of the heatmap. The default
#' is set to 4.
#' @param legend.direction Character. A string specifying the direction of the legend in the heatmap. The options are
#' 'horizontal' or 'vertical'. The default is set to 'horizontal'.
#' @param heatmap.legend.side Character. A string indicating the side of the legend in the heatmap. The options are
#' 'right', 'bottom', 'left', 'top'. The default is set to 'bottom'.
#' @param column.names.rot Numeric. A numeric value indicating the angle of the column labels in the heatmap. The default
#' is set to 25.
#' @param save.se.obj Logical. Indicates whether to save the study outline plot in the metadata of the SummarizedExperiment
#' object or output the result as a plot. The default is set to 'TRUE'. The plot will be saved in 'se.obj->metadata->StudyOutline'.
#' @param verbose Logical. If 'TRUE', shows the messages for different steps of the function.

#' @return Either a SummarizedExperiment object containing a study outline plot or just the plot.

#' @importFrom ComplexHeatmap Heatmap draw ht_opt
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr arrange pick
#' @import RColorBrewer
#' @import viridis
#' @export

plotStudyOutline <- function(
        se.obj,
        variables,
        variable.to.sort = NULL,
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
    printColoredMessage(message = '------------The plotStudyOutline function starts:',
                        color = 'white',
                        verbose = verbose)
    # checking the function input
    if (is.null(se.obj) | is.logical(se.obj)){
        stop('The "se.obj" cannot be NULL or logical. A SummarizedExperiment object must be provided.')
    }
    if (!class(se.obj)[1] %in% c('SummarizedExperiment', 'RangedSummarizedExperiment')){
        stop('The "se.obj" object is not a class of  SummarizedExperiment object.')
    }
    if (is.null(variables) | is.logical(variables) | is.list(variables)){
        stop('The "variables" cannot be NULL or logical or list. Pleas see the description of the argument.')
    }
    if(!sum(variables %in% colnames(colData(se.obj))) == length(variables)){
        not.found.variables <- variables[!variables %in% colnames(colData(se.obj))]
        stop(paste0(
            'The variable(s) ',
            paste(not.found.variables, collapse = ', '),
            ' cannot be found in the SummarizedExperiment object.'))
    }
    if (!variable.to.sort %in% colnames(colData(se.obj))){
        stop('The "variable.to.sort" cannot be found in the SummarizedExperiment object.')
    }
    if (remove.na %in% c('assays', 'sample.annotation', 'both', 'none')){
        stop('The "remove.na" must be one of the "assays", "sample.annotation", "both" or "none".')
    }
    if (legend.direction %in% c('horizontal', 'vertical')){
        stop('The "legend.direction" must be one of the "horizontal" or "vertical".')
    }
    if (heatmap.legend.side %in% c('right', 'top', 'bottom', 'left')){
        stop('The "legend.direction" must be one of the "right", "top", "bottom" or "left".')
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical.')
    }

    # Sample information ####
    sample.info <- as.data.frame(colData(x = se.obj)[ , variables, drop = FALSE])
    colnames(sample.info) <- variables
    if(!is.null(variable.to.sort))
        sample.info <- sample.info %>%
        arrange(dplyr::pick(variable.to.sort))

    var.class <- sapply(variables, function(x) class(sample.info[[x]]))
    cat.var <- var.class %in% c('factor', 'character')
    cont.var <- var.class %in% c('integer', 'numeric')

    # # Select colors ####
    # ## Categorical
    # palette.colors <- brewer.pal.info
    # palette.colors <- palette.colors[order(palette.colors$category, decreasing = FALSE) , ]
    # names(variables)[cat.var] <- seq(sum(cat.var))
    #
    # ## Continuous
    viridis.colors <- c('plasma', 'cividis', 'inferno', 'magma', 'mako','rocket', 'turbo', 'viridis')
    viridis.colors <- rep(viridis.colors, c(ceiling(sum(cont.var) /8 )))
    names(variables)[cont.var] <- seq(sum(cont.var))
    # # Generate heatmaps ####
    # ht.list <- NULL
    # for(i in 1:length(variables)){
    #     selected.variable <- variables[i]
    #     if(class(sample.info[ , selected.variable]) %in% c('factor', 'character') ){
    #         num.colors <- length(unique(sample.info[ , selected.variable]))
    #         colfunc <- grDevices::colorRampPalette(
    #             RColorBrewer::brewer.pal(
    #                 n = palette.colors$maxcolors[as.numeric(names(selected.variable))],
    #                 name = row.names(palette.colors)[as.numeric(names(selected.variable))]
    #                 )
    #             )
    #         color.plates <- colfunc(num.colors)
    #         color.plates <- color.plates[as.integer(as.factor(sample.info[ , selected.variable]))]
    #         names(color.plates) <- sample.info[ , selected.variable]
    #         ht.list <-  ht.list + ComplexHeatmap::Heatmap(
    #             sample.info[ , selected.variable],
    #             cluster_rows = FALSE,
    #             name = selected.variable,
    #             column_names_rot = column.names.rot,
    #             col = color.plates,
    #             heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend.font.size), by_row = TRUE, ncol = legend.ncol))
    #     } else if (class(sample.info[ , selected.variable])  %in% c('integer', 'numeric')){
    #         ht.list <-  ht.list + ComplexHeatmap::Heatmap(
    #             sample.info[ , selected.variable],
    #             cluster_rows = FALSE,
    #             name = selected.variable,
    #             column_names_rot = column.names.rot,
    #             col = viridis(n = 10, option = viridis.colors[as.numeric(names(selected.variable))]),
    #             heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend.font.size), legend_direction = legend.direction))
    #     }
    # }
    currentCols <-  c(
        RColorBrewer::brewer.pal(8, "Dark2")[-5],
        RColorBrewer::brewer.pal(10, "Paired"),
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(10, "Paired"))
    currentCols <- rep(currentCols, 3)

    # Generate heatmaps ####
    ht.list <- NULL
    for(i in 1:length(variables)){
        selected.variable <- variables[i]
        if(class(sample.info[ , selected.variable]) %in% c('factor', 'character') ){
            num.colors <- length(unique(sample.info[ , selected.variable]))
            # colfunc <- grDevices::colorRampPalette(
            #     RColorBrewer::brewer.pal(
            #         n = palette.colors$maxcolors[as.numeric(names(selected.variable))],
            #         name = row.names(palette.colors)[as.numeric(names(selected.variable))]
            #     )
            # )
            color.plates <- currentCols[seq(num.colors)]
            color.plates <- color.plates[as.integer(as.factor(sample.info[ , selected.variable]))]
            names(color.plates) <- sample.info[ , selected.variable]
            ht.list <-  ht.list + ComplexHeatmap::Heatmap(
                sample.info[ , selected.variable],
                cluster_rows = FALSE,
                name = selected.variable,
                column_names_rot = column.names.rot,
                col = color.plates,
                heatmap_legend_param = list(
                    title_gp = grid::gpar(fontsize = legend.font.size),
                    by_row = TRUE,
                    ncol = legend.ncol)
                )
            currentCols <- currentCols[-seq(num.colors)]
        } else if (class(sample.info[ , selected.variable])  %in% c('integer', 'numeric')){
            ht.list <-  ht.list + ComplexHeatmap::Heatmap(
                sample.info[ , selected.variable],
                cluster_rows = FALSE,
                name = selected.variable,
                column_names_rot = column.names.rot,
                col = viridis(n = 10, option = viridis.colors[as.numeric(names(selected.variable))]),
                heatmap_legend_param = list(
                    title_gp = grid::gpar(fontsize = legend.font.size),
                    legend_direction = legend.direction)
                )
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
        printColoredMessage(message = '------------The plotStudyOutline function finished.:',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The plotStudyOutline function finished.:',
                            color = 'white',
                            verbose = verbose)
        return(ht.list)
    }
}


