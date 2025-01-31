#' Create all possible assessment plots and metrics for the variables.

#' @author Ramyar Molania

#' @description
#' This functions provides the names of all possible assessment plots and metrics for the given variable(s). The list
#' will be used in the 'assessVariation' and 'assessNormalization' functions.

#' @param se.obj A summarized experiment object.
#' @param variables Symbol. A symbol and a vector of symbols indicating the columns names of variables in the samples
#' annotation in the SummarizedExperiment object. The 'variables' can be categorical and continuous variables.
#' @param plot.output Logical. Whether to print the plot of all possible assessment plots or not. The default is set to 'TRUE'.
#' @param save.se.obj Logical. Whether to save the results into the SummarizedExperiment object. The default is TRUE.
#'
#' @return A list of all possible assessment metrics for the variables.

#' @importFrom igraph vertex_attr layout_as_tree graph_from_data_frame
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble tibble
#' @export

getAssessmentMetrics <- function(
        se.obj ,
        variables,
        plot.output = TRUE,
        save.se.obj = TRUE) {
    categorical.var <- continuous.var <- NULL

    if(sum(variables %in% colnames(colData(se.obj))) != length(variables))
        stop('All or some "variables" cannot be found in the SummarizedExperiment object.')

    if (!is.null(variables)) {
        var.class <- sapply(
            variables,
            function(x)
                class(colData(se.obj)[[x]]))
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        continuous.var <- names(var.class[var.class %in% c('numeric', 'integer')])
    }
    cat.char <- lapply(
        categorical.var,
        function(x){
            all.char <- gregexpr(pattern ='_', text = x)[[1]]
            all.char[1:length(all.char)]
        })
    categorical.var <- gsub('_', '.', categorical.var)
    names(cat.char) <- categorical.var
    cont.char <- lapply(
        continuous.var,
        function(x){
            all.char <- gregexpr(pattern ='_', text = x)[[1]]
            all.char[1:length(all.char)]
        } )
    continuous.var <- gsub('_', '.', continuous.var)
    names(cont.char) <- continuous.var
    all.var.char <- c(cat.char, cont.char)

    # general assessment #####
    general.rle.med <- data.frame(
        Variables = 'General',
        Metrics = 'RLE',
        Factors = 'General',
        PlotTypes = 'rlePlot',
        Assessments = 'rleMed',
        AssessmentTypes = 'globalLevel'
        )
    general.rle.iqr <- data.frame(
        Variables = 'General',
        Metrics = 'RLE',
        Factors = 'General',
        PlotTypes = 'rlePlot',
        Assessments = 'rleIqr',
        AssessmentTypes = 'globalLevel'
        )
    # possible metrics for continuous variables #####
    if(length(continuous.var) > 0){
        metrics.for.cont.var <- c(
            'rleMedians||scatterPlot||RLE_corrCoeff||globalLevel',
            'rleIqr||scatterPlot||RLE_corrCoeff||globalLevel',
            'pcs||scatterPlot||PCA_averageCorrCoeff||globalLevel',
            'pcs.rseq||lineDotPlot||LRA_averageRseq||globalLevel',
            'corrCoeff||boxPlot||Correlation_corrCoeff||geneLevel',
            'corrCoeff||pvalHist||Correlation_pvalueDis||geneLevel',
            'corrCoeff||pvalHist||Correlation_pvalueNull||geneLevel',
            'corrCoeff||scatterPlot||PartialCorrelation_DA||geneLevel',
            'corrCoeff||barPlot||PartialCorrelation_corrCoeffDiff||geneLevel',
            'corrCoeff||histogram||PartialCorrelation_corrCoeff||geneLevel',
            'corrCoeff||scatterPlot||GeneSetScore_corrCoeff||globalLevel'
            )
        metrics.for.cont.var <- expand.grid(
            continuous.var,
            metrics.for.cont.var)
        colnames(metrics.for.cont.var) <- c('Variables', 'MetricsPlotsAssessment')
        metrics.for.cont.var <- metrics.for.cont.var[order(metrics.for.cont.var$Variables), ]

        # Metrics
        metrics.for.cont.var$Metrics <- unlist(lapply(
            metrics.for.cont.var$MetricsPlotsAssessment,
            function(x){
                d <- strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]
                strsplit(x = d, split = '_')[[1]][1]}))
        # Assessment types
        metrics.for.cont.var$AssessmentTypes <- unlist(lapply(
            metrics.for.cont.var$MetricsPlotsAssessment,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][4]))
        # Plot types
        metrics.for.cont.var$PlotTypes <- unlist(lapply(
            metrics.for.cont.var$MetricsPlotsAssessment,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        # Factors
        metrics.for.cont.var$Factors <- unlist(lapply(
            metrics.for.cont.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        # Assessments
        metrics.for.cont.var$Assessments <- unlist(lapply(
            metrics.for.cont.var$MetricsPlotsAssessment,
            function(x){
                d <- strsplit(x = as.character(x), split = '_')[[1]][2]
                strsplit(x = d, split = '\\|\\|')[[1]][1]}))
        # Final table
        metrics.for.cont.var.table <- metrics.for.cont.var[,c(1,3,4,6,5,7)]
        metrics.for.cont.var.list <- paste0(
            metrics.for.cont.var$Variables,
            '_',
            metrics.for.cont.var$MetricsPlots)

    } else {
        metrics.for.cont.var.table <- NULL
        metrics.for.cont.var.list <- NULL
    }

    # possible metrics for categorical variables #####
    if(length(categorical.var) > 0){
        metrics.for.cat.var <- c(
            'rle||coloredRLEplot||RLE_DA||globalLevel',
            'rleMedians||boxPlot||RLE_rleMed||globalLevel',
            'rleIqr||boxPlot||RLE_rleIqr||globalLevel',
            'pcs||boxPlot||PCA_DA||globalLevel',
            'pcs||scatterPlot||PCA_DA||globalLevel',
            'pcs.corr||lineDotPlot||VCA_averageCorr||globalLevel',
            'ariCoeff||barPlot||ARI_ari||globalLevel',
            'silhouetteCoeff||barPlot||Silhouette_silhouetteCoeff||globalLevel',
            'fStat||boxPlot||ANOVA_fStat||geneLevel',
            'fStat||pvalHist||ANOVA_pvalueDis||geneLevel',
            'fStat||pvalHist||ANOVA_pvalueNull||geneLevel',
            'pValue||pvalHist||DGE_pvalueNull||geneLevel',
            'pValue||pvalHist||DGE_pvalueDis||geneLevel')
        metrics.for.cat.var <- expand.grid(
            categorical.var,
            metrics.for.cat.var)
        colnames(metrics.for.cat.var) <- c('Variables', 'MetricsPlotsAssessments')
        metrics.for.cat.var <- metrics.for.cat.var[order(metrics.for.cat.var$Variables), ]
        # Metrics
        metrics.for.cat.var$Metrics <- unlist(lapply(
            metrics.for.cat.var$MetricsPlotsAssessment,
            function(x){
                d <- strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]
                strsplit(x = d, split = '_')[[1]][1]}))
        # Assessment types
        metrics.for.cat.var$AssessmentTypes <- unlist(lapply(
            metrics.for.cat.var$MetricsPlotsAssessment,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][4]))
        # Plot types
        metrics.for.cat.var$PlotTypes <- unlist(lapply(
            metrics.for.cat.var$MetricsPlotsAssessments,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        # Factors
        metrics.for.cat.var$Factors <- unlist(lapply(
            metrics.for.cat.var$MetricsPlotsAssessments,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        # Assessments
        metrics.for.cat.var$Assessments <- unlist(lapply(
            metrics.for.cat.var$MetricsPlotsAssessment,
            function(x){
                d <- strsplit(x = as.character(x), split = '_')[[1]][2]
                strsplit(d, split = '\\|\\|')[[1]][1]}))

        # Final table
        metrics.for.cat.var.table <- metrics.for.cat.var[,c(1,3,4,6,5,7)]
        metrics.for.cat.var.list <- paste0(
            metrics.for.cat.var$Variables,
            '_',
            metrics.for.cat.var$MetricsPlots)

    } else {
        metrics.for.cat.var.table <- NULL
        metrics.for.cat.var.list <- NULL
        }

    # combined ari and silhouette
    if(length(categorical.var) > 1){
        metrics.for.two.cat.var <- c(
            'ariCoeff||combinedPlot||ARI',
            'silhouetteCoeff||combinedPlot||Silhouette')
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        two.cat.vars <- apply(
            combn(x = categorical.var, m = 2),
            2,
            function(x) paste0(x, collapse = '&'))

        metrics.for.two.cat.var <- expand.grid(
            two.cat.vars,
            metrics.for.two.cat.var)
        colnames(metrics.for.two.cat.var) <- c('Variables', 'MetricsPlots')
        metrics.for.two.cat.var <- metrics.for.two.cat.var[order(metrics.for.two.cat.var$Variables), ]
        metrics.for.two.cat.var$Metrics <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]))
        metrics.for.two.cat.var$PlotTypes <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        metrics.for.two.cat.var$Factors <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        metrics.for.two.cat.var.table <- metrics.for.two.cat.var[,c(1,5,4,3)]
        metrics.for.two.cat.var.table$Assessments <- 'DA'
        metrics.for.two.cat.var.table <- metrics.for.two.cat.var.table[,c(1,4,2,3,5)]
        metrics.for.two.cat.var.table$AssessmentTypes <- 'globalLevel'
        metrics.for.two.cat.var.list <- paste0(
            metrics.for.two.cat.var$Variables,
            '_',
            metrics.for.two.cat.var$MetricsPlots)
    } else {
        metrics.for.two.cat.var.table <- NULL
        metrics.for.two.cat.var.list <- NULL
    }

    # final metric files ####
    final.metrics.list <- c(
        metrics.for.cont.var.list,
        metrics.for.cat.var.list,
        metrics.for.two.cat.var.list,
        'GeneralRLEplot')
    final.metrics.table <- as.data.frame(rbind(
        metrics.for.cont.var.table,
        metrics.for.cat.var.table,
        metrics.for.two.cat.var.table))
    final.metrics.table$Variables <- as.character(final.metrics.table$Variables)
    rownames(final.metrics.table) <- c(1:nrow(final.metrics.table))
    for(i in names(all.var.char)){
        if(all.var.char[[i]][1] != -1){
            cur.name <- names(all.var.char[i])
            for(j in 1:length(all.var.char[[i]]) ){
                pos.to.rep <- all.var.char[[i]][j]
                init.name <- paste0(
                    substring(cur.name, 1,  pos.to.rep - 1),
                    '_',
                    substring(cur.name, pos.to.rep + 1))
                cur.name <- init.name
                }
            final.metrics.table$Variables[final.metrics.table$Variables == names(all.var.char[i])] <- cur.name
        }
    }
    final.metrics.table.toplot <- final.metrics.table
    final.metrics.table <- rbind(
        final.metrics.table ,
        general.rle.med,
        general.rle.iqr
        )
    final.metrics.table$Code <- paste0('PA', 1:nrow(final.metrics.table))

    # plot ####
    steps <- y <- label <- xmin <- xmax <- type <- s_e <- ymin <- ymax <- id <- NULL
    plot.metrics <- lapply(
        variables,
        function(x){
            sub.final.metrics.table.toplot <- final.metrics.table[final.metrics.table$Variables == x, ]
            new.order <- mixedorder(x = sub.final.metrics.table.toplot$Code, decreasing = TRUE)
            sub.final.metrics.table.toplot <- sub.final.metrics.table.toplot[new.order , ]
            metrics.tree <- tibble::tibble(
                from = sub.final.metrics.table.toplot$Variables,
                to = paste(
                    paste0('M: ', sub.final.metrics.table.toplot$Metrics),
                    paste0('V: ', sub.final.metrics.table.toplot$Factors),
                    paste0('P: ', sub.final.metrics.table.toplot$PlotTypes),
                    paste0('A: ', sub.final.metrics.table.toplot$Assessments),
                    paste0('T: ', sub.final.metrics.table.toplot$AssessmentTypes),
                    paste0('C: ', sub.final.metrics.table.toplot$Code),
                    sep = '\n')
            )
            graph.tree <- igraph::graph_from_data_frame(metrics.tree, directed = TRUE)
            coords <- igraph::layout_as_tree(graph = graph.tree)
            colnames(coords) <- c("x", "y")
            output.df <- tibble::as_tibble(coords) %>%
                mutate(
                    steps = igraph::vertex_attr(graph.tree, "name"),
                    label = steps,
                    x = x * -1,
                    type = factor(1)
                )
            plot.nodes <- output.df %>%
                mutate(
                    xmin = x - 0.48,
                    xmax = x + 0.48
                )
            plot.nodes$ymax <- 0.4
            plot.nodes$ymin <- 0.1
            index.even.rows <- seq(from = 2, to = nrow(plot.nodes), 2)
            plot.nodes$ymax[index.even.rows] <-  0.3
            plot.nodes$ymin[index.even.rows] <- 0
            plot.edges <- metrics.tree %>%
                dplyr::mutate(id = dplyr::row_number()) %>%
                pivot_longer(
                    cols = c("from", "to"),
                    names_to = "s_e",
                    values_to = "steps") %>%
                dplyr::left_join(plot.nodes, by = "steps") %>%
                dplyr::select(-c(label, type, y, xmin, xmax)) %>%
                dplyr::mutate(y = ifelse(s_e == "from", ymin, ymax)) %>%
                dplyr::select(-c(ymin, ymax))
            plot.nodes$xmin[1] <- 0
            plot.nodes$ymin[1] <- 0
            plot.nodes$xmax[1] <- 0
            plot.nodes$ymax[1] <- 0
            p <- ggplot() +
                ylim(0,0.5) +
                geom_rect(
                    data = plot.nodes,
                    mapping = aes(
                        xmin = xmin,
                        ymin = ymin,
                        xmax = xmax,
                        ymax = ymax,
                        fill = type),
                    alpha = 0.5,
                    fill = 'grey80'
                )
            x2 <- y2 <- NULL
            plot.nodes$y2 <- c(plot.nodes$ymin + plot.nodes$ymax)/2
            plot.nodes$x2 <- c(plot.nodes$xmin + plot.nodes$xmax)/2
            p <- p + geom_text(
                data = plot.nodes[-1,],
                mapping = aes(x = x2, y = y2, label = label),
                color = "black", size = 3, stat = "identity", hjust = 0.5, vjust = 0.5
            ) +
                theme(
                    panel.background = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    plot.margin = unit(c(0,-.5,0,-.5), 'lines'),
                    legend.position = 'none') +
                geom_rect(
                    mapping = aes(
                        xmin = min(plot.nodes$xmin),
                        xmax = max(plot.nodes$xmax),
                        ymin = min(plot.nodes$y),
                        ymax = 0.5
                    ),
                    alpha = 0.005,
                    fill = 'grey',
                    color = 'black',
                    size = 1,
                ) + geom_text(
                    data = plot.nodes[1,],
                    mapping = aes(x = 0, y = 0.45, label = label),
                    hjust = 0.5,
                    vjust = 0.5,
                    color = "black",
                    size = 6,
                    fontface = 'bold'
                )
            p
        })
    plot.caption <- expression(atop(
        scriptstyle("M: metrics | V: variable | P: plot type | A: assessment | T: type | C: Code" ))
        )
    all.plots <- ggpubr::ggarrange(
        plotlist = plot.metrics,
        common.legend = TRUE,
        ncol = 1
        )
    all.plots <- ggpubr::annotate_figure(
        p = all.plots,
        bottom = ggpubr::text_grob(label = plot.caption, size = 20)
        )
    if(isTRUE(plot.output)) print(all.plots)
    # save plot
    if(isTRUE(save.se.obj)){
        if(!'AssessmentMetrics' %in% names(se.obj@metadata)){
            se.obj@metadata[['AssessmentMetrics']] <- list()
        }
        if(!'plot' %in% names(se.obj@metadata[['AssessmentMetrics']])){
            se.obj@metadata[['AssessmentMetrics']][['plot']] <- list()
        }
        se.obj@metadata[['AssessmentMetrics']][['plot']] <- all.plots
        if(!'metrics.list' %in% names(se.obj@metadata[['AssessmentMetrics']])){
            se.obj@metadata[['AssessmentMetrics']][['metrics.list']] <- list()
        }
        se.obj@metadata[['AssessmentMetrics']][['metrics.list']] <- unlist(final.metrics.list)
        if(!'metrics.table' %in% names(se.obj@metadata[['AssessmentMetrics']])){
            se.obj@metadata[['AssessmentMetrics']][['metrics.table']] <- list()
        }
        se.obj@metadata[['AssessmentMetrics']][['metrics.table']] <- final.metrics.table
        return(se.obj)

    }else {
        return(all.metrics = list(
            final.metrics.list = unlist(final.metrics.list),
            final.metrics.table = final.metrics.table))
    }
}


