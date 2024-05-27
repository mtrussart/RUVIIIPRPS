#' Generates Kaplan-Meier survival plots.

#' @author Ramyar Molania

#' @description
#' This function computes gene-level and variable-level survival analysis.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartile ranges. The default is set to "all, which
#' indicates all the assays of the SummarizedExperiment object will be selected.
#' @param genes Symbol. A symbol or a vector of symbols for the selection of genes for gene level survival analysis. The
#' default is set to 'all'.
#' @param variable Symbol. A symbol that indicates the column name of the SummarizedExperiment object that contains
#' a categorical variable e.g. tumour subtypes, ... for variable-level survival analysis. The default is set to 'NULL'.
#' @param survival.time Symbol.A symbol that indicates the column name of the SummarizedExperiment object that contains
#' survival time.
#' @param survival.events Symbol.A symbol that indicates the column name of the SummarizedExperiment object that contains
#' survival event (0 and 1).
#' @param expr.stratify Numeric. A numeric value that indicates how the gene expression of each genes should be stratified
#' for gene-level survival analysis. The stratification is based on the quantiles of gene expression. The default is set
#' to 2. This means expression level of each specified genes will be divvied in to group based on 50% quantiles.
#' @param return.p.values Logical. Indicates whether to calculate p.values in the survival analysis or not. The default
#' is set to 'TRUE'.
#' @param return.survival.plots Logical. Indicates whether to generate survival plot or not. The default
#' is set to 'FALSE'.
#' @param plot.output Logical. Indicates whether to plot variable-level survival analysis or not. The default is set to
#' 'TRUE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object or
#' to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom survival survfit Surv
#' @importFrom survminer ggsurvplot
#' @export

computeSurvival <- function(
        se.obj ,
        assay.names = 'all' ,
        genes = 'all',
        variable = NULL,
        survival.time,
        survival.events,
        expr.stratify = 2,
        return.p.values = TRUE,
        return.survival.plots = FALSE,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
        ){

    printColoredMessage(message = '------------The computeSurvival function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check inputs ####
    if(!survival.time %in% colnames(colData(se.obj))){
        stop('The "survival.time" cannot be found in the SummarizedExperiment object.')
    }
    if(!survival.events %in% colnames(colData(se.obj))){
        stop('The "survival.events" cannot be found in the SummarizedExperiment object.')
    }
    if(!is.null(variable)){
        if(!variable %in% colnames(colData(se.obj))){
            stop('The "variable" cannot be found in the SummarizedExperiment object.')
        }
    }
    if(!is.null(genes)){
        if (is.logical(expr.stratify)){
            stop('The "expr.stratify" must be a numeric value.')
        } else if(is.null(expr.stratify)){
            stop('The "expr.stratify" cannot be NULL.')
        } else if (expr.stratify == 1){
            stop('The "expr.stratify" cannot be set to 1.')
        }
    }
    if(isFALSE(return.p.values) & isFALSE(return.survival.plots)){
        stop('Both "return.p.values" and "return.survival.plots" cannot be set to "FALSE".')
    }

    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # Select colored ####
    selected.colores <-  c(
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

    # Gene level survival analysis ####
    if(!is.null(genes)){
        ## Assess genes ####
        if (length(genes) == 1 && genes == 'all')
            genes <- row.names(se.obj)
        if(!sum(genes %in% row.names(se.obj)) == length(genes)){
            stop('All or some of the  "genes" cannot be found in the SummarizedExperiment object.')
        }
        # Survival data ####
        survival.data <- as.data.frame(colData(se.obj))[ , c(survival.time, survival.events)]
        colnames(survival.data) <- c('time', 'status')

        # Gene grouping ####
        all.genes.grouping <- lapply(
            levels(assay.names),
            function(i){
                all.genes.grouping <- lapply(
                    genes,
                    function(x){
                        # Group gene expression
                        quantiles <- stats::quantile(
                            x = assay(x = se.obj, i = i)[x , ],
                            probs = seq(0, 1, c(1/expr.stratify)))
                        expr.lables <- round(x = seq(0, 1, c(1/expr.stratify)) , digits = 2)
                        expr.lables <- sapply(
                            1:c(length(expr.lables)-1),
                            function(x) paste(expr.lables[x], expr.lables[x+1], sep = ':'))
                        genes.group <- NULL
                        genes.group <- as.numeric(
                            cut(x = assay(x = se.obj, i = i)[x , ],
                                breaks = quantiles,
                                include.lowest = TRUE))
                        temp.data$gene <- paste0('group', expr.lables[genes.group])
                        return(temp.data)
                    })
                names(all.genes.grouping) <- genes
                all.genes.grouping
            })
        names(all.genes.grouping) <- levels(assay.names)

        # p.values ####
        if(isTRUE(return.p.values)){
            all.genes.p.values <- lapply(
                levels(assay.names),
                function(i){
                    all.genes.p.values <- lapply(
                        genes,
                        function(x){
                            diff.surv <- survival::survdiff(formula = survival::Surv(time , status) ~ gene, data = all.genes.grouping[[i]][[x]])
                            diff.surv$pvalue
                        })
                    names(all.genes.p.values) <- genes
                    all.genes.p.values <- t(as.data.frame(all.genes.p.values))
                    colnames(all.genes.p.values) <- 'p.value'
                    all.genes.p.values
                })
            names(all.genes.p.values) <- levels(assay.names)
        }
        # Survival plots ####
        if(isTRUE(return.survival.plots)){
            all.genes.survial.plots <- lapply(
                levels(assay.names),
                function(i){
                    all.genes.survial.plots <- lapply(
                        genes,
                        function(x){
                            fit.surv <- survival::survfit(formula = survival::Surv(time , status) ~ gene, data = all.genes.grouping[[i]][[x]])
                            p <- ggsurvplot(fit.surv,
                                            pval = TRUE,
                                            conf.int = TRUE,
                                            risk.table = FALSE,
                                            risk.table.col = "strata",
                                            linetype = 1,
                                            legend = "bottom",
                                            title = x,
                                            palette = selected.colores[seq(length(unique(genes.group)))])
                            p

                        })
                })
            names(all.genes.survial.plots) <- levels(assay.names)
        }
    }
    # Variable survival analysis ####
    if(!is.null(variable)){
        # Save the data
        ## add results to the SummarizedExperiment object ####
        printColoredMessage(
            message = '-- Save the correlation coefficients data:',
            color = 'magenta',
            verbose = verbose)
        # Survival data ####
        survival.data <- as.data.frame(colData(se.obj))[ , c(survival.time, survival.events, variable)]
        colnames(survival.data) <- c('time', 'status', 'variable')
        diff.surv <- survival::survdiff(formula = survival::Surv(time , status) ~ variable, data = survival.data)
        surv.pvalue.variable <- diff.surv$pvalue
        fit.surv <- survival::survfit(formula = survival::Surv(time , status) ~ variable, data = survival.data)
        surv.plot.variable <- ggsurvplot(fit.surv,
                        pval = TRUE,
                        conf.int = TRUE,
                        risk.table = FALSE,
                        risk.table.col = "strata",
                        linetype = 1,
                        legend = "bottom",
                        title = x,
                        palette = selected.colores[seq(length(unique(survival.data$variable)))])
        if(isTRUE(plot.output))
            print(surv.plot.variable)
    }

    # Save the data ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save the survival p.values and plots data:',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Save all the survival p.values and plots to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ### for all assays
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'Survival' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['Survival']] <- list()
            }
            if(!is.null(genes)){
                ## check if metadata metric already exist for this assay and this metric
                if (!'gene.level' %in% names(se.obj@metadata[['metric']][[x]][['Survival']])) {
                    se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']] <- list()
                }
                if(isTRUE(return.p.values)){
                    ## check if metadata metric already exist for this assay and this metric
                    if (!'p.values' %in% names(se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']])) {
                        se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']][['p.values']] <- list()
                    }
                    se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']][['p.values']] <- all.genes.p.values[[x]]
                }
                if(isTRUE(return.survival.plots)){
                    ## check if metadata metric already exist for this assay and this metric
                    if (!'plots' %in% names(se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']])) {
                        se.obj@metadata[['metric']][[x]][['Survival']][['gene.level']][['plots']] <- all.genes.survial.plots[[x]]
                    }
                }
            }
            if(!is.null(variable)){
                if (!'variable.level' %in% names(se.obj@metadata[['metric']][[x]][['Survival']])) {
                    se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']] <- list()
                }
                if(isTRUE(return.p.values)){
                    if (!'p.values' %in% names(se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']])) {
                        se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']][['p.values']] <- list()
                    }
                    se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']][['p.values']] <- surv.pvalue.variable
                }
                if(isTRUE(return.survival.plots)){
                    ## check if metadata metric already exist for this assay and this metric
                    if (!'plots' %in% names(se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']])) {
                        se.obj@metadata[['metric']][[x]][['Survival']][['variable.level']][['plots']] <- surv.plot.variable
                    }
                }
            }
        }
        printColoredMessage(
            message = paste0('The RLE data of individual assay is saved to the metadata@metric in SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The computeSurvival function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(corr.coeff)
    }

}


