#' Computes variable and gene-level Kaplan-Meier survival analysis. Finalized

#' @author Ramyar Molania

#' @description
#' This function computes gene-level and variable-level survival analysis.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character or a vector of characters specifying the name(s) of the assay(s) in the
#' SummarizedExperiment object for gene-level survival analysis. The default is set to "all", which indicates that all
#' assays in the SummarizedExperiment object will be selected.
#' @param genes Character. A character or a vector of gene names/IDs specifying the genes to be used for gene-level survival
#' analysis. The default is set to 'all'.
#' @param variable Character. A character specifying the column name of the SummarizedExperiment object that contains
#' a categorical variable (e.g., tumor subtypes) for variable-level survival analysis. The default is set to 'NULL'.
#' @param survival.time Character. A character specifying the column name of the SummarizedExperiment object that contains
#' survival time.
#' @param survival.events Character. A character specifying the column name of the SummarizedExperiment object that contains
#' the survival event (0 and 1).
#' @param expr.stratify Numeric. A numeric value indicating how the gene expression of each gene should be stratified
#' for gene-level survival analysis. Stratification is based on the quantiles of gene expression. The default is 2,
#' meaning the expression level of each specified gene will be divided into two groups based on the 50% quantile (median).
#' @param return.p.value Logical. Indicates whether to calculate p-values in the survival analysis. The default is 'TRUE'.
#' @param return.survival.plot Logical. Indicates whether to generate survival plots. The default is set to 'FALSE'.
#' @param plot.output Logical. Indicates whether to plot variable-level survival analysis. The default is 'TRUE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. If 'TRUE', the 'checkSeObj'
#' function will be applied. The default is 'TRUE'.
#' @param remove.na Character. Indicates whether to remove NA or missing values from the assays. Options are 'assays'
#' , 'sample annotation' or 'none'. The default is set to 'assays', which removes all NA or missing values from the assays
#' before performing the gene level analysis. See the 'checkSeObj' function for more details.
#' @param save.se.obj Logical. Indicates whether to save the results of the survival analysis in the metadata of the
#' SummarizedExperiment object or output the results as a list. The default is 'TRUE'.
#' @param verbose Logical. If 'TRUE', displays messages for each step of the function.

#' @importFrom survival survfit survdiff Surv
#' @importFrom survminer ggsurvplot
#' @import RColorBrewer
#' @export

computeSurvival <- function(
        se.obj ,
        assay.names = 'all' ,
        genes = 'all',
        variable = NULL,
        survival.time,
        survival.events,
        expr.stratify = 2,
        return.p.value = TRUE,
        return.survival.plot = FALSE,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
        ){

    printColoredMessage(message = '------------The computeSurvival function starts:',
                        color = 'white',
                        verbose = verbose)

    # Checking inputs ####
    if (!survival.time %in% colnames(colData(se.obj))){
        stop('The "survival.time" cannot be found in the SummarizedExperiment object.')
    }
    if (!survival.events %in% colnames(colData(se.obj))){
        stop('The "survival.events" cannot be found in the SummarizedExperiment object.')
    }
    if (!is.null(variable)){
        if(!variable %in% colnames(colData(se.obj))){
            stop('The "variable" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!is.null(genes)){
        if (is.logical(expr.stratify)){
            stop('The "expr.stratify" must be a numeric value.')
        } else if(is.null(expr.stratify)){
            stop('The "expr.stratify" cannot be NULL.')
        } else if (expr.stratify == 1){
            stop('The "expr.stratify" cannot be set to 1.')
        }
    }
    if (isFALSE(return.p.value) & isFALSE(return.survival.plot)){
        stop('Both "return.p.value" and "return.survival.plot" cannot be set to "FALSE".')
    }
    if (!is.logical(return.p.value)){
        stop('The "return.p.value" must be logical(TRUE or FALSE).')
    }
    if (!is.logical(return.survival.plot)){
        stop('The "return.survival.plot" must be logical(TRUE or FALSE).')
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

    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose
            )
    }

    # Selecting colores ####
    selected.colores <-  c(
        c("#E7B800", "#2E9FDF", 'red4'),
        RColorBrewer::brewer.pal(8, "Dark2")[-5],
        RColorBrewer::brewer.pal(10, "Paired")
        )

    # Gene level survival analysis ####
    if (!is.null(genes)){
        ## assess provided genes ####
        if (length(genes) == 1 && genes == 'all')
            genes <- row.names(se.obj)
        if (!sum(genes %in% row.names(se.obj)) == length(genes)){
            stop('All or some of the  "genes" cannot be found in the SummarizedExperiment object.')
        }
        survival.data <- as.data.frame(colData(se.obj))[ , c(survival.time, survival.events)]
        colnames(survival.data) <- c('time', 'status')

        # grouping genes based on their expression ####
        all.genes.grouping <- lapply(
            levels(assay.names),
            function(i){
                expr.data <- assay(x = se.obj, i = i)
                genes.grouping <- lapply(
                    genes,
                    function(x){
                        # Group gene expression
                        quantiles <- stats::quantile(
                            x = expr.data[x , ],
                            probs = seq(0, 1, c(1/expr.stratify))
                            )
                        expr.lables <- round(x = seq(0, 1, c(1/expr.stratify)) , digits = 2)
                        expr.lables <- sapply(
                            1:c(length(expr.lables)-1),
                            function(x) paste(expr.lables[x], expr.lables[x+1], sep = ':')
                            )
                        genes.group <- as.numeric(
                            cut(x = assay(x = se.obj, i = i)[x , ],
                                breaks = quantiles,
                                include.lowest = TRUE)
                            )
                        survival.data$gene <- paste0('group', expr.lables[genes.group])
                        return(survival.data)
                    })
                names(genes.grouping) <- genes
                genes.grouping
            })
        names(all.genes.grouping) <- levels(assay.names)

        # obtaining p.values ####
        if (isTRUE(return.p.value)){
            all.genes.p.values <- lapply(
                levels(assay.names),
                function(i){
                    survival.data <- all.genes.grouping[[i]]
                    all.genes.p.values <- lapply(
                        genes,
                        function(x){
                            diff.surv <- survival::survdiff(
                                formula = survival::Surv(time , status) ~ gene,
                                data = survival.data[[x]]
                                )
                            diff.surv$pvalue
                        })
                    names(all.genes.p.values) <- genes
                    all.genes.p.values <- t(as.data.frame(all.genes.p.values))
                    colnames(all.genes.p.values) <- 'p.value'
                    all.genes.p.values
                })
            names(all.genes.p.values) <- levels(assay.names)
        }
        # Obtaining survival plots ####
        if (isTRUE(return.survival.plot)){
            all.genes.survial.plots <- lapply(
                levels(assay.names),
                function(i){
                    survival.data <- all.genes.grouping[[i]]
                    all.genes.survial.plots <- lapply(
                        genes,
                        function(x){
                            fit.surv <- survival::survfit(
                                formula = survival::Surv(time , status) ~ gene,
                                data = survival.data[[x]]
                                )
                            survival.plot <- ggsurvplot(
                                fit.surv,
                                pval = TRUE,
                                conf.int = FALSE,
                                risk.table = FALSE,
                                risk.table.col = "strata",
                                linetype = 1,
                                legend = "bottom",
                                title = x,
                                palette = selected.colores[seq(length(unique(survival.data[[x]][['gene']])))]
                                )
                            return(survival.plot)

                        })
                })
            names(all.genes.survial.plots) <- levels(assay.names)
        }
    }
    # Variable survival analysis ####
    if (!is.null(variable)){
        printColoredMessage(
            message = '-- Applying survival analysis on the provided variable:',
            color = 'magenta',
            verbose = verbose
            )
        # Survival data ####
        survival.data <- as.data.frame(colData(se.obj))[ , c(survival.time, survival.events, variable)]
        colnames(survival.data) <- c('time', 'status', 'group')
        survdiff.fun <- function(survival.data) {
            result <- survdiff(Surv(time, status) ~ group, data = survival.data)
            return(result)
        }
        diff.surv <- survdiff.fun(survival.data)
        surv.pvalue.variable <- diff.surv$pvalue

        survfit.fun <- function(survival.data) {
            result <- survfit(Surv(time, status) ~ group, data = survival.data)
            return(result)
        }
        fit.surv <- survfit.fun(survival.data)
        ggsurvplot.fun <- function(fit.surv) {
            p.plot <- ggsurvplot(fit = fit.surv,
                       data = survival.data,
                       pval = TRUE,
                       conf.int = FALSE,
                       risk.table = FALSE,
                       risk.table.col = "strata",
                       linetype = 1,
                       legend = "bottom",
                       title = variable,
                       palette = selected.colores[seq(length(unique(survival.data$group)))])
            ggpubr::ggpar(p.plot,
                  font.main = c(12),
                  font.x = c(12),
                  font.y = c(12),
                  font.caption = c(20),
                  font.legend = c(12),
                  font.tickslab = c(12))
        }
        surv.plot.variable <- ggsurvplot.fun(fit.surv)
        if (isTRUE(plot.output))
            print(surv.plot.variable)
    }

    # Save the data ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Saving the p-values and plots of the survival analysis:',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the survival p.values and plots to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ### for all assays
        if (!is.null(genes)){
            if (isTRUE(return.p.value)){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = assay.names,
                    assessment.type = 'gene.level',
                    assessment = 'Survival',
                    method = 'Kaplan.Meier',
                    file.name = 'data',
                    variables = 'p.values',
                    results.data = all.genes.p.values
                )
            }
           if (isTRUE(return.survival.plot)){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = assay.names,
                    assessment.type = 'gene.level',
                    assessment = 'Survival',
                    method = 'Kaplan.Meier',
                    file.name = 'data',
                    variables = 'plots',
                    results.data = all.genes.p.values
                )
            }

        }
        if (!is.null(variable)){
            if(isTRUE(return.survival.plot)){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = assay.names,
                    assessment.type = 'gene.level',
                    assessment = 'Survival',
                    method = 'Kaplan.Meier',
                    file.name = 'data',
                    variables = 'p.values',
                    results.data = surv.plot.variable
                )
            }
            if (isTRUE(return.survival.plot)){
                se.obj <- addMetricToSeObj(
                    se.obj = se.obj,
                    slot = 'Metrics',
                    assay.names = assay.names,
                    assessment.type = 'gene.level',
                    assessment = 'Survival',
                    method = 'Kaplan.Meier',
                    file.name = 'data',
                    variables = 'plots',
                    results.data = surv.plot.variable
                )
            }

        }
        printColoredMessage(
            message = paste0('All the survival results of individual assays is saved to the in SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The computeSurvival function finished.',
                            color = 'white',
                            verbose = verbose
                            )
        return(se.obj = se.obj)
    }
    if (isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The computeSurvival function finished.',
                            color = 'white',
                            verbose = verbose
        )
        if (is.null(genes)){
            if (isTRUE(return.survival.plot) & isTRUE(return.p.value)){
                return(list(
                    all.genes.survial.plots = all.genes.survial.plots,
                    all.genes.p.values = all.genes.p.values)
                    )
            }
            if (isTRUE(return.survival.plot) & isFALSE(return.p.value)){
                return(list(all.genes.survial.plots = all.genes.survial.plots))
            }
            if (isFALSE(return.survival.plot) & isTRUE(return.p.value)){
                return(list(all.genes.p.values = all.genes.p.values))
            }
        }
        if (!is.null(variable)){
            return(surv.plot.variable)
        }
    }
}


