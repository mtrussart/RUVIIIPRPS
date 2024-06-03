#' Compute ANOVA between individual gene expression and a categorical variable.

#' @author Ramyar Molania

#' @description
#' This function calculates the ANOVA between individual gene expression of the assays in a SummarizedExperiment object
#' and a categorical variable as factor.

#' @details
#' ANOVA enables us to assess the effects of a given qualitative variable (which we call a factor) on gene expression
#' measurements across any set of groups (labeled by the levels of the factor) under study. We use ANOVA F-statistics
#' to summarize the effects of a qualitative source of unwanted variation (for example, batches) on the expression levels
#' of individual genes, where genes having large F-statistics are deemed to be affected by the unwanted variation.
#' We also use ANOVA tests (the aov() function in R) to assign P values to the association between tumor purity and
#' molecular subtypes.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the ANOVA. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param method A character string indicating which method
#' is to be used for the differential analysis: "aov" or "welch". By default "aov" will
#' be selected.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of the number of genes
#' from the top listing of anova by default it is set to FALSE.
#' @param nb.top.genes Defines the number of genes from the top or bottom listing of anova to plot,
#' by default is set to 3.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object
#' 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return SummarizedExperiment A SummarizedExperiment object containing the log2 F-statistics of ANOVA on the continuous
#'  variable and if requested the associated boxplot.

#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer %>%
#' @importFrom dplyr mutate
#' @import ggplot2

#' @export

computeGenesVariableAnova <- function(
        se.obj,
        assay.names = 'all',
        variable,
        method = 'aov',
        apply.log = TRUE,
        pseudo.count = 1,
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        apply.round = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computeGenesVariableAnova function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be null.')
    }
    if(length(assay.names) == 1 && assay.names!= 'all'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(sum(!assay.names %in% names(assays(se.obj))) > 0 )
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(variable)) {
        stop('The variable cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, ' should be a categorical variable.'))
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one level. ANOVA cannot be performed.'))
    }
    if (!method %in% c('aov', 'welch')) {
        stop('The method should be one of the "aov" or "welch".')
    }
    if (apply.log){
        if(pseudo.count < 0 )
            stop('The value of pseudo.count cannot be negative.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # data transformation ####
    printColoredMessage(
        message = '-- Data log transformation:',
        color = 'magenta',
        verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (apply.log & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('- Apply log2 + ', pseudo.count,  ' (pseudo.count) transformation on the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (apply.log & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- Apply log2 transformation on the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x))
            } else {
                printColoredMessage(
                    message = paste0('- The "',x, '" data will be used without log transformation.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = 'Please note, the assay should be in log scale before computing ANOVA.',
                    color = 'red',
                    verbose = verbose)
                expr <- assay(x = se.obj, i = x)
            }
        })
    names(all.assays) <- levels(assay.names)

    # anova ####
    printColoredMessage(
        message = paste0('Perform ANOVA between individual genes expression of the assay(s) and the "',
                         variable, '" variable.'),
        color = 'magenta',
        verbose = verbose)
    all.aov <- lapply(
        levels(assay.names),
        function(x) {
            # anova ####
            if (method == 'aov') {
                printColoredMessage(
                    message = paste0('- Perform the ANOVA with equal variance for "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                anova.genes.var <- row_oneway_equalvar(
                    x = all.assays[[x]],
                    g = se.obj@colData[, variable]
                    )
            } else if (method == 'welch') {
                printColoredMessage(
                    message = paste0(
                        '- Perform ANOVA with Welch correction for the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose)
                anova.genes.var <- row_oneway_welch(
                    x =  all.assays[[x]],
                    g = se.obj@colData[, variable]
                    )
            }
            row.names(anova.genes.var) <- row.names(se.obj)

            # round the anova statistic obtained to 2 digits ####
            if (isTRUE(apply.round)) {
                anova.genes.var <- cbind(
                    round(anova.genes.var[, 1:9], digits = 3),
                    anova.genes.var[, 10, drop = FALSE])
            }
            # plot top genes ####
            if (plot.top.genes) {
                temp.anova <- anova.genes.var[order(anova.genes.var[, 'statistic'],decreasing = TRUE, na.last = TRUE) , ]
                var <- NULL
                p.high <- as.data.frame(t(all.assays[[x]][row.names(temp.anova)[c(1:nb.top.genes)],]))
                p.high <- mutate(p.high , var = se.obj@colData[, variable])
                p.high <- pivot_longer(
                    data = p.high,
                    cols =  -var,
                    names_to = 'genes',
                    values_to = 'expr')
                p.high <- ggplot(p.high, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    facet_wrap( ~ genes) +
                    ggtitle(paste0(nb.top.genes," Top affected genes by the variable ", variable, " for ", x)) +
                    theme(panel.background = element_blank(),
                          axis.line = element_line(colour = 'black', linewidth = 1),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(
                              size = 10,
                              angle = 45,
                              hjust = 1,
                              vjust = 1),
                          axis.text.y = element_text(size = 12),
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 14),
                          strip.text.x = element_text(size = 10),
                          plot.title = element_text(size = 12))
                print(p.high)
                rm(temp.anova)
            }
            results <- NULL
            results <- list(anova.genes.var =  anova.genes.var)
            return(results)
        })
    names(all.aov) <- levels(assay.names)
    # save the results ####
    printColoredMessage(
        message = '-- Save the ANOVA results :',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- The ANOVA results for the indiviaul assay(s) are saved to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'ANOVA' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['ANOVA']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!method %in% names(se.obj@metadata[['metric']][[x]][['ANOVA']])) {
                se.obj@metadata[['metric']][[x]][['ANOVA']][[method]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['ANOVA']][[method]][[variable]][['F.values']] <-
                all.aov[[x]][['anova.genes.var']][, c('pvalue', 'statistic')]
        }
        printColoredMessage(
            message = paste0('- The ANOVA results of all the assays are saved to the ',
                             ' "se.obj@metadata$metric$AssayName$ANOVA$', method, '$library.size" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    ## return the results as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- The ANOVA results for indiviaul assay are saved as list.',
            color = 'blue',
            verbose = verbose)
        all.aov <- lapply(
            levels(assay.names),
            function(x) log2(all.aov[[x]][['anova.genes.var']][, 'statistic']))
        names(all.aov) <- levels(assay.names)
        printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        return(genes.var.anova = all.aov)
    }
}
