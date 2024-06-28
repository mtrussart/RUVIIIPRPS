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
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment object to compute the correlation. By default all the assays of the SummarizedExperiment class
#' object will be selected.
#' @param variable Symbol. A symbol indicating a column name of the SummarizedExperiment object that contains a
#' categorical variable such as experimental batches, ... .
#' @param method Symbol. A symbol indicating which method is to be used for the ANOVA analysis. The options are "aov" and
#' "welch". The default is set to "aov".
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of "nb.top.genes" number genes with highest
#' and lowest F-statistics across the specified variable. The default is set to 'FALSE'.
#' @param nb.top.genes Numeric. A numeric value Defining the number of genes from the top or bottom listing of ANOVA
#' F-statistics to plot. The default is set to 3.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before applying ANOVA. The default
#' is set to 'TRUE'.
#' @param pseudo.count Numeric. A numeric value as a pseudo count to be added to all measurements before applying log transformation.
#' The default is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object. The default is set to
#' 'TRUE'.
#' @param remove.na Symbol. A symbol that specifies whether to eliminate missing values from either 'assays', 'sample.annotation',
# 'both', or 'none'. When 'assays' is chosen, genes containing missing values will be omitted. If 'sample.annotation'
# is selected, samples with NA or missing values for each 'variables' will be excluded. The default is 'both'.
#' @param apply.round Logical. Indicates whether to round the values of F-statistics. The default is set to 'TRUE'.
#' @param save.se.obj Logical. Indicates whether to save the results, ANOVA F-statistics and p-values, in the metadata
#' of the SummarizedExperiment object or to output the result as a list or vector. The default is set to 'TRUE'.
#' @param override.check Logical. When set to 'TRUE', the function verifies the current SummarizedExperiment object to
#' determine if the ANOVA has already been computed for the current parameters. If it has, the metric will not be recalculated.
#' The default is set to FALSE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object containing the log2 F-statistics of ANOVA on the continuous
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
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeGenesVariableAnova function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if(isTRUE(override.check)){
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'gene.level',
            assessment = 'ANOVA',
            method = method,
            variable = variable,
            file.name = 'fstatistics.pvalues'
        )
        if(is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if(isTRUE(compute.metric)){
        # Check the inputs ####
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
        if (isTRUE(apply.log)){
            if(pseudo.count < 0 )
                stop('The value of pseudo.count cannot be negative.')
        }

        # Check the assays ####
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

        # Data transformation ####
        printColoredMessage(
            message = '-- Data log transformation:',
            color = 'magenta',
            verbose = verbose)
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assessment = 'ANOVA',
            verbose = verbose
        )

        # Apply ANOVA ####
        printColoredMessage(
            message = paste0('-- Perform ANOVA between individual genes expression of the assay(s) and the "',
                             variable, '" variable.'),
            color = 'magenta',
            verbose = verbose
        )
        all.aov <- lapply(
            levels(assay.names),
            function(x) {
                # compute anova ####
                if (method == 'aov') {
                    printColoredMessage(
                        message = paste0('- Perform the ANOVA with equal variance for "', x, '" data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    anova.genes.var <- row_oneway_equalvar(
                        x = all.assays[[x]],
                        g = se.obj@colData[[variable]]
                    )
                } else if (method == 'welch') {
                    printColoredMessage(
                        message = paste0(
                            '- Perform ANOVA with Welch correction for the "', x, '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    anova.genes.var <- row_oneway_welch(
                        x =  all.assays[[x]],
                        g = se.obj@colData[[variable]]
                    )
                }
                row.names(anova.genes.var) <- row.names(se.obj)

                # round the anova statistic obtained to 2 digits ####
                if (isTRUE(apply.round)) {
                    anova.genes.var <- cbind(
                        round(anova.genes.var[, 1:9], digits = 3),
                        anova.genes.var[, 10, drop = FALSE])
                }
                anova.genes.var <- anova.genes.var[ , c('pvalue', 'statistic')]
                # plot top genes ####
                if (isTRUE(plot.top.genes)) {
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
                return(anova.genes.var)
            })
        names(all.aov) <- levels(assay.names)

        # Save the results ####
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
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'ANOVA',
                method = method,
                variables = variable,
                file.name = 'fstatistics.pvalues',
                results.data = all.aov
            )
            printColoredMessage(
                message = paste0('- The ANOVA results of all the assays are saved to the ',
                                 ' "se.obj@metadata$metric$AssayName$ANOVA$', method, '$library.size" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
            printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                                color = 'white',
                                verbose = verbose
            )
            return(se.obj = se.obj)
        }
        ## return the results as a list ####
        if (isFALSE(save.se.obj)) {
            printColoredMessage(
                message = '- The ANOVA results for indiviaul assay are saved as list.',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                                color = 'white',
                                verbose = verbose
                                )
            return(genes.var.anova = all.aov)
        }
    } else {
        printColoredMessage(message = '------------The computeGenesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose
        )
        return(se.obj = se.obj)
    }
}
