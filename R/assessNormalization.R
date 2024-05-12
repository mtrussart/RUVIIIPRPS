#' Assess the performance of normalization methods.

#' @author Ramyar Molania

#' @description
#' This function summarizes a range of global and gene level metrics obtained using the 'assessVariation' function for
#' individual biological and unwanted variables. The functions returns numerical assessments in a table to assess the
#' performance of difference normalization. Refer to details for more information.

#' @details
#' It is essential to assess the performance of a normalization as a remover of unwanted variation and a preserver of
#' biological variation in the data. An ideal normalization should preserve all known sources of biological variation
#' and reduce or remove the impact of all known sources of unwanted variation. To rank methods according to their
#' performance, we aggregate all evaluation metrics into a overall score as follows.
#' Several assessment will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#' For each continous variable:
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Correlation between gene expression and continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined Silhouette plot of the combined pair of all categorical variables.
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Boxplot of the correlation between gene expression and continuous variable.
#' - It will also output the RLE plot distribution.

#' @references
#' **Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#' **Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023


#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is set to 'all', so all the assays in the SummarizedExperiment object will
#' be selected.
#' @param bio.variables Symbols. Indicating the column names that contain known biological variable(s) in the
#' SummarizedExperiment object. These biological variables can be categorical or continuous.
#' @param uv.variables Symbols. Indicating the column names that contain unwanted variable(s) in the SummarizedExperiment
#' object. These unwanted variables can be categorical or continuous.

#' @param metrics.to.exclude Symbols. A symbol or vector of symbols indication which metrics to be excluded.
#' The metrics can be obtained using the  'getAssessmentMetrics' function. The default is set to NULL. Refer
#' to the details for more information.
#' @param assessments.to.exclude Symbols. A symbol or vector of symbols indicating which assessment metrics to exclude.
#' These metrics can be acquired using the 'getAssessmentMetrics' function. Default is set to NULL. See the details for
#' further information.
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param sil.dist.measure A character string indicating which method
#' is to be used for the differential analysis: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' By default 'euclidean' will be selected.
#' @param ari.clustering.method Symbol. Indicates which clustering methods should be applied on the PCs calculate the ARI.
#' The function provides the 'mclust' or 'hclust' methods. The default is 'hclust'.
#' @param ari.hclust.method Symbol. Indicate the agglomeration method to be used for the 'hclust' method. This should be
#' one of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or
#' 'centroid' (= UPGMC). See the hclust function for more details.
#' @param ari.hclust.dist.measure Symbol. Indicates the distance measure to be used in the dist function. This must be
#' one of euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'. See the dist function for more details.
#' @param fstat.cutoff Numeric. Indicates the number of PCs that should be used to compute the ARI.
#' @param corr.coef.cutoff Numeric. Indicates the number of PCs that should be used to compute the ARI.
#' @param corr.coef.cutoff Numeric. Indicates the number of PCs that should be used to compute the ARI.
#' @param corr.method Symbol. Indicates which correlation methods should be used. Options are The default is 'spearman'.
#' @param anova.method Symbol. Indicates which ANOVA methods should be used. Options are ... . The default is 'aov'.
#' @param plot.output BB
#' @param output.file Path and name of the output file to save the assessments plots in a pdf format.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.
#' @param fstat.cutoff Numeric. Specifies the threshold for selecting genes that have a F-statistics lower than this value
#'  with categorical variables. The default value is 1.
#' @param corr.coef.cutoff Numeric. Specifies the threshold for selecting genes that have an absolute correlation coefficient
#'  lower than this value with continuous variables. The default value is 0.2.

#' @return A SummarizedExperiment object containing all the assessments plots and metrics. If specified it will generate
#' a pdf containing the assessments plots and metrics used for the assessment.

#' @importFrom ggh4x strip_nested elem_list_text elem_list_rect facet_nested
#' @importFrom dplyr summarise_at case_match row_number
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom SummarizedExperiment assays colData
#' @importFrom gridExtra grid.arrange grid.table
#' @importFrom ggforestplot geom_stripes
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot.new text
#' @importFrom stats ks.test IQR
#' @importFrom qvalue qvalue
#' @export


assessNormalization <- function(
        se.obj,
        assay.names = 'all',
        bio.variables,
        uv.variables,
        metrics.to.exclude = NULL,
        assessments.to.exclude = NULL,
        fast.pca = TRUE,
        sil.dist.measure = 'euclidian',
        ari.clustering.method = "hclust",
        ari.hclust.method = "complete",
        ari.hclust.dist.measure = "euclidian",
        fstat.cutoff = 1,
        corr.coef.cutoff = 0.2,
        corr.method = 'spearman',
        anova.method = 'aov',
        output.file = NULL,
        plot.output = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The assessVariation function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs of function ####
    if (length(assay.names) == 1 && assay.names != 'all') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (length(setdiff(assay.names, names(assays(se.obj)))) > 0)
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if(is.null(bio.variables) | is.null(uv.variables)){
        stop('To performe "overall.performance" both "uv.variables" and "bio.variables" must be provided.')
    }


    # Assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # All possible metrics for each variable #####
    printColoredMessage(
        message = 'Find all possible assessment metrics:',
        color = 'magenta',
        verbose = verbose)
    se.obj <- getAssessmentMetrics(
        se.obj = se.obj,
        variables = c(uv.variables, bio.variables),
        plot.output = FALSE,
        save.se.obj = TRUE)
    # Metrics and plots to generate #####
    metrics.table <- se.obj@metadata$AssessmentMetrics$metrics.table
    if(!is.null(metrics.to.exclude)){
        metrics.table <- metrics.table[!metrics.table$Metrics %in% metrics.to.exclude, ]
    }
    if(!is.null(assessments.to.exclude)){
        assessments.table <- metrics.table[!metrics.table$Code %in% assessments.to.exclude, ]
    } else assessments.table <- metrics.table
    printColoredMessage(
        message = paste0(
            nrow(se.obj@metadata$AssessmentMetrics$metrics.table),
            ' assessment will be generated.'),
        color = 'blue',
        verbose = verbose)

    ## rle scores #####
    if('rleMed||rleIqr' %in% assessments.table$Assessments){
        rle.meds <- unlist(lapply(
            assay.names,
            function(x){
                if(is.null(se.obj@metadata$metric[[x]]$RLE$rle.data$rle.med))
                   return(x)
            }))
        if(!is.null(rle.meds)){
            stop(paste0('The RLE medians for the ',
                        paste0(rle.meds, collapse = ' & '),
                        ' data cannot be found in the SummarizedExperiment object.',
                        'Please run the assessVariation() function.'))
        }

        rle.iqrs <- unlist(lapply(
            assay.names,
            function(x){
                if(is.null(se.obj@metadata$metric[[x]]$RLE$rle.data$rle.iqr))
                    return(x)
            }))
        if(!is.null(rle.iqrs)){
            stop(paste0('The RLE IQRs for the ',
                        paste0(rle.meds, collapse = ' & '),
                        ' data cannot be found in the SummarizedExperiment object.',
                        'Please run the assessVariation() function.'))
        }

        # rle medians
        all.rle.meds <- sapply(
            levels(assay.names),
            function(x) abs(se.obj@metadata$metric[[x]]$RLE$rle.data$rle.med))
        max.rle.meds <- max(all.rle.meds) * ncol(se.obj)
        rle.med.scores <- sapply(
            assay.names,
            function(x){
                med <- abs(se.obj@metadata$metric[[x]]$RLE$rle.data$rle.med)
                1 - sum(med)/max.rle.meds
            })
        names(rle.med.scores) <- levels(assay.names)
        # rle iqr
        rle.iqr.scores <- sapply(
            levels(assay.names),
            function(x){
                iqr <- se.obj@metadata$metric[[x]]$RLE$rle.data$rle.iqr
                iqr.qun.0.5 <- quantile(x = iqr, probs = .5)
                iqr.qun.0.1_0.9 <- quantile(x = iqr, probs = c(0.1,0.9))
                1 - c(iqr.qun.0.1_0.9[[2]] - iqr.qun.0.1_0.9[[1]])/iqr.qun.0.5
            })
        names(rle.iqr.scores) <- levels(assay.names)

    } else {
        rle.med.scores <- NULL
        rle.iqr.scores <- NULL
    }

    ## rle medians and variable scores #####
    rle.med.variables <- assessments.table$Variables[metrics.table$Factors == 'rleMedians']
    if(!is.null(rle.med.variables)){
        # check
        rle.meds <- unlist(lapply(
            assay.names,
            function(x){
                if(is.null(se.obj@metadata$metric[[x]]$RLE$rle.data$rle.med))
                    return(x)
            }))
        if(!is.null(rle.meds)){
            stop(paste0('The RLE medians for the ',
                        paste0(rle.meds, collapse = ' & '),
                        ' data cannot be found in the SummarizedExperiment object.',
                        'Please run the assessVariation() function.'))
        }

        # scores
        rle.var.aov.scores <- NULL
        rle.var.corr.scores <- NULL
        for(i in rle.med.variables ){
            for(j in levels(assay.names)){
                if(class(colData(se.obj)[[i]]) %in% c('factor', 'character')){
                    iqr <- unlist(lapply(
                        unique(colData(se.obj)[[i]]),
                        function(x) {
                            index <- colData(se.obj)[[i]] == x
                            stats::IQR(x = se.obj@metadata$metric[[j]]$RLE$rle.data$rle.med[index])
                        }))
                    iqr.qun.0.5 <- stats::quantile(x = iqr, probs = .5)
                    iqr.qun.0.1_0.9 <- stats::quantile(x = iqr, probs = c(0.1,0.9))
                    rle.var.aov.scores[[i]][j] <- 1 - c(c(iqr.qun.0.1_0.9[[2]] - iqr.qun.0.1_0.9[[1]])/iqr.qun.0.5)
                }
                if (class(colData(se.obj)[[i]]) %in% c('numeric', 'integer')){
                    rle.var.corr.scores[[i]][j] <- abs(suppressWarnings(stats::cor.test(
                        x = se.obj@metadata$metric[[j]]$RLE$rle.data$rle.med,
                        y = colData(se.obj)[[i]], method = 'spearman')[[4]]))
                }
            }
        }
    } else if (is.null(rle.med.variables)){
        rle.var.aov.scores <- NULL
        rle.var.corr.scores <- NULL
    }


    ## vector correlation scores ####
    pc.vec.corr.vars <- assessments.table$Variables[assessments.table$Metrics == 'VCA']
    if(!is.null(pc.vec.corr.vars)){
        # check
        for(i in pc.vec.corr.vars){
            vca <- unlist(lapply(
                assay.names,
                function(x){
                    if(is.null(se.obj@metadata$metric[[x]]$VCA[[i]]$vec.cor))
                        return(x)
                }))
            if(!is.null(vca)){
                stop(paste0('The vector correlation for the variable ', i, ' for the ',
                            paste0(vca, collapse = ' & '),
                            ' data cannot be found in the SummarizedExperiment object.',
                            'Please run the assessVariation() function.'))
            }
        }
        # scores
        pc.vec.corr.scores <- NULL
        for(i in pc.vec.corr.vars){
            for(j in levels(assay.names))
                pc.vec.corr.scores[[i]][j] <- mean(se.obj@metadata$metric[[j]]$VCA[[i]]$vec.cor)
        }
    } else pc.vec.corr.scores <- NULL


    ## linear regression scores ####
    pc.reg.vars <- assessments.table$Variables[assessments.table$Metrics == 'LRA']
    if(!is.null(pc.reg.vars)){
        # check
        for(i in pc.reg.vars){
            lra <- unlist(lapply(
                assay.names,
                function(x){
                    if(is.null(se.obj@metadata$metric[[j]]$LRA[[i]]$rseq))
                        return(x)
                }))
            if(!is.null(lra)){
                stop(paste0('The linear regression for the variable ', i, ' for the ',
                            paste0(vca, collapse = ' & '),
                            ' data cannot be found in the SummarizedExperiment object.',
                            'Please run the assessVariation() function.'))
            }
        }
        # scores
        pc.reg.scores <- NULL
        for(i in pc.reg.vars){
            for(j in levels(assay.names))
                pc.reg.scores[[i]][j] <- mean(se.obj@metadata$metric[[j]]$LRA[[i]]$rseq)
        }
    } else pc.reg.scores <- NULL

    ## silhouette scores ####
    all.sil.vars <- assessments.table$Metrics == 'Silhouette' & assessments.table$PlotTypes != 'combinedPlot'
    all.sil.vars <- assessments.table$Variables[all.sil.vars]
    if(!is.null(all.sil.vars)){
        # check
        for(i in all.sil.vars){
            sil.coeff <- unlist(lapply(
                assay.names,
                function(x){
                    if(is.null(se.obj@metadata$metric[[j]]$Silhouette[[paste0('sil.', sil.dist.measure)]][[i]]$sil.coef))
                        return(x)
                }))
            if(!is.null(sil.coeff)){
                stop(paste0('The silhouette coefficient for the variable ', i, ' for the ',
                            paste0(vca, collapse = ' & '),
                            ' data cannot be found in the SummarizedExperiment object.',
                            'Please run the assessVariation() function.'))
            }
        }
        # scores
        sil.scores <- NULL
        for(i in all.sil.vars){
            for(j in assay.names){
                x <- se.obj@metadata$metric[[j]]$Silhouette[[paste0('sil.', sil.dist.measure)]][[i]]$sil.coef
                sil.scores[[i]][j] <- c(x+1)/2
            }
        }
    } else sil.scores <- NULL


    ## ari scores ####
    all.ari.vars <- assessments.table$Metrics == 'ARI' & assessments.table$PlotTypes != 'combinedPlot'
    all.ari.vars <- assessments.table$Variables[ all.ari.vars]
    if(!is.null(all.ari.vars)){
        ari.scores <- NULL
        if(ari.clustering.method == 'mclust'){
            out.put.name <- 'mclust'
        } else out.put.name <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)

        for(i in all.sil.vars){
            # check
            for(i in all.sil.vars){
                ari.coeff <- unlist(lapply(
                    assay.names,
                    function(x){
                        if(is.null(se.obj@metadata$metric[[j]]$ARI[[out.put.name]][[i]]$ari))
                            return(x)
                    }))
                if(!is.null(ari.coeff)){
                    stop(paste0('The ARI for the variable ', i, ' for the ',
                                paste0(vca, collapse = ' & '),
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the assessVariation() function.'))
                }
            }
            # scores
            for(j in assay.names)
                ari.scores[[i]][j] <- se.obj@metadata$metric[[j]]$ARI[[out.put.name]][[i]]$ari
        }
    } else ari.scores <- NULL


    ## gene variable correlations scores ####
    gene.var.corr.vars <- assessments.table$Variables[assessments.table$Factors == 'geneCorr']
    if(!is.null(gene.var.corr.vars)){
        gene.var.corr.coef.scores <- NULL
        gene.var.corr.pvalue.scores <- NULL
        gene.var.corr.qvalue.scores <- NULL
        for(i in gene.var.corr.vars ){
            # check
            for(v in gene.var.corr.vars){
                corr.coeff <- unlist(lapply(
                    assay.names,
                    function(x){
                        if(is.null(se.obj@metadata$metric[[x]]$Correlation[[corr.method]][[v]]$cor.coef))
                            return(x)
                    }))
                if(!is.null(corr.coeff)){
                    stop(paste0('The gene-level correlations for the variable ', i, ' for the ',
                                paste0(vca, collapse = ' & '),
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the assessVariation() function.'))
                }
            }
            # scores
            for(j in assay.names){
                cor.coef <- abs(se.obj@metadata$metric[[j]]$Correlation[[corr.method]][[i]]$cor.coef[,2])
                gene.var.corr.coef.scores[[i]][j] <- sum(cor.coef < corr.coef.cutoff)/nrow(se.obj)
            }
            for(j in assay.names){
                p.values <- se.obj@metadata$metric[[j]]$Correlation[[corr.method]][[i]]$cor.coef[,1]
                gene.var.corr.pvalue.scores[[i]][j] <- ks.test(x = p.values[p.values > 0.05], "punif")$statistic
            }
            for(j in assay.names){
                p.values <- se.obj@metadata$metric[[j]]$Correlation[[corr.method]][[i]]$cor.coef[,1]
                gene.var.corr.qvalue.scores[[i]][j] <- qvalue::qvalue(p = p.values)$pi0
            }
        }
    } else {
        gene.var.corr.coef.scores <- NULL
        gene.var.corr.pvalue.scores <- NULL
        gene.var.corr.qvalue.scores <- NULL
    }

    ## gene variable anova scores ####
    index <- assessments.table$Factors == 'geneAnova'
    gene.var.anova.vars <- assessments.table$Variables[index]
    if(!is.null(gene.var.anova.vars)){
        gene.var.anova.coef.scores <- NULL
        gene.var.anova.pvalue.scores <- NULL
        gene.var.anova.qvalue.scores <- NULL
        for(i in gene.var.anova.vars ){
            # check
            for(v in gene.var.anova.vars){
                corr.coeff <- unlist(lapply(
                    assay.names,
                    function(x){
                        if(is.null(se.obj@metadata$metric[[j]]$ANOVA[[anova.method]][[i]]$F.values))
                            return(x)
                    }))
                if(!is.null(corr.coeff)){
                    stop(paste0('The gene-level ANOVA for the variable ', i, ' for the ',
                                paste0(vca, collapse = ' & '),
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the assessVariation() function.'))
                }
            }
            # scores
            for(j in assay.names){
                aov.fvalues <- se.obj@metadata$metric[[j]]$ANOVA[[anova.method]][[i]]$F.values$statistic
                gene.var.anova.fvalue.scores[[i]][j] <- sum(aov.fvalues < fstat.cutoff)/nrow(se.obj)
            }
            for(j in assay.names){
                aov.pvalues <- se.obj@metadata$metric[[j]]$ANOVA[[anova.method]][[i]]$F.values$pvalue
                gene.var.anova.pvalue.scores[[i]][j] <- ks.test(x = aov.pvalues[aov.pvalues > 0.05], "punif")$statistic
            }
            for(j in assay.names){
                aov.pvalues <- se.obj@metadata$metric[[j]]$ANOVA[[anova.method]][[i]]$F.values$pvalue
                gene.var.anova.qvalue.scores[[i]][j] <- qvalue::qvalue(p = aov.pvalues)$pi0
            }
        }
    } else {
        gene.var.anova.fvalue.scores <- NULL
        gene.var.anova.pvalue.scores <- NULL
        gene.var.anova.qvalue.scores <- NULL
    }


    ## dge scores  ####
    dge.qvalue.scores <- NULL
    dge.pvalue.scores <- NULL
    dge.vars <- assessments.table$Variables[metrics.table$Metrics == 'DGE' & assessments.table$Assessments != 'Exc']
    if(!is.null(dge.vars)){
        for(i in dge.vars){
            # check
            for(v in dge.vars){
                corr.coeff <- unlist(lapply(
                    assay.names,
                    function(x){
                        if(is.null(se.obj@metadata$metric[[x]]$DGE[[i]]$p.values))
                            return(x)
                    }))
                if(!is.null(corr.coeff)){
                    stop(paste0('The DEG for the variable ', i, ' for the ',
                                paste0(vca, collapse = ' & '),
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the assessVariation() function.'))
                }
            }
            # scores
            for(j in assay.names){
                p.values <- se.obj@metadata$metric[[j]]$DGE[[i]]$p.values
                dge.qvalue.scores[[i]][j] <- mean(sapply(
                    names(p.values),
                    function(x) qvalue::qvalue(p = p.values[[x]]$pvalue)$pi0))
            }
            for(j in assay.names){
                p.values <- se.obj@metadata$metric[[j]]$DGE[[i]]$p.values
                dge.pvalue.scores[[i]][j] <- mean(sapply(
                    names(p.values),
                    function(x)
                        ks.test(x = p.values[[x]]$pvalue[p.values[[x]]$pvalue > .05], "punif")$statistic
                ))
            }
        }
    } else {
        dge.qvalue.scores <- NULL
        dge.pvalue.scores <- NULL
    }

    ## bio genes ####

    # Final performance assessments ####
    test <- data <- NULL
    cols.names <- c('data', 'variable', 'test', 'measurements')
    ## RLE ####
    ### rle medians ####
    if(!is.null(rle.med.scores)){
        fa.rle.med.scores <- as.data.frame(rle.med.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            dplyr::mutate(variable = 'RLE') %>%
            dplyr::mutate(test = 'RLE medians') %>%
            data.frame(.)
        row.names(fa.rle.med.scores) <- c(1:nrow(fa.rle.med.scores))
        colnames(fa.rle.med.scores)[1] <- 'measurements'
        fa.rle.med.scores <- fa.rle.med.scores[ , cols.names]
    } else fa.rle.med.scores <- NULL

    ### rle IQR ####
    if(!is.null(rle.iqr.scores)){
        fa.rle.iqr.scores <- as.data.frame(rle.iqr.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            dplyr::mutate(variable = 'RLE') %>%
            dplyr::mutate(test = 'RLE IQRs') %>%
            data.frame(.)
        row.names(fa.rle.iqr.scores) <- c(1:nrow(fa.rle.iqr.scores))
        colnames(fa.rle.iqr.scores)[1] <- 'measurements'
        fa.rle.iqr.scores <- fa.rle.iqr.scores[ , cols.names]
    } else fa.rle.iqr.scores <- NULL

    ### rle variable correlation ####
    if(!is.null(rle.var.corr.scores)){
        fa.rle.var.corr.scores <- as.data.frame(rle.var.corr.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Correlation with RLE medains') %>%
            data.frame(.)
        fa.rle.var.corr.scores <- fa.rle.var.corr.scores[ , cols.names]
    } else fa.rle.var.corr.scores <- NULL

    ### rle variable ANOVA ####
    if(!is.null(rle.var.aov.scores)){
        fa.rle.var.aov.scores <- as.data.frame(rle.var.aov.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Association with RLE IQRs') %>%
            data.frame(.)
        fa.rle.var.aov.scores <- fa.rle.var.aov.scores[ , cols.names]
    } else fa.rle.var.aov.scores <- NULL


    ## PCA ####
    ### pcs regression ####
    if(!is.null(pc.reg.scores)){
        fa.pc.reg.scores <- as.data.frame(pc.reg.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'LRA') %>%
            data.frame(.)
        fa.pc.reg.scores <- fa.pc.reg.scores[ , cols.names]
    } else fa.pc.reg.scores <- NULL

    ### pcs vector correlation ####
    if(!is.null(pc.vec.corr.scores)){
        fa.pc.vec.corr.scores <- as.data.frame(pc.vec.corr.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'VCA') %>%
            data.frame(.)
        fa.pc.vec.corr.scores <- fa.pc.vec.corr.scores[ , cols.names]
    } else fa.pc.vec.corr.scores <- NULL

    ## Silhouette ####
    if(!is.null(sil.scores)){
        fa.sil.scores <- as.data.frame(sil.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Silhouette coefficient') %>%
            data.frame(.)
        fa.sil.scores <- fa.sil.scores[ , cols.names]
    } else fa.sil.scores <- NULL

    ## ARI ####
    if(!is.null(ari.scores)){
        fa.ari.scores <- as.data.frame(ari.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'ARI') %>%
            data.frame(.)
        fa.ari.scores <- fa.ari.scores[ , cols.names]
    } else fa.ari.scores <- NULL

    ## Gene variable correlation ####
    ### correlation coefficients ####
    if(!is.null(gene.var.corr.coef.scores)){
        fa.gene.var.corr.coef.scores <- as.data.frame(gene.var.corr.coef.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (correlation cutoff)') %>%
            data.frame()
        fa.gene.var.corr.coef.scores <- fa.gene.var.corr.coef.scores[ , cols.names]
    } else fa.gene.var.corr.coef.scores <- NULL
    ### p-values ####
    if(!is.null(gene.var.corr.pvalue.scores)){
        fa.gene.var.corr.pvalue.scores <- as.data.frame(gene.var.corr.pvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (p-values distribution)') %>%
            data.frame(.)
        fa.gene.var.corr.pvalue.scores <- fa.gene.var.corr.pvalue.scores[ ,cols.names ]
    } else fa.gene.var.corr.pvalue.scores <- NULL
    ### q-values ####
    if(!is.null(gene.var.corr.qvalue.scores)){
        fa.gene.var.corr.qvalue.scores <- as.data.frame(gene.var.corr.qvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (q-values)') %>%
            data.frame(.)
        fa.gene.var.corr.qvalue.scores <- fa.gene.var.corr.qvalue.scores[ ,cols.names ]
    } else fa.gene.var.corr.qvalue.scores <- NULL

    ## Gene variable ANOVA ####
    ### f-values ####
    if(!is.null(gene.var.anova.fvalue.scores)){
        fa.gene.var.anova.fvalue.scores <- as.data.frame(gene.var.anova.fvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (F-values cutoff)') %>%
            data.frame(.)
        fa.gene.var.anova.fvalue.scores <- fa.gene.var.anova.fvalue.scores[ , cols.names]
    } else fa.gene.var.anova.fvalue.scores <- NULL
    ### p-values ####
    if(!is.null(gene.var.anova.pvalue.scores)){
        fa.gene.var.anova.pvalue.scores <- as.data.frame(gene.var.anova.pvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (p-values distribution)') %>%
            data.frame(.)
        fa.gene.var.anova.pvalue.scores <- fa.gene.var.anova.pvalue.scores[ , cols.names]
    } else fa.gene.var.anova.pvalue.scores <- NULL
    ### q-values ####
    if(!is.null(gene.var.anova.qvalue.scores)){
        fa.gene.var.anova.qvalue.scores <- as.data.frame(gene.var.anova.qvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (q-values)') %>%
            data.frame(.)
        fa.gene.var.anova.qvalue.scores <- fa.gene.var.anova.qvalue.scores[ , cols.names]
    } else fa.gene.var.anova.qvalue.scores <- NULL

    ## DGE ####
    if(!is.null(dge.qvalue.scores)){
        fa.dge.qvalue.scores <- as.data.frame(dge.qvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'DGE (qvalue)') %>%
            data.frame()
        fa.dge.qvalue.scores <- fa.dge.qvalue.scores[ , cols.names]
    } else fa.dge.qvalue.scores <- NULL

    if(!is.null(dge.pvalue.scores)){
        fa.dge.pvalue.scores <- as.data.frame(dge.pvalue.scores) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'DGE (p-values distribution)')
        fa.dge.pvalue.scores <- fa.dge.pvalue.scores[ , cols.names]
    } else fa.dge.pvalue.scores <- NULL

    ## Put all measurements together ####
    all.measurements <- rbind(
        fa.rle.med.scores,
        fa.rle.iqr.scores,
        fa.rle.var.corr.scores,
        fa.rle.var.aov.scores,
        fa.pc.reg.scores,
        fa.pc.vec.corr.scores,
        fa.sil.scores,
        fa.ari.scores,
        fa.gene.var.corr.coef.scores,
        fa.gene.var.corr.pvalue.scores,
        fa.gene.var.corr.qvalue.scores,
        fa.gene.var.anova.fvalue.scores,
        fa.gene.var.anova.pvalue.scores,
        fa.gene.var.anova.qvalue.scores,
        fa.dge.qvalue.scores,
        fa.dge.pvalue.scores)
    all.measurements[,1] <- gsub(':', ".", all.measurements[,1])
    assay.names <- gsub(':', ".", assay.names)
    group <- NULL
    ### label biological and unwanted variables ####
    all.measurements$group <- 'Removal of unwanted variation'
    all.measurements$group[all.measurements$variable %in% bio.variables] <- 'Preservation of biological variation'

    ### change scores based on the biological and unwanted variables ####
    #### find categorical and continuous variables ####
    uv.categorical.vars <- uv.continuous.vars <- NULL
    if (!is.null(uv.variables)) {
        vars.class <- sapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]]))
        uv.categorical.vars <- names(vars.class[vars.class %in% c('character', 'factor')])
        uv.continuous.vars <- names(vars.class[vars.class %in% c('numeric', 'integer')])
        #### change PC RLA and correlation ####
        if(!is.null(uv.continuous.vars)){
            for(i in uv.continuous.vars){
                index <- all.measurements$variable == i & all.measurements$test == 'LRA'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
                index <- all.measurements$variable == i & all.measurements$test == 'Correlation with RLE medains'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
                index <- all.measurements$variable == i & all.measurements$test == 'Gene-level correlation (p-values distribution)'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
            }
        }
        #### change PC VCA, ARI and Silhouette ####
        if(!is.null(uv.categorical.vars)){
            for(i in uv.categorical.vars){
                index <- all.measurements$variable == i & all.measurements$test == 'Silhouette coefficient'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
                index <- all.measurements$variable == i & all.measurements$test == 'ARI'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
                index <- all.measurements$variable == i & all.measurements$test == 'VCA'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
                index <- all.measurements$variable == i & all.measurements$test == 'Gene-level ANOVA (p-values distribution)'
                all.measurements$measurements[index] <- 1 - all.measurements$measurements[index]
            }
        }
    }

    ## average of biological and unwanted variables scores ######
    bio.uv.overall.scores <- all.measurements %>%
        group_by(data, group) %>%
        summarise_at(vars("measurements"), mean) %>%
        dplyr::mutate(test = 'Score') %>%
        dplyr::mutate(variable = 'Score') %>%
        data.frame(.)
    bio.uv.overall.scores <- bio.uv.overall.scores [ , colnames(all.measurements)]

    ### calculate final scores ####
    final.overall.scores <- lapply(
        assay.names,
        function(x){
            index.a <- bio.uv.overall.scores$data == x & bio.uv.overall.scores$group == 'Removal of unwanted variation'
            index.b <- bio.uv.overall.scores$data == x & bio.uv.overall.scores$group == 'Preservation of biological variation'
            c(.4*bio.uv.overall.scores$measurements[index.a] + 0.6*bio.uv.overall.scores$measurements[index.b])
        })
    names(final.overall.scores) <- assay.names
    final.overall.scores <- final.overall.scores %>%
        data.frame(.) %>%
        pivot_longer(everything(), names_to = 'data', values_to = 'measurements') %>%
        dplyr::mutate(variable = 'Score') %>%
        dplyr::mutate(group = 'Final performance') %>%
        dplyr::mutate(test = 'Score') %>%
        data.frame()
    final.overall.scores <- final.overall.scores[ , colnames(all.measurements)]

    ###
    group.test.variable <- measurements <- test <- NULL
    all.measurements <- rbind(
        all.measurements,
        bio.uv.overall.scores,
        final.overall.scores)
    final.overall.scores <- final.overall.scores[order(final.overall.scores$measurements, decreasing = F) , ]
    all.measurements$data <- factor(x = all.measurements$data, levels = final.overall.scores$data)
    all.measurements$variable <- factor(
        x = all.measurements$variable,
        levels = unique(all.measurements$variable))
    all.measurements$group.test.variable <- paste0(
        all.measurements$group,
        all.measurements$test,
        all.measurements$variable)
    all.measurements <- all.measurements %>%
        group_by(group.test.variable) %>%
        mutate(rank = row_number(-measurements)) %>%
        data.frame(.)
    all.measurements$rank <- as.factor(all.measurements$rank )
    strip_background <- strip_nested(
        text_x = elem_list_text(colour = "white", face = "bold"),
        background_x =
            elem_list_rect(
                fill = c(
                    # level1 colors
                    case_match(
                        unique(all.measurements$group),
                        "Final performance" ~ "purple",
                        "Removal of unwanted variation" ~ "orange",
                        .default = "blue4"
                    ),
                    # level2 colors
                    case_match(
                        all.measurements$group,
                        "RLE" ~ "grey",
                        .default = "grey")
                ))
    )
    p.overall <- ggplot(data = all.measurements, aes(x = test, y = data)) +
        geom_point(aes(size = measurements, color = rank)) +
        scale_color_manual(values = RColorBrewer::brewer.pal(n = 8, name = 'Oranges')[8:1]) +
        theme_bw()  +
        xlab('') +
        ylab('') +
        facet_nested(. ~ group + variable, scales = "free_x", strip = strip_background) +
        theme(panel.spacing = unit(0, "lines"),
              panel.grid = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = 'black', linewidth = 1),
              strip.background = element_rect(color = "grey30", fill = "grey90"),
              strip.text = element_text(size = c(12)),
              panel.border = element_rect(color = "grey90"),
              axis.line.y = element_line(colour = 'white', linewidth = 1),
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_stripes(odd = "#00000000", even = "#33333333") +
        guides(color = guide_legend(override.aes = list(size = 5)))
    if(isTRUE(plot.output)) print(p.overall)
    printColoredMessage(message = '------------The normAssessment function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
