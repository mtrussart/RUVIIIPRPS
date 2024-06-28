#' Plot the the variation of biological and unwanted variables.

#' @param se.obj #' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is set to 'all', so all the assays in the SummarizedExperiment object will
#' be selected.
#' @param variables Symbol. A symbol or vector of symbols indicating the name(s) of the column(s) in the SummarizedExperiment
#' object for which the specified metrics below have been calculated.

#' @param fast.pca Logical. Indicates whether to the fast SVD approach has been used in the 'computePCA' function or not.
#' The default is set to 'TRUE'.
#' @param anova.method Symbol. A symbol indicating which anova method has been used to compute association between
#' gene-level expression and a categorical variable. The options are 'aov' and 'welch'. The default is set to 'aov'.
#' @param corr.method Symbol. A symbol indicating which correlation method should be used to compute correlation between
#' gene-level expression and a continuous variable. The options are 'pearson', 'kendall', or "spearman". The default is
#' set to 'spearman'.
#' @param pcorr.method Symbol. A symbol indicating which correlation method should be used to compute gene-gene partial
#' correlation. The options are 'pearson', 'kendall', or "spearman". The default is set to 'spearman'.
#' @param sil.dist.measure Symbol. A symbol indicating which ditsance measure to be applied on the PCs to calculate silhouette.
#' The options are 'euclidean', maximum', manhattan', 'canberra', 'binary' or 'minkowski'. The default is set to 'euclidean'.
#' Refer to the function 'computeSilhouette' for more details.
#' @param ari.clustering.method Symbol. A symbol that indicates which clustering methods should be applied on the PCs to
#' calculate the ARI. The options are 'mclust' or 'hclust' methods. The default is set to 'hclust'. Refer to the 'computeARI'
#' for more details.
#' @param ari.hclust.method Symbol. A symbol specifying the agglomeration method to be used when the 'clustering.method'
#' is set to 'hclust' method for computing ARI. The options are: 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA),
#' 'mcquitty' (= WPGMA), median' (= WPGMC) or centroid' (= UPGMC). The default is set to 'complete'. Refer to the 'computeARI'
#' for more details.
#' @param ari.hclust.dist.measure Symbol. A symbol indicating the distance measure to be used in the 'dist' function when
#' the 'clustering.method' is set to 'hclust' method for computing ARI. The options are 'euclidean', 'maximum', 'manhattan',
#' 'canberra', 'binary' or 'minkowski'. The default is set to 'euclidean'. Refer to the 'computeARI' for more details.
#' @param output.file.name Path and name of the output file to save the assessments plots in a pdf format.
#' @param pdf.width TTT
#' @param pdf.height TTT
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

plotAssessVariation <- function(
        se.obj,
        assay.names = 'all',
        variables,
        fast.pca = TRUE,
        anova.method = 'aov',
        corr.method = 'spearman',
        pcorr.method = 'spearman',
        sil.dist.measure = 'euclidian',
        ari.clustering.method = "hclust",
        ari.hclust.method = "complete",
        ari.hclust.dist.measure = "euclidian",
        output.file.name = NULL,
        pdf.width = 12,
        pdf.height = 12,
        verbose = TRUE
    ){

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    metrics.table <- se.obj@metadata$AssessmentMetrics$metrics.table

    # Put all plots together ####
    ## find classes of different variables ####
    categorical.vars <- continuous.vars <- NULL
    if (!is.null(variables)) {
        vars.class <- sapply(
            variables,
            function(x) class(colData(se.obj)[[x]]))
        categorical.vars <- names(vars.class[vars.class %in% c('character', 'factor')])
        continuous.vars <- names(vars.class[vars.class %in% c('numeric', 'integer')])
    }
    ## select output file names ####
    if(is.null(output.file.name)){
        output.file.name <- 'RUVIIIPRPS_AssessVariation'
    }
    pdf(paste0(output.file.name, '.pdf'),
        width = pdf.width,
        height = pdf.height
    )
    plot.new()
    text(.5, .7, "Assess variation", font = 2, cex = 2.5)
    text(.5, .6, "variables:", font = 2, cex = 2.5)
    text(.5, .4, paste0(variables, collapse = '\n'), font = 2, cex = 2)
    print(se.obj@metadata$AssessmentMetrics$plot)

    ## general RLE plot ####
    if('rlePlot' %in% metrics.table$PlotTypes) {
        if(length(assay.names) > 1){
            print(
                se.obj@metadata$Plots$global.level$RLE$gene.median.center$general$un.colored
            )
        }
    }
    ## continuous variables ####
    for(i in continuous.vars){
        plot.new()
        text(.5, .7, paste0("Assess variation \n in the variable: \n ", i ), font = 2, cex = 2.5)
        metrics.table.var <- metrics.table[metrics.table$Variables == i, ]
        ### scatter plot between RLE medians and variable ####
        if('rleMedians' %in% metrics.table.var$Factors){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.medians.variable[[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.medians.variable[[i]]$scatter.plot
            )
        }
        ## scatter plot between RLE iqrs and variable ####
        if('rleIqr' %in% metrics.table.var$Factors ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.iqrs.variable[[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.iqr.variable[[i]]$scatter.plot
            )
        }
        ## scatter plot between PCs and variable ####
        if(isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        if('PCA' %in% metrics.table.var$Metrics){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$scatter.plot
            )
        }
        ## line-dot plot for linear regression ####
        if('LRA' %in% metrics.table.var$Metrics){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$LRA[[svd.method]][[i]]$line.dotplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$LRA[[svd.method]][[i]]$plot
            )
        }
        ## boxplot of gene variable correlation coefficients ####
        if('CorrelationboxPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
            )
        }
        ## histograms of gene variable correlation coefficients ####
        if('CorrelationpvalHist' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
            )
        }
        ## scatter plot of gene variable partial correlation coefficients ####
        if('PartialCorrelationscatterPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$scatter.plot
            )
        }
        ## barplot of gene variable partial correlation coefficients ####
        if('PartialCorrelationbarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$barplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$barplot
            )
        }
        ## histogram of gene variable partial correlation coefficients ####
        if('PartialCorrelationhistogram' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$histogram
            )
        }
        ## histogram of gene variable partial correlation coefficients ####
        if('GeneSetScore' %in% paste0(metrics.table.var$Metrics)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$GeneSetSocore$singscore[[i]]$general
                )
            }
        }
    }

    ## categorical variables ####
    for(i in categorical.vars){
        plot.new()
        text(.5, .7, paste0("Assess variation \n in the variable: \n ", i ), font = 2, cex = 2.5)
        metrics.table.var <- metrics.table[metrics.table$Variables == i, ]
        ### rle plots colored by the variable ####
        if('coloredRLEplot' %in% metrics.table.var$PlotTypes){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$gene.median.center[[i]]$colored
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$gene.median.center[[i]]$colored
            )
        }
        ## boxplot between RLE medians and variable ####
        if('rleMedians' %in% metrics.table.var$Factors){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.medians.variable[[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.medians.variable[[i]]$boxplot
            )
        }
        ## boxplot between RLE iqrs and variable ####
        if('rleIqr' %in% metrics.table.var$Factors){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.iqrs.variable[[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.iqr.variable[[i]]$boxplot
            )
        }
        ## scatter plot of PCs colored by the variable ####
        if(isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        if('pcsscatterPlot' %in% paste0(metrics.table.var$Factors, metrics.table.var$PlotTypes)){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$boxplot.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$boxplot.plot
            )
        }
        ## boxplot plot of PCs colored by the variable ####
        if('pcsboxPlot' %in% paste0(metrics.table.var$Factors, metrics.table.var$PlotTypes) ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$scatter.plot
            )
        }
        ## line-dot plot for vector correlation ####
        if('VCA' %in% metrics.table.var$Metrics ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$VCA[[svd.method]][[i]]$line.dotplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$VCA[[svd.method]][[i]]$plot
            )
        }
        ## barplot of ari ####
        if('ARIbarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)) {
            if(ari.clustering.method == 'mclust'){
                ari.method <- 'mclust'
            } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$ARI[[ari.method]][[i]]$single.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$ARI[[ari.method]][[i]]$single.plot
            )
        }
        ## barplot of silhouette ####
        if('SilhouettebarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            silhouette.method <- paste0('sil.', sil.dist.measure)
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$Silhouette[[silhouette.method]][[i]]$single.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$Silhouette[[silhouette.method]][[i]]$single.plot
            )
        }
        ## boxplot of gene variable ANOVA F-stat ####
        if('ANOVAboxPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes) ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$ANOVA[[anova.method]][[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$ANOVA[[anova.method]][[i]]$boxplot
            )
        }
        ## histograms of gene variable ANOVA p-values ####
        if('ANOVApvalHist' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes) ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$ANOVA[[anova.method]][[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$ANOVA[[anova.method]][[i]]$histogram
            )
        }
        ## histograms of differential gene expression p-values  ####
        if('DGE' %in% metrics.table.var$Metrics ){
            if(length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$DGE$Wilcoxon[[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$DGE$Wilcoxon$time.interval$plot
            )
        }
    }
    dev.off()
}


