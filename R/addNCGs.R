#' Add pre-selected set of negative control genes to SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function adds a pre-selected set of negative control genes to SummarizedExperiment object.

#' @param se.obj A SummarizedExperiment object.
#' @param ncg A logical or vector of gene names of pre-selected genes as NCGs.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not. This involves principal
#' component analysis on the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs' first principal
#' components and with biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected genes as NCG needs to be evaluated. The default is NULL. This
#' means all the variables in the 'uv.variables' will be assessed.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the
#' performance of NCGs. The defaut is 5.
#' @param center Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the
#' performance of NCGs. The defaut is 5.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Symbol. Indicates
#' @param output.name Symbol. Indicates the name of the output in the meta data of the SummarizedExperiment object. The
#' default is 'NULL'. This means the function will create a name based the specified argument.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix rowSums colSums
#' @importFrom ruv replicate.matrix
#' @export

addNCGs <- function(
        se.obj,
        ncg,
        assay.name = NULL,
        assess.ncg = FALSE,
        variables.to.assess.ncg = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        output.name = NULL,
        verbose = TRUE
        ){

    # Check inputs ####
    if(is.logical(ncg)){
        if(length(ncg) > nrow(se.obj)){
            stop('The length of the "ncg" is larger than the number of rows in the SummarizedExperiment object.')
        } else if (sum(ncg) == 0){
            stop('The "ncg" does not contain any "TRUE" value.')
        }
        printColoredMessage(
            message = paste0(sum(ncg), 'NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose)
    } else if (is.character(ncg) | is.factor(ncg)){
        ncg <- intersect(unique(ncg), row.names(se.obj))
        if(length(ncg) == 0){
            stop('None of the "ncg" can be found in the SummarizedExperiment object. ')
        }
        printColoredMessage(
            message = paste0(length(ncg), 'NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose)
        ncg <- row.names(se.obj) %in% ncg
    } else if (is.numeric(ncg)){
        if(max(ncg) > nrow(se.obj)){
            stop('The "ncg" contains some number(s) that is larger than the number of rows in the SummarizedExperiment object.')
        }
        printColoredMessage(
            message = paste0(length(ncg), ' NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose)
        ncg.log <- rep(FALSE, nrow(se.obj))
        ncg.log[ncg] <- TRUE
        ncg <- ncg.log
    }

    if(isTRUE(assess.ncg)){
        if(is.null(variables.to.assess.ncg)){
            stop('The "variables.to.assess.ncg" must be provided to assess the performance of the NCG')
        }
        if(!variables.to.assess.ncg %in% colnames(colData(se.obj)) ){
            stop('Some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
        }
        if(is.null(assay.name)){
            stop('The "assay.name" must be provided to assess the performance of the NCG')
        }
    }

    # assessment of performance of the of NCGs  ####
    if(isTRUE(assess.ncg)){
        ## apply log ####
        if (isTRUE(apply.log) & !is.null(pseudo.count)){
            printColoredMessage(
                message = paste0(
                    'Applying log2 + ',
                    pseudo.count,
                    ' (pseudo.count) on the ',
                    assay.name,
                    ' data.'),
                color = 'blue',
                verbose = verbose)
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log) & is.null(pseudo.count)){
            printColoredMessage(
                message = paste0(
                    'Applying log2 on the ',
                    assay.name,
                    ' data.'),
                color = 'blue',
                verbose = verbose)
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose)
            expr.data <- assay(x = se.obj, i = assay.name)
        }

        # assessment of performance of the of NCGs  ####
        printColoredMessage(
            message = '-- Assess the performance of selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = paste0('Perform PCA on only selected genes as NCG.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Explore the association of the first ',
                nb.pcs,
                '  PCs with the ',
                paste0(variables.to.assess.ncg, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale)$u
        all.corr <- lapply(
            variables.to.assess.ncg,
            function(x){
                if(class(se.obj[[x]]) %in% c('numeric', 'integer')){
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared)
                } else if(class(se.obj[[x]]) %in% c('factor', 'character')){
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y){ cca <- cancor(
                            x = pca.data[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)
                        })
                }
            })
        names(all.corr) <- variables.to.assess.ncg
        pca.ncg <- as.data.frame(do.call(cbind, all.corr))
        pcs <- Groups <- NULL
        pca.ncg['pcs'] <- c(1:nb.pcs)
        pca.ncg <- tidyr::pivot_longer(
            data = pca.ncg,
            -pcs,
            names_to = 'Groups',
            values_to = 'ls')
        pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), linewidth = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab (expression("Correlations")) +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                plot.title = element_text(size = 16)
            )
        if(verbose)  print(pca.ncg)
    }

    # save the results ####
    if(is.null(output.name)){
        output.name <- paste0(sum(ncg), ' genes')
    }
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if(!'pre.selected' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['pre.selected']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
            se.obj@metadata[['NCG']][['pre.selected']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['pre.selected']][[output.name]] <- ncg
        printColoredMessage(
            message = 'The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The supervisedFindNcgTWAnova function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    }
}
