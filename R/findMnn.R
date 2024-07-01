#' Find mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function finds mutual nearest neighbors between all pairs of batches in RNA-seq data. The mutual nearest neighbors
#' will be used to find and create pseudo samples and eventually pseudo-replicates for RUV-III normalization.

#' @details
#' Additional details...

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol that indicates the name of the assay in the SummarizedExperiment object. This data
#' should be the one that will be used as input for RUV-III normalization.
#' @param uv.variable Symbol. A symbol that indicates the name of the column in the SummarizedExperiment object. The
#' 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
#' divided into 'nb.clusters' groups using the 'clustering.method' methdo.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find mutual nearest
#' neighbors samples across the batches. The default is set to 'NULL'. The 'findBioGenes' function can be used for specify
#' a set of genes.
#' @param clustering.method Symbol. A symbol that indicates the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param mnn Numeric. The maximum number of nearest neighbors to compute. The default is set 1.
#' @param normalization Symbol. A Symbol that indicates which normalization methods should be applied before finding the
#' knn. The option are 'CPM', 'TMM', 'VST'. The default is 'CPM'. Refer to the 'applyOtherNormalization' for more details.
#' @param regress.out.bio.variables Symbol. A symbol or vector of symbols indicating the column name(s) that contain
#' biological variable(s) in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding genes that are highly affected by unwanted variation variable. The default is set to 'NULL', indicates the
#' regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is SerialParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '||' , assay.name)".
#' @param prps.group TT
#' @param plot.output Logical. If 'TRUE', the function plots the distribution of MNN across the batches.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN
#' @importFrom BiocParallel SerialParam
#' @importFrom utils setTxtProgressBar
#' @export

findMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        hvg = NULL,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        mnn = 1,
        normalization = 'CPM',
        regress.out.bio.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        plot.output = TRUE,
        output.name = NULL,
        prps.group = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findMnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs #####
    if (is.list(assay.name)) {
        stop('The "assay.name" cannot be a list.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be the name of signle assay in the SummarizedExperiment object.')
    } else if (is.null(uv.variable)) {
        stop('The "uv.variable" variable cannot be empty.')
    } else if (length(uv.variable) > 1) {
        stop('The "uv.variable" must contain the name of signle variable in the SummarizedExperiment object.')
    } else if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" variable cannot be found in the SummarizedExperiment object.')
    }
    if (!is.null(hvg)) {
        if (sum(hvg %in% row.names(se.obj)) != length(hvg))
            stop('All the hvg genes are not found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment object #####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Keep the original sample orders and the variable ####
    all.samples.index <- c(1:ncol(se.obj))
    ini.variable <- se.obj[[uv.variable]]
    ini.sample.names <- colnames(se.obj)
    colnames(se.obj) <- paste0('sample_', seq(ncol(se.obj)))

    # Assess and group the unwanted variable ####
    printColoredMessage(
        message = '-- Assess and group the unwanted variable:',
        color = 'magenta',
        verbose = verbose
    )
    if (class(se.obj[[uv.variable]]) %in% c('integer', 'numeric')) {
        ## continuous variable ####
        printColoredMessage(
            message = paste0(
                '- The "', uv.variable, '" is a continouse variable, then it will be divided into ',
                nb.clusters, ' groups using the ', clustering.method,' clustering.'),
            color = 'blue',
            verbose = verbose
        )
        ## kmeans ####
        if (clustering.method == 'kmeans') {
            set.seed(3456)
            uv.cont.clusters <- stats::kmeans(
                x = colData(se.obj)[[uv.variable]],
                centers = nb.clusters,
                iter.max = 1000
            )
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters$cluster))
        }
        ## cut ####
        if (clustering.method == 'cut') {
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = nb.clusters,
                include.lowest = TRUE
            ))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters))
        }
        ## quantile ####
        if (clustering.method == 'quantile') {
            quantiles <- stats::quantile(
                x = colData(se.obj)[[uv.variable]],
                probs = seq(0, 1, 1 / nb.clusters)
                )
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = quantiles,
                include.lowest = TRUE
            ))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters))
        }
    }
    ## categorical variable ####
    if (is.factor(se.obj[[uv.variable]]) | is.character(se.obj[[uv.variable]])) {
        if(length(unique(se.obj[[uv.variable]])) == 1){
            stop('To create MNN, the "uv.variable" must have at least two groups/levels.')
        } else if (length(unique(se.obj[[uv.variable]])) > 1){
            printColoredMessage(
                message = paste0('- The ', uv.variable, ' is a categorical variable with  ',
                                 length(unique(se.obj[[uv.variable]])), ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
        }
    }

    # Data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
        )
    groups <- levels(se.obj[[uv.variable]])
    all.norm.data <- lapply(
        groups,
        function(x) {
            selected.samples <- colData(se.obj)[[uv.variable]] == x
            if (!is.null(normalization) & is.null(regress.out.bio.variables)) {
                printColoredMessage(
                    message = paste0('- Apply the ', normalization, ' within the samples from the "', x, '" group.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- applyOtherNormalizations(
                        se.obj = se.obj[, selected.samples],
                        assay.name = assay.name,
                        method = normalization,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        assess.se.obj = FALSE,
                        save.se.obj = FALSE,
                        remove.na = 'none',
                        verbose = FALSE
                    )
                norm.data
            } else if (!is.null(normalization) & !is.null(regress.out.bio.variables)) {
                printColoredMessage(
                    message = paste0('- Apply the ', normalization, ' within the samples from "', x,
                        '" group and then regressing out ', paste0(regress.out.bio.variables, collapse = '&'), ' from the data.'),
                    color = 'blue',
                    verbose = verbose
                )
                ## normalization ####
                norm.data <- applyOtherNormalizations(
                        se.obj = se.obj[, selected.samples],
                        assay.name = assay.name,
                        method = normalization,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        assess.se.obj = FALSE,
                        save.se.obj = FALSE,
                        remove.na = 'none',
                        verbose = FALSE
                    )
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                # regression ####
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.bio.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)
                norm.data
            } else if (is.null(normalization) & !is.null(regress.out.bio.variables)){
                if(isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0('- Apply log and then regress out ',
                                         paste0(regress.out.bio.variables, collapse = '&'), x, '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if(!is.null(pseudo.count)){
                        norm.data <- log2(assay(se.obj[, selected.samples], assay.name) + pseudo.count)
                        norm.data
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0('- Regress out ',
                                         paste0(regress.out.bio.variables, collapse = '&'), x, '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                    norm.data
                }
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                # regression ####
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.bio.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)
                norm.data

            } else if (is.null(normalization) & is.null(regress.out.bio.variables)) {
                if(isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0('- Apply the log2 within the samples from "', x,
                                         '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if(!is.null(pseudo.count)){
                        norm.data <- log2(assay(se.obj[, selected.samples], assay.name) + pseudo.count)
                        norm.data
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0('No any library size normalization and transformation ', 'within samples from ',
                                         unique(colData(se.obj)[[uv.variable]][x]), ', and just finding ', ' mnn for individual samples.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                    norm.data
                }
            }
        })
    names(all.norm.data) <- groups

    # Find MNN between batches ####
    printColoredMessage(
        message = '-- Find MNN across all possible pair of batches, this may take few minuets:',
        color = 'magenta',
        verbose = verbose
        )
    pairs.batch <- combn(x = groups, m = 2)
    printColoredMessage(
        message = paste0(
            '- All MNN between all ', ncol(pairs.batch),
            ' pairs of the sub-groups of the "', uv.variable, '" variable will be found:'),
        color = 'orange',
        verbose = verbose
        )
    pb <- utils::txtProgressBar(
        min = 0,
        max = ncol(pairs.batch),
        style = 3
        )
    all.mnn <- lapply(
        1:ncol(pairs.batch),
        function(x) {
            message(' ')
            printColoredMessage(
                message = paste0(
                    '* Find MNN between the "', pairs.batch[1, x], '" and the "', pairs.batch[2, x], '" sub-groups.'),
                color = 'blue',
                verbose = verbose
            )
            if (is.null(hvg)) {
                mnn.samples <- BiocNeighbors::findMutualNN(
                    data1 = t(all.norm.data[[pairs.batch[1, x]]]),
                    data2 = t(all.norm.data[[pairs.batch[2, x]]]),
                    k1 = mnn,
                    k2 = mnn,
                    BPPARAM = mnn.bpparam
                )
            } else if (!is.null(hvg)) {
                mnn.samples <- findMutualNN(
                    data1 = t(all.norm.data[[pairs.batch[1, x]]][hvg , ]),
                    data2 = t(all.norm.data[[pairs.batch[2, x]]][hvg , ]),
                    k1 = mnn,
                    k2 = mnn,
                    BPPARAM = mnn.bpparam
                )
            }
            if (is.null(mnn.samples)){
                printColoredMessage(
                    message = '* MNN are not found.',
                    color = 'blue',
                    verbose = verbose
                )
            }
            group1 <- group2 <- NULL
            df <- data.frame(
                group1 = rep(pairs.batch[1, x], length(mnn.samples$first)),
                group2 = rep(pairs.batch[2, x], length(mnn.samples$second)),
                sample.no.1 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[1, x]][mnn.samples$first],
                sample.no.2 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[2, x]][mnn.samples$second],
                sample.no.a = as.numeric(gsub('sample_', '', colnames(all.norm.data[[pairs.batch[1, x]]])[mnn.samples$first])),
                sample.no.b = as.numeric(gsub('sample_', '', colnames(all.norm.data[[pairs.batch[2, x]]])[mnn.samples$second]))
            )
            # Sanity check ###
            if(!all.equal(df$sample.no.1, df$sample.no.1.1)){
                stop('There something wrong with the MNN.')
            }
            if(!all.equal(df$sample.no.2, df$sample.no.2.2)){
                stop('There something wrong with the MNN.')
            }

            printColoredMessage(
                message = paste0('* ', nrow(df), ' mnn are found.'),
                color = 'blue',
                verbose = verbose
            )
            setTxtProgressBar(pb, x)
            return(df)
        })
    all.mnn <- do.call(rbind, all.mnn)
    if (is.null(all.mnn)){
        stop('MNN are not found across any batches. You may want to increase the valeu of the mnn.')
    }

    # Plot the distribution of MNN  ####
    p.mnn <- ggplot(all.mnn, aes(x = group1, y = group2)) +
        geom_point() +
        geom_count() +
        ggtitle('Distribution of MNN across batches') +
        xlab('') +
        ylab('') +
        theme_bw() +
        theme(
            axis.line = element_line(colour = 'black', linewidth = 1),
            plot.title = element_text(size = 12),
            axis.text.x = element_text(
                size = 12,
                angle = 35,
                vjust = 1,
                hjust = 1),
            axis.text.y = element_text(
                size = 12,
                angle = 35,
                vjust = 1,
                hjust = 1
            )
        )
    if (isTRUE(plot.output)) print(p.mnn)
    message(' ')
    printColoredMessage(
        message = paste0('- in total ' , nrow(all.mnn), ' mnn are found.'),
        color = 'blue',
        verbose = verbose
    )
    if (isTRUE(sum(is.na(all.mnn))))
        stop('There are NA in the MNN.')
    se.obj[[uv.variable]] <- ini.variable
    colnames(se.obj) <- ini.sample.names
    # Save the results ####
    ## select prps name ####
    if(is.null(prps.group)){
        prps.group <- paste0('prps|mnn|', uv.variable)
    }
    ## select output name ####
    if (is.null(output.name))
        output.name <- paste0(uv.variable, '|' , assay.name)

    ## save the results in the SummarizedExperiment object ####
    message(' ')
    printColoredMessage(message = '-- Save the results.',
                        color = 'magenta',
                        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        ## check
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## check
        if (!'un.supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['un.supervised']] <- list()
        }
        ## check
        if (!prps.group %in% names(se.obj@metadata[['PRPS']][['un.supervised']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]] <- list()
        }
        ## check
        if (!'KnnMnn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']] <- list()
        }
        ## check
        if (!'knn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['mnn']] <- list()
        }
        ## check
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['knn']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['mnn']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['mnn']][[output.name]] <- all.mnn
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['mnn']][['plot']] <- p.mnn
        printColoredMessage(
            message = '- All the mnn results are saved in the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The findMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the results as matrix  ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- All the mnn results are outputed as matrix.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The findMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(
            mnn = all.mnn,
            mnn.plot = p.mnn)
            )
    }
}
