#' Find mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function finds mutual nearest neighbors between all pairs of batches in RNA-seq data. The mutual nearest neighbors
#' will be used to find and create pseudo samples and eventually pseudo-replicates for RUV-III normalization.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol that indicates the name of the assay in the SummarizedExperiment object. This data
#' should be the one that will be used as input for RUV-III normalization.
#' @param uv.variable Symbol. A symbol indicating the name of the column in the SummarizedExperiment object. The
#' 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
#' divided into 'nb.clusters' groups using the 'clustering.method' method.
#' @param mnn Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute. The default
#' is set 1.
#' @param clustering.method Symbol. A symbol that indicates the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param normalization Symbol. A Symbol that indicates which normalization methods should be applied before finding the
#' knn. The option are 'CPM', 'TMM', 'VST'. The default is 'CPM'. Refer to the 'applyOtherNormalization' for more details.
#' @param apply.cosine.norm Logical. Idicates whether to apply cosine normalization on the data or not. The default is
#' set to 'FALSE'.
#' @param regress.out.variables Symbol. A symbol or vector of symbols indicating the column name(s) that contain
#' biological variable(s) in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding genes that are highly affected by unwanted variation variable. The default is set to 'NULL', indicates the
#' regression will not be applied.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find mutual nearest
#' neighbors samples across the batches. The default is set to 'NULL'. The 'findBioGenes' function can be used for specify
#' a set of genes.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is SerialParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param mnn.nbparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is KmknnParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param plot.output Logical. If 'TRUE', the function plots the distribution of MNN across the batches.
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '||' , assay.name)".
#' @param prps.group Symbol. A symbol specifying the name of the PRPS group. If is 'NULL', the function will select a name
#' based on "paste0('prps|mnn|', uv.variable)".
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom utils setTxtProgressBar
#' @export

findMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        nb.mnn = 1,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        normalization = 'CPM',
        apply.cosine.norm = TRUE,
        regress.out.variables = NULL,
        hvg = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
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
        if(is.logical(hvg)){
            if(sum(hvg) == 0){
                stop('The "hvg" does not contain any genes.')
            } else if (sum(hvg) <= 3){
                stop('The number of "hvg" must be at least 3')
            } else if (length(hvg) > nrow(se.obj)){
                stop('The length of the "hvg" cannot be larger than the row numbers of the SammarizedExperiment object.')
            }
        }
        if(is.character(hvg)){
            if (sum(hvg %in% row.names(se.obj)) != length(hvg)){
                stop('All the hvg genes are not found in the SummarizedExperiment object.')
            } else if (length(hvg) > nrow(se.obj)){
                stop('The length of the "hvg" cannot be larger than the row numbers of the SammarizedExperiment object.')
            } else if (length(hvg) <= 3){
                stop('The number of "hvg" must be at least 3')
            }
            hvg <- row.names(se.obj) %in% hvgs
        }
    } else if (is.null(hvg)){
        hvg <- rep(TRUE, nrow(se.obj))
    }

    # Assess the SummarizedExperiment object #####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Keeping the original sample orders and the unwanted variable ####
    all.samples.index <- c(1:ncol(se.obj))
    initial.variable <- se.obj[[uv.variable]]
    initial.sample.names <- colnames(se.obj)
    colnames(se.obj) <- paste0('sample_', seq(ncol(se.obj)))

    # Assess and group the unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
    )
    if(is.numeric(initial.variable)){
        se.obj[[uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_group',
            verbose = verbose
        )
    }
    if(!is.numeric(initial.variable)){
        length.variable <- length(unique(initial.variable))
        if( length.variable == 1){
            stop('To create MNN, the "uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
        }
    }

    # Checking sample sizes of each sub group ####
    printColoredMessage(
        message = '-- Checking the sample size of each sub-group of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    sub.group.sample.size <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = nb.mnn + 1
        )
    if (length(sub.group.sample.size) == 0){
        stop(paste0(
            'All subgroups of the unwanted variable have less than ',
            mnn + 1,
            ' samples. MNN cannot be found.'))
    } else if (length(sub.group.sample.size) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the unwanted variable have less than ',
                mnn + 1,
                ' samples. '),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the subgroups of the unwanted variable have at least, ',
                nb.mnn + 1,
                ' samples'),
            color = 'blue',
            verbose = verbose
        )
    }

    # Data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
    )
    all.norm.data <- lapply(
        sub.group.sample.size,
        function(x) {
            selected.samples <- colData(se.obj)[[uv.variable]] == x
            ## normalization ####
            if (!is.null(normalization) & is.null(regress.out.variables)) {
                printColoredMessage(
                    message = paste0(
                        '- applying the ',
                        normalization,
                        ' on the samples from the "',
                        x,
                        '" sub-group.'),
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
            }
            ## normalization and regression ####
            if (!is.null(normalization) & !is.null(regress.out.variables)) {
                printColoredMessage(
                    message = paste0(
                        '- applying the ',
                        normalization,
                        ' on the samples from "',
                        x,
                        '" group and then regressing out the ',
                        paste0(regress.out.variables, collapse = '&'),
                        ' variable(s) from the data.'),
                    color = 'blue',
                    verbose = verbose
                )
                ### normalization ####
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
                ## regression ####
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)

            }
            ## regression ####
            if (is.null(normalization) & !is.null(regress.out.variables)){
                if(isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- applying log transformation and then regressing out the ',
                            paste0(regress.out.variables, collapse = '&'),
                            ' variable(s) on the samples ',
                            x,
                            '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if(!is.null(pseudo.count)){
                        norm.data <- log2(assay(se.obj[, selected.samples], assay.name) + pseudo.count)
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- regress out ', paste0(regress.out.variables, collapse = '&'), x,
                            '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                }
                ### regression ####
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)
            }
            ## log transformation ####
            if (is.null(normalization) & is.null(regress.out.variables)) {
                if(isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- applying the log2 within the samples from "',
                            x,
                            '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if(!is.null(pseudo.count)){
                        norm.data <- log2(assay(x = se.obj[, selected.samples], i = assay.name) + pseudo.count)
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- no library size normalization and transformation is applied on samples from ', x, '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                }
            }
            return(norm.data)
        })
    names(all.norm.data) <- sub.group.sample.size

    # Selecting input data for knn analysis ####
    printColoredMessage(
        message = '-- Selecting the data input for knn analysis',
        color = 'magenta',
        verbose = verbose
    )
    all.data.input <- lapply(
        sub.group.sample.size,
        function(x){
            norm.data <- all.norm.data[[x]]
            ## data input: expression matrix with hvg #####
            if (data.input == 'expr' & !is.null(hvg)) {
                printColoredMessage(
                    message = '- selecting the gene expression matrix with the highly variable genes as the data input.',
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- t(norm.data[hvg,])
            }
            ## data input: expression matrix with all genes #####
            if (data.input == 'expr' & is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- selecting the gene expression matrix with all genes as the data input for',
                        'the sub-group: ',
                        x,
                        '.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- t(norm.data)
            }
            ## data input: PCA with hvg #####
            if (data.input == 'pcs' & !is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- performing PCA on the gene expression matrix using highly',
                        ' variable genes and using PCs as the data input.'),
                    color = 'blue',
                    verbose = verbose
                )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data[hvg,]),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
                )
                norm.data <- sv.dec$u
            }
            ## data input: PCA all genes #####
            if (data.input == 'pcs' & is.null(hvg)) {
                printColoredMessage(
                    message = '- performing PCA on the gene expression matrix and select PCs as the data input.',
                    color = 'blue',
                    verbose = verbose
                )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
                )
                norm.data <- sv.dec$u
            }
            return(norm.data)
        })
    names(all.data.input) <- sub.group.sample.size

    # Finding MNN between batches ####
    printColoredMessage(
        message = '-- Finding MNN across all possible pairs of sub-groups of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    pairs.batch <- combn(
        x = sub.group.sample.size,
        m = 2
        )
    printColoredMessage(
        message = paste0(
            '- All MNN between all ',
            ncol(pairs.batch),
            ' pairs of the sub-groups of the "',
            uv.variable,
            '" variable will be identified:'),
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
                    '* Finding MNN between the "',
                    pairs.batch[1, x],
                    '" and the "',
                    pairs.batch[2, x],
                    '" sub-groups.'),
                color = 'blue',
                verbose = verbose
                )
            # select data
            if(isTRUE(apply.cosine.norm)){
                printColoredMessage(
                    message = '--applying cosine normalization of the data.',
                    color = 'blue',
                    verbose = verbose
                    )
                data1 = cosineNorm(
                    x = all.data.input[[pairs.batch[1, x]]],
                    mode = 'matrix'
                    )
                data2 = cosineNorm(
                    x = all.data.input[[pairs.batch[2, x]]],
                    mode = 'matrix'
                    )
            } else {
                data1 = all.data.input[[pairs.batch[1, x]]]
                data2 = all.data.input[[pairs.batch[2, x]]]
            }
            if(!is.null(hvg)){
                printColoredMessage(
                    message = '--using HVG.',
                    color = 'blue',
                    verbose = verbose
                )
                data1 <- data1[ , hvg]
                data2 <- data2[ , hvg]
            }
            mnn.samples <- BiocNeighbors::findMutualNN(
                data1 = data1,
                data2 = data2,
                k1 = nb.mnn,
                k2 = nb.mnn,
                BPPARAM = mnn.bpparam,
                mnn.nbparam = KmknnParam()
                )
            if (is.null(mnn.samples)){
                stop('* MNN cannot be found.')
            }
            group1 <- group2 <- NULL
            mnn.data <- data.frame(
                group1 = rep(pairs.batch[1, x], length(mnn.samples$first)),
                group2 = rep(pairs.batch[2, x], length(mnn.samples$second)),
                sample.no.1 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[1, x]][mnn.samples$first],
                sample.no.2 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[2, x]][mnn.samples$second],
                sample.no.a = as.numeric(gsub('sample_', '', row.names(all.data.input[[pairs.batch[1, x]]])[mnn.samples$first])),
                sample.no.b = as.numeric(gsub('sample_', '', row.names(all.data.input[[pairs.batch[2, x]]])[mnn.samples$second]))
            )
            # Sanity check ###
            if(!all.equal(mnn.data$sample.no.1, mnn.data$sample.no.a)){
                stop('There something wrong with the MNN.')
            }
            if(!all.equal(mnn.data$sample.no.2, mnn.data$sample.no.b)){
                stop('There something wrong with the MNN.')
            }

            printColoredMessage(
                message = paste0('* ', nrow(mnn.data), ' mnn are found.'),
                color = 'blue',
                verbose = verbose
            )
            setTxtProgressBar(pb, x)
            return(mnn.data)
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
    se.obj[[uv.variable]] <- initial.variable
    colnames(se.obj) <- initial.sample.names

    # Save the results ####
    ## select prps name ####
    if(is.null(prps.group)){
        prps.group <- paste0('prps|knnMnn|', uv.variable)
    }
    ## select output name ####
    if (is.null(output.name))
        output.name <- paste0(uv.variable, '|' , assay.name)

    ## save the results in the SummarizedExperiment object ####
    message(' ')
    printColoredMessage(message = '-- Saving the results.',
                        color = 'magenta',
                        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        if (!'un.supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['un.supervised']] <- list()
        }
        if (!prps.group %in% names(se.obj@metadata[['PRPS']][['un.supervised']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]] <- list()
        }
        if (!'KnnMnn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']] <- list()
        }
        if (!'knn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['mnn']] <- list()
        }
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


