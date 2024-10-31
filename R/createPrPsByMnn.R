#' Create PRPS sets using k and mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses the k and mutual nearest neighbors approaches to create PRPS in the RNA-seq data. This function
#' can be used in situation that the biological variation are entirely unknown. The function applies the 'findKnn' function
#' to find similar samples per batch and then average them to create pseudo-samples. Then, function uses the 'findMnn' to
#' match up pseudo samples across batches to create pseudo-replicates.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol indicating the name of the assay in the SummarizedExperiment object. This assay will
#' be used to first find k nearest neighbors and them mutual nearest neighbors data. This data must the one that will be
#' used for the RUV-III normalization.
#' @param uv.variable Symbol. Indicates the name of a column in the sample annotation of the SummarizedExperiment object.
#' The 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
#' divided into 'nb.clusters' groups using the 'clustering.method'.
#' @param filter.prps.sets Logical. If 'TRUE', the number of PRPS sets across each pair of batches will be filtered if they are
#' higher than the 'max.prps.sets' value. The default is 'TRUE'. The high number of PRPS sets will just increase the
#' computational time for the RUV-III normalization.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS set across each pair of batches.
#' The default is set to 10.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param mnn Numeric. The maximum number of nearest neighbors to compute. The default is set 3.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find the anchors
#' samples across the batches. The default is NULL.
#' @param normalization Symbol. Indicates which normalization methods should be applied before finding the knn. The default
#' is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param apply.cosine.norm TTT
#' @param regress.out.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed. The default
#' is SerialParam(). We refer to the 'findMutualNN' function from the BiocNeighbors R package for more details.
#' @param mnn.nbparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is KmknnParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#' or to output the result as list. By default it is set to TRUE.
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '|', 'anchor', '|', assay.name))".
#' @param prps.group Symbol. A symbol specifying the name of the output file. If is 'NULL', the function will select a name
#' based on "paste0('prps_mnn_', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return The SummarizedExperiment object that contain all the PPRS data, knn, mnn and plot results in the metadata, or
#' a list of the results.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN
#' @importFrom batchelor cosineNorm
#' @importFrom RANN nn2
#' @export

createPrPsByMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        mnn = 3,
        hvg = NULL,
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsByMnn function starts:',
                        color = 'white',
                        verbose = verbose
                        )
    # Assess SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # Assess and group the unwanted variable ####
    printColoredMessage(
        message = '- Assess and group the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[uv.variable]]
    if(is.numeric(initial.variable)){
        se.obj[[uv.variable]] <- groupContiunousVariable(
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
                    '- The "', uv.variable, '" is a categorical variable with ',
                    length(unique(se.obj[[uv.variable]])), ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
        }
    }

    ## Sanity check on variable groups ####
    subgroups.size <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = mnn
    )
    if (length(subgroups.size) != length(unique(se.obj[[uv.variable]])) ){
        stop(paste0(
            'Some sub-groups of the variable "', uv.variable, '" have less than ', mnn,
            ' (mnn) samples. Then, knn cannot be found.'))
    } else {
        printColoredMessage(
            message = paste0('- all the subgroups of the unwanted variable have at least ', mnn, ' samples.'),
            color = 'blue',
            verbose = verbose
        )
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
            ## normalization ####
            if (!is.null(normalization) & is.null(regress.out.variables)) {
                printColoredMessage(
                    message = paste0('- apply ', normalization, ' on the samples from the "', x, '" subgroup.'),
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
                        '- apply ', normalization, ' on the samples from "', x,
                        '" subgroup and then regressing out ', paste0(regress.out.variables, collapse = '&'),
                        ' from the data.'),
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
                            '- apply log and then regress out ', paste0(regress.out.variables, collapse = '&'),
                            ' on the ',x, '" group from the data.'),
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
                norm.data

            }
            ## log transformation ####
            if (is.null(normalization) & is.null(regress.out.variables)) {
                if(isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- apply the log2 within the samples from "', x, '" subgroup.'),
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
                            '- no library size normalization and transformation is applied on samples from ',
                            x, '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], i = assay.name)
                }
            }
            return(norm.data)
        })
    names(all.norm.data) <- groups

    # Create PRPS data ####
    ## find mutual nearest neighbor ####
    printColoredMessage(
        message = '-- Create PRPS data across all pairs of subgroups:',
        color = 'magenta',
        verbose = verbose
        )
    pairs.batch <- combn(
        x = levels(se.obj[[uv.variable]]),
        m = 2
        )
    prps.data <- assay(x = se.obj, i = assay.name)
    all.prps.data <- lapply(
        1:ncol(pairs.batch),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Create PRPS data between "', pairs.batch[1 , x],
                    '" and "' , pairs.batch[2 , x], '" subgroups:'),
                color = 'orange',
                verbose = verbose
                )
            data.a <- all.norm.data[[pairs.batch[1 , x]]][hvg , ]
            data.b <- all.norm.data[[pairs.batch[2 , x]]][hvg , ]
            if(isTRUE(apply.cosine.norm)){
                printColoredMessage(
                    message = '- Apply cosine normalization',
                    color = 'blue',
                    verbose = verbose
                    )
                data.a <- cosineNorm(x = data.a, mode = 'matrix')
                data.b <- cosineNorm(x = data.b, mode = 'matrix')
            }
            ## find mnn for data b ####
            printColoredMessage(
                message = paste0(
                    '* find knn between the "', pairs.batch[1 , x],
                    '" and "' , pairs.batch[2 , x], '" subgroups:'),
                color = 'blue',
                verbose = verbose
                )
            mnn.data.a.b <- RANN::nn2(
                data = t(data.a),
                query = t(data.b),
                k = mnn
                )
            ### index
            printColoredMessage(
                message = '** obtain knn indexs.',
                color = 'blue',
                verbose = verbose
                )
            mnn.data.a.b.index <- mnn.data.a.b$nn.idx
            colnames(mnn.data.a.b.index) <- paste0('mnn', seq(mnn))
            row.names(mnn.data.a.b.index) <- c(1:ncol(data.b))
            ### distance
            printColoredMessage(
                message = '** obtain distances.',
                color = 'blue',
                verbose = verbose
                )
            mnn.data.a.b.distance <- as.data.frame(mnn.data.a.b$nn.dists)
            colnames(mnn.data.a.b.distance) <- paste0('mnn', seq(mnn))
            row.names(mnn.data.a.b.distance) <- c(1:ncol(data.b))
            mnn.data.a.b.distance$aver.dist <- rowMeans(mnn.data.a.b.distance)

            ## find knn for data a ####
            printColoredMessage(
                message = paste0(
                    '* find knn between the "', pairs.batch[2 , x],
                    '" using "' , pairs.batch[1 , x], '" subgroups:'),
                color = 'blue',
                verbose = verbose
            )
            mnn.data.b.a <- RANN::nn2(
                data = t(data.b),
                query = t(data.a),
                k = mnn
            )
            printColoredMessage(
                message = '** obtain knn indexs.',
                color = 'blue',
                verbose = verbose
            )
            mnn.data.b.a.index <- mnn.data.b.a$nn.idx
            colnames(mnn.data.b.a.index) <- paste0('mnn', seq(mnn))
            row.names(mnn.data.b.a.index) <- c(1:ncol(data.a))
            ### distance
            printColoredMessage(
                message = '** obtain distances.',
                color = 'blue',
                verbose = verbose
                )
            mnn.data.b.a.distance <- as.data.frame(mnn.data.b.a$nn.dists)
            colnames(mnn.data.b.a.distance) <- paste0('mnn', seq(mnn))
            row.names(mnn.data.b.a.distance) <- c(1:ncol(data.a))
            mnn.data.b.a.distance$aver.dist <- rowMeans(mnn.data.b.a.distance)

            ## final mnn ####
            printColoredMessage(
                message = paste0(
                    '* find MNN between the "', pairs.batch[1 , x],
                    '" using "' , pairs.batch[2 , x], '" subgroups:'),
                color = 'blue',
                verbose = verbose
                )
            all.mnn <- BiocNeighbors::findMutualNN(
                data1 = t(data.a),
                data2 = t(data.b),
                k1 = mnn,
                BNPARAM = mnn.nbparam,
                BPPARAM = mnn.bpparam
                )
            printColoredMessage(
                message = paste0('** ', length(all.mnn$first), ' MNN are found.'),
                color = 'blue',
                verbose = verbose
                )

            ## create prps index ####
            printColoredMessage(
                message = paste0('* create PRPS indexs and score.'),
                color = 'blue',
                verbose = verbose
                )
            all.prps.index <- lapply(
                1:length(all.mnn$first),
                function(x){
                    set.a <- mnn.data.b.a.index[all.mnn$first[x] , ]
                    set.b <- mnn.data.a.b.index[all.mnn$second[x] , ]
                    prps.index <- data.frame(set.a = set.b, set.b = set.a)
                    avre.dist.a <- mnn.data.b.a.distance[all.mnn$first[x] , 'aver.dist' ]
                    avre.dist.b <- mnn.data.a.b.distance[all.mnn$second[x] , 'aver.dist' ]
                    list(prps.index = prps.index, aver.dist = c(avre.dist.a+avre.dist.b)/2)
                })
            if(isTRUE(filter.prps.sets)){
                printColoredMessage(
                    message = '* filter the number of PRPS sets.',
                    color = 'blue',
                    verbose = verbose
                )
                if(length(all.prps.index)  >=  max.prps.sets){
                    aver.dists <- sapply(
                        1:length(all.prps.index),
                        function(x) all.prps.index[[x]]$aver.dist
                    )
                    names(aver.dists) <- 1:length(all.prps.index)
                    aver.dists <- aver.dists[order(aver.dists, decreasing = FALSE)]
                    select.prps <- names(aver.dists)[1:max.prps.sets]
                    all.prps.index <- all.prps.index[as.numeric(select.prps)]
                    printColoredMessage(
                        message = paste0('** ', length(all.prps.index), ' PRPS set are kept.'),
                        color = 'blue',
                        verbose = verbose
                    )
                }
            }
            ## prps data
            # Find PRPS sets ####
            printColoredMessage(
                message = '* crete PRPS data:',
                color = 'blue',
                verbose = verbose
                )
            data.a <- prps.data[ , colnames(data.a)]
            data.b <- prps.data[ , colnames(data.b)]
            if(isTRUE(apply.log)){
                printColoredMessage(
                    message = '** apply log on the data before creating PRPS. ',
                    color = 'blue',
                    verbose = verbose
                )
                data.a <- log2(data.a + pseudo.count)
                data.b <- log2(data.b + pseudo.count)
            } else{
                printColoredMessage(
                    message = '** the data will be used without any transformation. ',
                    color = 'blue',
                    verbose = verbose
                )
            }
            all.prps.data <- lapply(
                1:length(all.prps.index),
                function(x){
                    prps.set <- all.prps.index[[x]]
                    ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                    ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                    prps <- data.frame(ps.a, ps.b)
                    prps
                })
            all.prps.data <- do.call(cbind, all.prps.data)
            colnames(all.prps.data) <- paste0(
                pairs.batch[1 , x],
                pairs.batch[2 , x],
                rep(1:c(ncol(all.prps.data)/2), each = 2)
                )
            ### sanity check ####
            all.prps.data
        })
    all.prps.data <- do.call(cbind, all.prps.data)
    se.obj[[uv.variable]] <- initial.variable

    # Save the results ####
    ## select output name ####
    if(is.null(output.name))
        output.name <- paste0(uv.variable, '|', 'mnn', '|', assay.name)
    if (is.null(prps.group))
        prps.group <- paste0('prps|mnn|', uv.variable)

    printColoredMessage(
        message = '-- Save the PRPS data',
        color = 'magenta',
        verbose = verbose
        )
    ## save the PRPS data in the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = 'Save all the PRPS data into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        if (!'un.supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['un.supervised']] <- list()
        }
        if (!prps.group %in% names(se.obj@metadata[['PRPS']][['un.supervised']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]] <- list()
        }
        if (!'prps.data' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- all.prps.data

        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.data = all.prps.data)
    }
}

