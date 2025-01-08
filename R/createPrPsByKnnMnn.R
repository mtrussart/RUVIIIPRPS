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
#' @param data.input Symbol. Indicates which data should be used an input for finding the k nearest neighbors data. Options
#' include: 'expr' and 'pcs'. If 'pcs' is selected, the first PCs of the data will be used as input. If 'expr' is selected,
#' the data will be selected as input.
#' @param filter.prps.sets Logical. If 'TRUE', the number of PRPS sets across each pair of batches will be filtered if they are
#' higher than the 'max.prps.sets' value. The default is 'TRUE'. The high number of PRPS sets will just increase the
#' computational time for the RUV-III normalization.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS set across each pair of batches.
#' The default is set to 20.
#' @param nb.pcs Numeric. A numeric value indicating the number PCs should be used as data input for finding the k nearest
#' neighbors. The nb.pcs' must be set when the "data.input = PCs". The default is 2.
#' @param center Logical. Indicates whether to scale the data or not. If center is TRUE, then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param k.nn Numeric. The maximum number of nearest neighbors to compute. The default is set 3.
#' @param mnn Numeric. The maximum number of nearest neighbors to compute. The default is set 1.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find the anchors
#' samples across the batches. The default is NULL.
#' @param normalization Symbol. Indicates which normalization methods should be applied before finding the knn. The default
#' is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param apply.cosine.norm Logical. Idicates whether to apply cosine normalization on the data or not. The default is
#' set to 'FALSE'.
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
#' based on "paste0(uv.variable, '||' , assay.name)".
#' @param prps.group Symbol. A symbol specifying the name of the output file. If is 'NULL', the function will select a name
#' based on "paste0('prps_mnn_', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.
#' @return The SummarizedExperiment object that contain all the PPRS data, knn, mnn and plot results in the metadata, or
#' a list of the results.

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom BiocNeighbors findMutualNN KmknnParam
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel SerialParam
#' @importFrom stats dist
#' @importFrom RANN nn2
#' @export

createPrPsByKnnMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        filter.prps.sets = TRUE,
        max.prps.sets = 10,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        clustering.method = 'kmeans',
        nb.clusters = 3,
        k.nn = 2,
        mnn = 1,
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
    printColoredMessage(message = '------------The createPrPsByKnnMnn function starts:',
                        color = 'white',
                        verbose = verbose)

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
    # Assessing and grouping the unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[uv.variable]]
    if (is.numeric(initial.variable)){
        se.obj[[uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_group',
            verbose = verbose
        )
    }
    if (!is.numeric(initial.variable)){
        length.variable <- length(unique(initial.variable))
        if (length.variable == 1){
            stop('To create PRPS, the "uv.variable" must have at least two groups/levels.')
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
    ## KNN ####
    sub.group.sample.size.knn <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = nb.knn + 1
        )
    if (length(sub.group.sample.size.knn) == 0){
        stop(paste0(
            'All subgroups of the unwanted variable have less than ',
            nb.knn + 1,
            ' samples. KNN cannot be found.'))
    } else if (length(sub.group.sample.size.knn) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the unwanted variable have less than ',
                nb.knn + 1,
                ' samples. Then KNN for those sub-groups cannot be created.'),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the sub-groups of the unwanted variable have at least ',
                nb.knn + 1,
                ' samples.'),
            color = 'blue',
            verbose = verbose
        )
    }
    ## MNN ####
    sub.group.sample.size.mnn <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = nb.mnn + 1
    )
    if (length(sub.group.sample.size.mnn) == 0){
        stop(paste0(
            'All subgroups of the unwanted variable have less than ',
            nb.mnn + 1,
            ' samples. MNN cannot be found.'))
    } else if (length(sub.group.sample.size.mnn) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the unwanted variable have less than ',
                nb.mnn + 1,
                ' samples. '),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the subgroups of the unwanted variable have at least, ',
                nb.mnn + 1,
                ' samples.'),
            color = 'blue',
            verbose = verbose
        )
    }

    # Finding common sub-groups ####
    common.sub.groups <- intersect(
        sub.group.sample.size.knn,
        sub.group.sample.size.mnn
        )
    if(length(common.sub.groups) == 0){
        stop('There is not')
    } else if (length(common.sub.groups) == 1){
        stop('')
    } else if (length(common.sub.groups) > 1){
        printColoredMessage()
    }
    uv.variable

    se.obj <- se.obj[ , se.obj[[uv.variable]] %in% common.sub.groups]
    # Finding k nearest neighbor ####
    printColoredMessage(
        message = '-- Finding k nearest neighbor, applying the findKnn function:',
        color = 'magenta',
        verbose = verbose
        )
    se.obj <- findKnn(
        se.obj = se.obj,
        assay.name = assay.name,
        uv.variable = uv.variable,
        data.input = data.input,
        nb.pcs = nb.pcs,
        center = center,
        scale = scale,
        svd.bsparam = svd.bsparam,
        clustering.method = clustering.method,
        nb.clusters = nb.clusters,
        nb.knn = nb.knn,
        hvg = hvg,
        normalization = normalization,
        regress.out.variables = regress.out.variables,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        assess.se.obj = assess.se.obj,
        remove.na = remove.na,
        output.name = output.name,
        prps.group = prps.group,
        save.se.obj = save.se.obj,
        verbose = verbose
    )
    if (is.null(output.name)) {
        output.name.knn <- paste0(uv.variable, '|' , assay.name)
    } else output.name.knn <- output.name
    if (is.null(prps.group)){
        prps.group.mnn <- paste0('prps|knnMnn|', uv.variable)
    } else prps.group.mnn <- prps.group
    all.knn <- se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$knn[[output.name.knn]]
    if (isFALSE(save.se.obj)) {
        se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$knn[[output.name.knn]] <- NULL
    }

    # Find mutual nearest neighbor ####
    printColoredMessage(
        message = '-- Finding mutual nearest neighbors by applying the findMnn function:',
        color = 'magenta',
        verbose = verbose
    )
    se.obj <- findMnn(
        se.obj = se.obj,
        assay.name = assay.name,
        uv.variable = uv.variable,
        clustering.method = clustering.method,
        nb.clusters = nb.clusters,
        nb.mnn = nb.mnn,
        hvg = hvg,
        normalization = normalization,
        apply.cosine.norm = apply.cosine.norm,
        regress.out.variables = regress.out.variables,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        assess.se.obj = assess.se.obj,
        mnn.bpparam = mnn.bpparam,
        mnn.nbparam = mnn.nbparam,
        remove.na = remove.na,
        output.name = output.name,
        prps.group = prps.group,
        save.se.obj = save.se.obj,
        verbose = verbose
    )
    if (is.null(output.name)) {
        output.name.mnn <- paste0(uv.variable, '|' , assay.name)
    } else output.name.mnn <- output.name
    if (is.null(prps.group)){
        prps.group.mnn <- paste0('prps|knnMnn|', uv.variable)
    } else prps.group.mnn <- prps.group
    all.mnn <- se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$mnn[[output.name.mnn]]
    if (isFALSE(save.se.obj)) {
        se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$mnn[[output.name.mnn]] <- NULL
    }

    # Finding similar sample sets ####
    printColoredMessage(
        message = '-- Finding all possible similar samples across batches:',
        color = 'magenta',
        verbose = verbose
        )
    ## find the knn for each mnn set ####
    printColoredMessage(
        message = '- Finding similar samples using the KNN and MNN data:',
        color = 'orange',
        verbose = verbose
    )
    printColoredMessage(
        message = '* matching the MNN sets with the corresponding KNN sets:',
        color = 'blue',
        verbose = verbose
    )
    mnn.sets <- NULL
    all.prps.sets <- lapply(
        1:nrow(all.mnn),
        function(x) {
            # ps set 1
            ps.set.1 <- all.knn[ , c(1:c(nb.knn + 1))] == all.mnn$sample.no.1[x]
            ps.set.1 <- all.knn[rowSums(ps.set.1) > 0 , ]
            ps.set.1$mnn.sets <- paste0(sort(c(all.mnn[x , 3], all.mnn[x , 4])), collapse = '_')
            ps.set.1$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '_')
            if (nrow(ps.set.1) > 1) {
                ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
            }
            # ps set 2
            ps.set.2 <- all.knn[, c(1:c(nb.knn + 1))] == all.mnn$sample.no.2[x]
            ps.set.2 <- all.knn[rowSums(ps.set.2) > 0 , ]
            ps.set.2$mnn.sets <- paste0(sort(c(all.mnn[x , 3], all.mnn[x , 4])), collapse = '_')
            ps.set.2$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '_')
            if (nrow(ps.set.2) > 1) {
                ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
            }
            prps.set <- rbind(ps.set.1, ps.set.2)
            prps.set
        })
    all.prps.sets <- do.call(rbind, all.prps.sets)

    if (is.null(all.prps.sets)) {
        stop('PRPS cannot be created. You may want to increase the value of the mnn.')
    }
    ### sanity check ####
    if(nrow(all.prps.sets) == 2*nrow(all.mnn)){
        printColoredMessage(
            message = paste0('* The nrow of the matched MNN and KNN is ', nrow(all.prps.sets), '.'),
            color = 'blue',
            verbose = verbose
        )
    } else stop('For individual MNN set, the corresponding KNN sets cannot be found. Check the the input.')

    ## add the average of the knn sets for each PRPS set and then rank them ####
    printColoredMessage(
        message = '* Average the knn sets for each MNN set and then rank them:',
        color = 'blue',
        verbose = verbose
    )
    aver.mnn.sets <- NULL
    all.prps.sets$aver.mnn.sets <- unlist(lapply(
        seq(1, nrow(all.prps.sets), 2),
        function(x)
            rep(mean(all.prps.sets$aver.dist[x:(x + 1)]), 2))
        )
    set.seed(2233)
    all.prps.sets$rank.aver.mnn.sets <- rank(x = all.prps.sets$aver.mnn.sets, ties.method = 'random')

    ## filter PRPS sets ####
    if (isTRUE(filter.prps.sets)) {
        printColoredMessage(
            message = '- Filter the PRPS sets across each pair of batches:',
            color = 'orange',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0( '- The maximum number of PRPS sets per batch is set to ', max.prps.sets, '.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- The PRPS sets will be filtered based on the distances between each knn sets.'),
            color = 'blue',
            verbose = verbose
        )
        all.prps.sets <- lapply(
            unique(all.prps.sets$mnn.sets.data),
            function(x) {
                temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                if (length(unique(temp.prps.set$mnn.sets)) >= max.prps.sets) {
                    printColoredMessage(
                        message = paste0(
                            '* The number of PRPS sets across the batches "', x, '" is ',
                            length(unique(temp.prps.set$mnn.sets)), '.' ),
                        color = 'blue',
                        verbose = verbose
                    )
                    temp.prps.set <- arrange(temp.prps.set, aver.mnn.sets, mnn.sets)
                    printColoredMessage(
                        message = paste0('* ', length(unique(temp.prps.set$mnn.sets)) - max.prps.sets, ' PRPS sets are removed.' ),
                        color = 'blue',
                        verbose = verbose
                    )
                    temp.prps.set <- temp.prps.set[1:c(2 * max.prps.sets) , ]
                } else if (length(unique(temp.prps.set$mnn.sets))  < max.prps.sets) {
                    printColoredMessage(
                        message = paste0(
                            '* The number of PRPS sets across the batches "', x, '" is ',
                            length(unique(temp.prps.set$mnn.sets)) , '.' ),
                        color = 'blue',
                        verbose = verbose
                    )
                    temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                }
                return(temp.prps.set)
            })
        all.prps.sets <- do.call(rbind, all.prps.sets)
    }
    printColoredMessage(
        message = paste0('- ', length(unique(all.prps.sets$mnn.sets)), ' PRPS stes are found in total.'),
        color = 'blue',
        verbose = verbose
    )

    # Create PRPS data ####
    printColoredMessage(
        message = '-- Create PRPS expression data:',
        color = 'magenta',
        verbose = verbose
        )
    ## apply log ####
    printColoredMessage(
        message = '- Data log transformation before creating the PRPS expression data:',
        color = 'blue',
        verbose = verbose
        )
    if (isTRUE(apply.log) & !is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0(
                '- Applying log2 on the "',
                assay.name, '" + ',
                pseudo.count,
                ' (pseudo.count)  data.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0(
                'Applying log2 on the "',
                assay.name,
                '" data.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The "',
                assay.name,
                '" data will be used without any log transformation.' ),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    printColoredMessage(
        message = '- Aeveraging samples to create pseudo samples:',
        color = 'blue',
        verbose = verbose
        )
    prps.data <- lapply(
        unique(all.prps.sets$mnn.sets),
        function(x) {
            temp.prps <- all.prps.sets[all.prps.sets$mnn.sets == x, ]
            index.a <- unlist(unname(temp.prps[1, grep('overal', colnames(temp.prps))]))
            index.b <- unlist(unname(temp.prps[2, grep('overal', colnames(temp.prps))]))
            prps.a <- rowMeans(expr.data[, index.a])
            prps.b <- rowMeans(expr.data[, index.b])
            prps <- cbind(prps.a, prps.b)
            colnames(prps) <- paste(uv.variable, temp.prps$mnn.sets, sep = '_')
            return(prps)
        })
    prps.data <- do.call(cbind, prps.data)
    ## sanity check ####
    if (!sum(table(colnames(prps.data)) == 2) == ncol(prps.data) / 2) {
        stop('There someting wrong with PRPS sets.')
    }

    se.obj[[uv.variable]] <- initial.variable

    # Save the results ####
    ## select output name ####
    output.name <- paste0(uv.variable, '|', 'mnn', '|', assay.name)
    if (is.null(prps.group)) {
        prps.group <- paste0('prps|knnMnn|', uv.variable)
    }
    printColoredMessage(message = '-- Saving the PRPS data',
                        color = 'magenta',
                        verbose = verbose)
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
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- prps.data

        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.data = prps.data)
    }
}
