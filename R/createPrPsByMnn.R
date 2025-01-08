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
        main.uv.variable,
        other.uv.variables = NULL,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        nb.other.uv.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        cont.cor.coef = c(.9, .9),
        cat.cor.coef = c(.9, .9),
        min.sample.for.ps = 3,
        min.mnn = 3,
        min.batch.number = 'all',
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
    if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single assay name.')
    } else if (is.null(assay.name)) {
        stop('The "assay.name" cannot be empty.')
    } else if(is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty')
    } else if (length(main.uv.variable) > 1) {
        stop('The function can only take a single categorical uv variable for the "main.uv.variable" argument.')
    } else if (is.null(main.uv.variable)) {
        stop('The main.uv.variable cannot be empty.')
    } else if (!class(se.obj[[main.uv.variable]]) %in% c('character', 'factor')) {
        stop('The "main.uv.variable" must be a categorical source of unwanted variation, e.g., platform effects')
    } else if (length(unique(se.obj[[main.uv.variable]])) == 1) {
        stop('The "main.uv.variable" must have at least two levels or groups.')
    } else if (min.sample.for.ps <= 1) {
        stop('The minimum value for the "min.sample.for.ps" is 2.')
    } else if(main.uv.variable %in% other.uv.variables){
        stop('The "main.uv.variable" must not be in the "other.uv.variables".')
    } else if(!is.logical(apply.log)){
        stop('The "apply.log" must be logical.')
    } else if(!is.logical(check.prps.connectedness)){
        stop('The "check.prps.connectedness" must be logical.')
    } else if(!is.logical(assess.se.obj)){
        stop('The "assess.se.obj" must be logical.')
    } else if (pseudo.count < 0){
        stop('The value for "pseudo.count" can not be negative.')
    } else if (max(cat.cor.coef) > 1 | max(cont.cor.coef) > 1){
        stop('The maximum value for "cat.cor.coef" or "cont.cor.coef" cannot be more than 1.')
    } else if (max(cat.cor.coef) < 0 | max(cont.cor.coef) < 0){
        stop('The maximum value for "cat.cor.coef" or "cont.cor.coef" cannot be more negative.')
    } else if(!is.logical(plot.prps.map)){
        stop('The "plot.prps.map" must be logical.')
    } else if(!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical.')
    } else if(is.logical(output.name)){
        stop('The "output.name" must be a character or NULL.')
    } else if(is.logical(prps.group)){
        stop('The "prps.group" must be a character or NULL.')
    }
    # Assessing SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, other.uv.variables, regress.out.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Applying data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
    )
    if (!is.null(normalization) & is.null(regress.out.variables)) {
        printColoredMessage(
            message = paste0(
                '- applying the ',
                normalization,
                ' normalization on the data before finding MNN.'),
            color = 'blue',
            verbose = verbose
        )
        norm.data <- applyOtherNormalizations(
            se.obj = se.obj,
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
                ' normalization and then regressing out the ',
                paste0(regress.out.variables, collapse = '&'),
                ' variable(s) from the data before finding MNN.'),
            color = 'blue',
            verbose = verbose
        )
        ### normalization ####
        norm.data <- applyOtherNormalizations(
            se.obj = se.obj,
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
        sample.info <- as.data.frame(colData(se.obj))
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
                    '- applying log2 transformation and then regressing out the ',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variable(s) from the data before finding MNN. '),
                color = 'blue',
                verbose = verbose
            )
            if(!is.null(pseudo.count)){
                norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- regressing out the ',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variable(s) from the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            norm.data <- assay(x = se.obj, i = assay)
        }
        ### regression ####
        sample.info <- as.data.frame(colData(se.obj))
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
                    '- apply the log2 transformation on the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            if(!is.null(pseudo.count)){
                norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- no library size normalization and transformation is applied on the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            norm.data <- assay(x = se.obj, i = assay.name)
        }
    }
    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[main.uv.variable]]
    if(is.numeric(initial.variable)){
        se.obj[[main.uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = main.uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_group',
            verbose = verbose
        )
    }
    if(!is.numeric(initial.variable)){
        length.variable <- length(unique(initial.variable))
        if( length.variable == 1){
            stop('To create MNN, the "main.uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    main.uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[main.uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[main.uv.variable]] <- factor(x = se.obj[[main.uv.variable]])
        }
    }

    # Assessing and grouping the other unwanted variable ####
    if(!is.null(other.uv.variables)){
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.uv.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            cont.cor.coef = cont.cor.coef,
            cat.cor.coef = cat.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
        )
        all.uv.groups <- data.frame(
            main.uv = se.obj[[main.uv.variable]],
            other.uv = homo.uv.groups
        )
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = max(min.sample.for.ps, min.mnn)
                )
            })
        if(min.batch.number == 'all') {
            covered.batches.length <- sapply(
                covered.batches,
                function(x){
                    length(x)
                })
            if (sum(covered.batches.length == length(unique(all.uv.groups$main.uv))) == 0 ){
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
        if(is.numeric(min.batch.number)) {
            if (sum(covered.batches.length >= min.batch.number) > 0){
                printColoredMessage(
                    message = paste0(
                        '- At least ',
                        min.batch.number,
                        ' sub-groups of the variable ',
                        main.uv.variable,
                        ' have at least ',
                        max(min.sample.for.ps, min.mnn),
                        ' samples.'),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
    }
    if(is.null(other.uv.variables)){
        # Check the sample size of each group in the variable ####
        subgroups.size <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = max(min.sample.for.ps, min.mnn)
        )
        if(min.batch.number == 'all') {
            if (length(subgroups.size) != length(unique(se.obj[[main.uv.variable]])) ){
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
        if(is.numeric(min.batch.number)) {
            if (length(subgroups.size) >= min.batch.number){
                printColoredMessage(
                    message = paste0(
                        '- At least ',
                        min.batch.number,
                        ' sub-groups of the variable ',
                        main.uv.variable,
                        ' have at least ',
                        max(min.sample.for.ps, min.mnn),
                        ' samples.'),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
    }

    # Create PRPS data ####
    ## find mutual nearest neighbor ####
    printColoredMessage(
        message = '-- Creating PRPS data across all pairs of the subgroups:',
        color = 'magenta',
        verbose = verbose
        )
    all.possible.batches <- lapply(
        unique(all.uv.groups$other.uv),
        function(x){
            combn(
                x = findRepeatingPatterns(
                    vec = all.uv.groups[all.uv.groups$other.uv == x, ]$main.uv,
                    n.repeat = max(min.sample.for.ps, min.mnn)),
                m = 2
            )
        })
    names(all.possible.batches) <- unique(all.uv.groups$other.uv)
    prps.data <- assay(x = se.obj, i = assay.name)
    sample.annotation <- as.data.frame(colData(x = se.obj))
    all.prps.data <- lapply(
        1:length(all.possible.batches),
        function(x){
            sample.annotation.all <- sample.annotation[all.uv.groups$other.uv == names(all.possible.batches)[x] , ]
            pairs.batch <- all.possible.batches[[x]]
            sub.prps.data <- lapply(
                1:ncol(pairs.batch),
                function(y){
                    printColoredMessage(
                        message = paste0(
                            '- Creating PRPS data between the "',
                            pairs.batch[1 , y],
                            '" and "' ,
                            pairs.batch[2 , y],
                            '" subgroups:'),
                        color = 'orange',
                        verbose = verbose
                    )
                    ## sample annotation
                    sample.annot.a <- sample.annotation.all[sample.annotation.all[[main.uv.variable]] == pairs.batch[1 , y] , , drop = FALSE]
                    sample.annot.b <- sample.annotation.all[sample.annotation.all[[main.uv.variable]] == pairs.batch[2 , y] , , drop = FALSE]
                    if(is.null(hvg)){
                        printColoredMessage(
                            message = '* using all genes.',
                            color = 'blue',
                            verbose = verbose
                        )
                        data.a <- norm.data[ , row.names(sample.annot.a)]
                        data.b <- norm.data[ , row.names(sample.annot.b)]
                    }
                    if(!is.null(hvg)){
                        printColoredMessage(
                            message = '* using the specified highly variable genes.',
                            color = 'blue',
                            verbose = verbose
                        )
                        data.a <- norm.data[hvg , row.names(sample.annot.a)]
                        data.b <- norm.data[hvg , row.names(sample.annot.b)]
                    }
                    if(isTRUE(apply.cosine.norm)){
                        printColoredMessage(
                            message = '- Applying cosine normalization on the data:',
                            color = 'blue',
                            verbose = verbose
                        )
                        data.a <- cosineNorm(x = data.a, mode = 'matrix')
                        data.b <- cosineNorm(x = data.b, mode = 'matrix')
                    }
                    ## find mnn for data b ####
                    printColoredMessage(
                        message = paste0(
                            '* finding ',
                            min.sample.for.ps,
                            ' number of near neighbours for each point of the "',
                            pairs.batch[2 , y],
                            '" using "' ,
                            pairs.batch[1 , y],
                            '" subgroups:'),
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.a.b <- RANN::nn2(
                        data = t(data.a),
                        query = t(data.b),
                        k = min.sample.for.ps
                    )
                    ### index
                    printColoredMessage(
                        message = '** obtaining the knn indexs.',
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.a.b.index <- knn.data.a.b$nn.idx
                    colnames(knn.data.a.b.index) <- paste0('knn', seq(min.sample.for.ps))
                    row.names(knn.data.a.b.index) <- c(1:ncol(data.b))
                    ### distance
                    printColoredMessage(
                        message = '** obtaining the distances.',
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.a.b.distance <- as.data.frame(knn.data.a.b$nn.dists)
                    colnames(knn.data.a.b.distance) <- paste0('knn', seq(min.sample.for.ps))
                    row.names(knn.data.a.b.distance) <- c(1:ncol(data.b))
                    knn.data.a.b.distance$aver.dist <- rowMeans(knn.data.a.b.distance)

                    ## find knn for data a ####
                    printColoredMessage(
                        message = paste0(
                            '* finding ',
                            min.sample.for.ps,
                            ' number of near neighbours for each point of the "',
                            pairs.batch[1 , y],
                            '" using "' ,
                            pairs.batch[2 , y],
                            '" subgroups:'),
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.b.a <- RANN::nn2(
                        data = t(data.b),
                        query = t(data.a),
                        k = min.sample.for.ps
                    )
                    printColoredMessage(
                        message = '** obtaining the knn indexs.',
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.b.a.index <- knn.data.b.a$nn.idx
                    colnames(knn.data.b.a.index) <- paste0('knn', seq(min.sample.for.ps))
                    row.names(knn.data.b.a.index) <- c(1:ncol(data.a))
                    ### distance
                    printColoredMessage(
                        message = '** obtaining the distances.',
                        color = 'blue',
                        verbose = verbose
                    )
                    knn.data.b.a.distance <- as.data.frame(knn.data.b.a$nn.dists)
                    colnames(knn.data.b.a.distance) <- paste0('knn', seq(min.sample.for.ps))
                    row.names(knn.data.b.a.distance) <- c(1:ncol(data.a))
                    knn.data.b.a.distance$aver.dist <- rowMeans(knn.data.b.a.distance)

                    ## final mnn ####
                    printColoredMessage(
                        message = paste0(
                            '* finding MNN between the "',
                            pairs.batch[1 , y],
                            '" using "' ,
                            pairs.batch[2 , y],
                            '" subgroups:'),
                        color = 'blue',
                        verbose = verbose
                    )
                    all.mnn <- BiocNeighbors::findMutualNN(
                        data1 = t(data.a),
                        data2 = t(data.b),
                        k1 = min.mnn,
                        BNPARAM = mnn.nbparam,
                        BPPARAM = mnn.bpparam
                    )
                    printColoredMessage(
                        message = paste0(
                            '** ',
                            length(all.mnn$first),
                            ' MNN are found.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    ## create prps index ####
                    printColoredMessage(
                        message = paste0('* creating PRPS indexs and score.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    all.prps.index <- lapply(
                        1:length(all.mnn$first),
                        function(z){
                            sample.to.ave.b <- knn.data.b.a.index[all.mnn$first[z] , ]
                            sample.to.ave.a <- knn.data.a.b.index[all.mnn$second[z] , ]
                            prps.index <- data.frame(
                                set.a = sample.to.ave.a,
                                set.b = sample.to.ave.b
                            )
                            samples.ids <- data.frame(
                                set.a = row.names(sample.annot.a)[sample.to.ave.a],
                                set.b = row.names(sample.annot.b)[sample.to.ave.b]
                            )
                            avre.dist.a <- knn.data.b.a.distance[all.mnn$first[z] , 'aver.dist' ]
                            avre.dist.b <- knn.data.a.b.distance[all.mnn$second[z] , 'aver.dist' ]
                            list(
                                prps.index = prps.index,
                                aver.dist = c(avre.dist.a+avre.dist.b)/2,
                                samples.ids = samples.ids
                            )
                        })
                    names(all.prps.index) <- paste0('prps.set', 1:length(all.mnn$first))
                    if(isTRUE(filter.prps.sets)){
                        printColoredMessage(
                            message = '* filtering the number of PRPS sets.',
                            color = 'blue',
                            verbose = verbose
                        )
                        if(length(all.prps.index)  >=  max.prps.sets){
                            aver.dists <- sapply(
                                1:length(all.prps.index),
                                function(p) all.prps.index[[p]]$aver.dist
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
                        message = '* creating the PRPS data:',
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- prps.data[ , colnames(data.a)]
                    data.b <- prps.data[ , colnames(data.b)]
                    if(isTRUE(apply.log)){
                        printColoredMessage(
                            message = '** applying log transformation on the data before creating PRPS. ',
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
                        function(a){
                            prps.set <- all.prps.index[[a]]
                            ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                            ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                            prps <- data.frame(ps.a, ps.b)
                            prps
                        })
                    all.prps.data <- do.call(cbind, all.prps.data)
                    colnames(all.prps.data) <- paste0(
                        pairs.batch[1 , y],
                        pairs.batch[2 , y],
                        rep(1:c(ncol(all.prps.data)/2), each = 2)
                    )
                    all.prps.sample.annot <- lapply(
                        1:length(all.prps.index),
                        function(f){
                            prps.set <- all.prps.index[[f]]$samples.ids
                        })
                    list(
                        all.prps.data = all.prps.data,
                        all.prps.sample.annot = all.prps.sample.annot
                    )
                })
            sub.prps.epxr.data <- lapply(
                1:length(sub.prps.data),
                function(t) sub.prps.data[[t]]$all.prps.data)
            sub.prps.epxr.data <- do.call(cbind, sub.prps.epxr.data)
            sub.prps.info.data <- lapply(
                1:length(sub.prps.data),
                function(t) sub.prps.data[[t]]$all.prps.sample.annot)
            return(list(
                sub.prps.epxr.data = sub.prps.epxr.data,
                sub.prps.info.data = sub.prps.info.data)
                )
        })

    all.prps.expr.data <- lapply(
        1:length(all.prps.data),
        function(x){
            all.prps.data[[x]]$sub.prps.epxr.data
        })
    all.prps.expr.data <- do.call(cbind, all.prps.expr.data)
    se.obj[[ main.uv.variable]] <- initial.variable
    sample.annotation <- as.data.frame(colData(x = se.obj))

    # Plot PRPS map ####
    all.prps.plot <- lapply(
        1:length(all.prps.data),
        function(x){
            prps.index <- all.prps.data[[x]]$all.prps.sample.annot
            expr.data <- lapply(
                1:length(prps.index),
                function(y){
                    data.frame(
                        group1 = sample.annotation[prps.index[[y]]$set.a , ][[ main.uv.variable]],
                        group2 =  sample.annotation[prps.index[[y]]$set.b , ][[ main.uv.variable]],
                        set.name = rep(
                            paste0('prps.set.', y),
                            length(sample.annotation[prps.index[[y]]$set.a , ][[ main.uv.variable]])),
                        group.name = rep(
                            paste0(pairs.batch[1 , x],'||' ,  pairs.batch[2 , x]),
                            length(sample.annotation[prps.index[[y]]$set.a , ][[ main.uv.variable]]))
                        )
                })
            expr.data <- do.call(rbind,expr.data )
        })
    all.prps.plot <- do.call(rbind, all.prps.plot) %>%
        pivot_longer(-c(set.name, group.name), values_to = 'var', names_to = 'rep')
    all.prps.plot$new.g <- paste0(all.prps.plot$group.name, all.prps.plot$set.name)
    ggplot(all.prps.plot, aes(x = new.g, y = var)) +
        geom_point() +
        facet_grid(~new.g, scales = 'free', space = 'free') +
        scale_x_discrete(expand = c(0, 0.5)) +
        # xlab(main.uv.variable) +
        ylab('Homogeneous groups') +
        # xlim(c(
        #     min(se.obj[[main.uv.variable]]),
        #     max(se.obj[[main.uv.variable]])
        # )) +
        # geom_hline(yintercept = c(
        #     min(se.obj[[main.uv.variable]]),
        #     max(se.obj[[main.uv.variable]])), color = 'gray70') +
        theme_bw() +
        theme(
            legend.key = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 18),
            strip.text.y = element_text(size = 0)
        )


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

