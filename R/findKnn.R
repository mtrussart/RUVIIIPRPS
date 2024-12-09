#' Find k-nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function finds k nearest neighbors samples in RNA-seq data. The k nearest neighbors will be used to create pseudo
#' sample within individual batches.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol for the selection of the name of the assay in the SummarizedExperiment object to be
#' used to find k nearest neighbors.
#' @param uv.variable Symbol. A symbol that indicates the name of the column in the sample annotation of the
#' SummarizedExperiment object. The 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous
#' variable, this will be divided into 'nb.clusters' groups using the 'clustering.method'.
#' @param data.input Symbol. A symbol that indicates which data format should be used as input for finding the k nearest
#' neighbors data. Options include: 'expr' and 'pcs'. If 'pcs' is selected, the first 'nb.pcs' of PCs of the data will be
#' used as input. If 'expr' is selected, the expression data will be used as input. The default is set to 'expr'.
#' @param nb.pcs Numeric. A numeric valuse that indicates the number PCs should be used as data input for finding the k
#' nearest neighbors. The nb.pcs' must be set when the "data.input = PCs". The default is set to 2.
#' @param center Logical. Indicates whether to scale the data or not. If center is TRUE, then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE
#' @param svd.bsparam Symbol. A BiocParallelParam object specifying how parallelization should be performed. The default
#' is bsparam(). We refer to the 'runSVD' function from the BiocSingular R package.
#' @param clustering.method Symbol. A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param k.nn Numeric. A numeric number that indicates the maximum number of nearest neighbors to compute. The default
#' is set 3.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to select the input
#' data. The default is NULL.
#' @param normalization Symbol. A symbol that indicates which normalization method should be applied  on the dara before
#' finding the knn. The default is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param regress.out.variables Symbol. A symbol or symbols that indicates the column name(s) in the sample annotation in
#' the SummarizedExperiment. These variables will be regressed out from the data before
#' finding KNN. The default is set to NULL, indicates the regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is set to 'TRUE'.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay before applying log transformation to avoid -Inf for measurements that are equal to 0. The default is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to set to 'TRUE'. See the checkSeObj() function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the KNN results in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to 'TRUE'.
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '|' , assay.name)".
#' @param prps.group Symbol. A symbol specifying the name of the PRPS group. If is 'NULL', the function will select a name
#' based on "paste0('prps|mnn|', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom utils txtProgressBar
#' @importFrom stats dist
#' @importFrom RANN nn2
#' @export


findKnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        clustering.method = 'kmeans',
        nb.clusters = 3,
        k.nn = 3,
        hvg = NULL,
        normalization = 'CPM',
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        output.name = NULL,
        prps.group = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findKnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs #####
    if (is.list(assay.name)) {
        stop('The "assay.name" cannot be a list.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be the name of signle assay in the SummarizedExperiment object.')
    } else if (is.null(uv.variable)) {
        stop('The "uv.variable" variable cannot be empty.')
    } else if (length(uv.variable) > 1) {
        stop('The "uv.variable" must contain the name of signle variable in the SummarizedExperiment object.')
    } else if (!data.input %in% c('expr', 'pcs')) {
        stop('The "data.input" must be one of the "expr" or "pcs".')
    } else if (data.input == 'pcs' & is.null(nb.pcs)) {
        stop('The valuse of "nb.pcs" must be sepcified when the data.input = pcs.')
    } else if (k.nn == 0) {
        stop('The k.nn cannot be 0.')
    } else if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" variable cannot be found in the SummarizedExperiment object.')
    }
    if (!is.null(hvg)) {
        if (sum(hvg %in% row.names(se.obj)) != length(hvg))
            stop('All the hvg genes are not found in the SummarizedExperiment object.')
    }

    # Assessing SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Defining and saving the unwanted variable ####
    all.samples.index <- c(1:ncol(se.obj))
    colnames.seobj <- colnames(se.obj)
    colnames(se.obj) <- paste0(
        colnames.seobj,
        '_',
        all.samples.index
        )

    # Assessing and grouping the unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[uv.variable]]
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
        message = '-- Checking the sample size of each subgroup of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    sub.group.sample.size <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = k.nn + 1
        )
    if (length(sub.group.sample.size) == 0){
        stop(paste0(
            'All subgroups of the unwanted variable have less than ',
            k.nn + 1,
            ' samples. KNN cannot be found.'))
    } else if (length(sub.group.sample.size) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the unwanted variable have less than ',
                k.nn + 1,
                ' samples. Then KNN for those sub-groups cannot be created.'),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the sub-groups of the unwanted variable have at least ',
                k.nn + 1,
                ' samples.'),
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
    if (!is.null(normalization) & is.null(regress.out.variables)) {
        printColoredMessage(
            message = paste0(
                '- applying the ',
                normalization,
                ' normalization on the data.'),
            color = 'blue',
            verbose = verbose
        )
        all.norm.data <- applyOtherNormalizations(
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
                ' normalization and then regressing out ',
                paste0(regress.out.variables, collapse = '&'),
                ' variable(s) from the data.'),
            color = 'blue',
            verbose = verbose
        )
        ### normalization ####
        all.norm.data <- applyOtherNormalizations(
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
        norm.data <- t(all.norm.data)
        lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
        norm.data <- lm(as.formula(paste(
            'norm.data',
            paste0(lm.formua, collapse = '+') ,
            sep = '~'
        )))
        all.norm.data <- t(all.norm.data$residuals)
        colnames(all.norm.data) <- colnames(all.norm.data)
        row.names(all.norm.data) <- row.names(all.norm.data)

    }
    ## regression ####
    if (is.null(normalization) & !is.null(regress.out.variables)){
        if(isTRUE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- applying log2 and then regressing out the',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variables from the data.'),
                color = 'blue',
                verbose = verbose
            )
            if(!is.null(pseudo.count)){
                all.norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                all.norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- regressing out the ',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variable(s) from the data.'),
                color = 'blue',
                verbose = verbose
            )
            all.norm.data <- assay(x = se.obj, i = assay)
        }
        ### regression ####
        sample.info <- as.data.frame(colData(se.obj))
        all.norm.data <- t(all.norm.data)
        lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
        all.norm.data <- lm(as.formula(paste(
            'all.norm.data',
            paste0(lm.formua, collapse = '+') ,
            sep = '~'
        )))
        all.norm.data <- t(all.norm.data$residuals)
        colnames(all.norm.data) <- colnames(all.norm.data)
        row.names(all.norm.data) <- row.names(all.norm.data)
    }
    ## log transformation ####
    if (is.null(normalization) & is.null(regress.out.variables)) {
        if(isTRUE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- applying the log2 on the daya.'),
                color = 'blue',
                verbose = verbose
            )
            if(!is.null(pseudo.count)){
                all.norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                all.norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- no library size normalization and transformation is applied on the data.'),
                color = 'blue',
                verbose = verbose
            )
            all.norm.data <- assay(x = se.obj, i = assay)
        }
    }

    # Select input data for knn analysis ####
    printColoredMessage(
        message = '-- Selecting the data input for knn analysis',
        color = 'magenta',
        verbose = verbose
        )
    all.data.input <- lapply(
        sub.group.sample.size,
        function(x){
            norm.data <- all.norm.data[ , se.obj[[uv.variable]] == x]
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
                        '.'
                        ),
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

    # Finding k nearest neighbors ####
    printColoredMessage(
        message = '-- Finding k nearest neighbors within each sub-group of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    printColoredMessage(
        message = paste0(
            '- For individual samples within the selected sub-group(s) of the variable "',
            uv.variable,
            '", k = ',
            k.nn,
            ' nearest neighbours will be found.'),
        color = 'orange',
        verbose = verbose
        )
    pb <- utils::txtProgressBar(
        min = 0,
        max = length(sub.group.sample.size),
        style = 3
        )
    all.knn <- lapply(
        1:length(sub.group.sample.size),
        function(x) {
            norm.data <- all.data.input[[sub.group.sample.size[x]]]
            ## find knn ####
            printColoredMessage(
                message = paste0(
                    '* Finding k nearest neighbors within the "',
                    sub.group.sample.size[x],
                    '" sub-group.'),
                color = 'blue',
                verbose = verbose
                )
            knn.samples <- RANN::nn2(
                data = norm.data,
                query = norm.data,
                k = c(k.nn + 1),
                treetype = 'bd'
                )
            knn.index <- as.data.frame(knn.samples$nn.idx)
            colnames(knn.index) <- paste0('dataset.index', c(1:c(k.nn + 1)))
            selected.samples <- se.obj[[uv.variable]] == sub.group.sample.size[x]
            ovral.cell.no <- sapply(
                1:ncol(knn.index),
                function(x) all.samples.index[selected.samples][knn.index[ , x]])
            colnames(ovral.cell.no) <- paste0('overal.index', 1:c(k.nn + 1))
            knn.index <- as.data.frame(cbind(ovral.cell.no , knn.index))
            # distance between all knn ####
            ## find knn ####
            printColoredMessage(
                message = paste0(
                    '* Calculating all pairwise distances between the k nearest neighbors within the "',
                    sub.group.sample.size[x],
                    '" subgroup.'),
                color = 'blue',
                verbose = verbose
            )
            knn.dis <- round(as.data.frame(knn.samples$nn.dists), digits = 3)
            knn.dis <- knn.dis[,-1, drop = FALSE]
            colnames(knn.dis) <- paste0('distance1_', 2:c(k.nn + 1))
            knn.index.dist <- cbind(knn.index, knn.dis)
            if (isTRUE(k.nn > 1)) {
                all.comb <- combn(x = paste0('dataset.index', 2:c(k.nn + 1)), m = 2)
                all.comb.names <- combn(x = 2:c(k.nn + 1), m = 2)
                for (z in 1:ncol(all.comb)) {
                    pair.dist <- unlist(lapply(
                        1:nrow(knn.index.dist),
                        function(y) {
                            col1 <- knn.index.dist[, all.comb[, z][1]][y]
                            col2 <- knn.index.dist[, all.comb[, z][2]][y]
                            stats::dist(
                                x = norm.data[c(col1, col2) ,],
                                method = "euclidean",
                                diag = FALSE,
                                upper = FALSE
                            )
                        }))
                    name <- paste0('dist',
                               all.comb.names[, z][1],
                               '_',
                               all.comb.names[, z][2])
                    knn.index.dist[[name]] <- pair.dist
                }
            }
            cols.index <- grep('distance', colnames(knn.index.dist))
            knn.index.dist$aver.dist <- rowMeans(knn.index.dist[, cols.index, drop = FALSE])
            set.seed(4589)
            knn.index.dist$rank.aver.dist <- rank(
                x = knn.index.dist$aver.dist,
                ties.method = 'random'
                )
            knn.index.dist$group <- sub.group.sample.size[x]
            setTxtProgressBar(pb, x)
            message(' ')
            return(knn.index.dist)
        })
    all.knn <- do.call(rbind, all.knn)
    se.obj[[uv.variable]] <- initial.variable

    # Saving the results ####
    ## select prps name ####
    if(is.null(prps.group)){
        prps.group <- paste0('prps|knnMnn|', uv.variable)
    }
    ## select output name ####
    if (is.null(output.name))
        output.name <- paste0(uv.variable, '|' , assay.name)

    ## save the results in the SummarizedExperiment object ####
    message(' ')
    printColoredMessage(
        message = '-- Saving the results:',
        color = 'magenta',
        verbose = verbose
        )
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
        if (!'KnnMnn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]] )) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']] <- list()
        }
        if (!'knn' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['knn']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['knn']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['knn']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['KnnMnn']][['knn']][[output.name]] <- all.knn

        printColoredMessage(
            message = '- All the knn results are saved in the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The findKnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the results as matrix  ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '- All the knn results are outputed as matrix.',
                            color = 'blue',
                            verbose = verbose)
        printColoredMessage(message = '------------The findKnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.knn)
    }
}



