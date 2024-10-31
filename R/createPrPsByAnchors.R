#' Create PRPS sets across batches using integration anchors.

#' @author Ramyar Molania

#' @description
#' This functions employs the FindIntegrationAnchors function from the Seurat R package to create PRPS sets for the
#' RUV-III normalization of RNA-seq data. This function can be used in situations in which biological variation are unknown.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol indicating the assay name within the SummarizedExperiment object for the creation
#' of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variable Symbol. A symbol specifying the name of columns in the sample annotation of the SummarizedExperiment
#' object. This variable can to be categorical or continuous variable. If a continuous variable is provide, this will be
#' divided into groups using the clusteintf method.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param hvg Vector. A vector containing the names of highly variable genes. These genes will be utilized to identify
#' anchor samples across different batches. The default value is set to 'NULL'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay before applying
#' log transformation to avoid -Inf for raw counts that are equal to 0. The default is 1.
#' @param anchor.features Numeric. A numeric value indicating the provided number of features to be used in anchor finding.
#' The default is 2000. We refer to the FindIntegrationAnchors R function for more details.
#' @param scale Logical. Whether or not to scale the features provided. Only set to FALSE if you have previously scaled
#' the features you want to use for each object in the object.list.
#' @param min.ps.samples Numeric. The minimum number of samples to be averaged to create a pseudo-sample. The default
#' is 3. The minimum value is 2.
#' @param max.ps.samples Numeric value indicating the maximum number of samples to be averaged for creating a pseudo-sample.
#' The default is 'inf'. Please note that averaging a high number of samples may lead to inaccurate PRPS estimates.
#' @param max.prps.sets Numeric. The maximum number of PRPS sets across batches. The default in 10.
#' @param normalization.method Symbol. Indicate which normalization methods should be used before finding the anchors.
#' The options are "LogNormalize" or "SCT". The default is "LogNormalize".
#' @param sct.clip.range Numeric. Numeric of length two specifying the min and max values the Pearson residual will be
#' clipped to. The default is 'NULL'. We refer to the FindIntegrationAnchors R function for more details.
#' @param reduction Symbol. Indicates which dimensional reduction to perform when finding anchors. The options are "cca":
#' canonical correlation analysis, "rpca": reciprocal PCA and "rlsi": Reciprocal LSI. The default is "cca".
#' @param l2.norm Logical. Indicates whether to perform L2 normalization on the CCA samples embeddings after dimensional
#' reduction or not. The default is 'TRUE'.
#' @param dims Numeric. Indicates which dimensions to use from the CCA to specify the neighbor search space. Th default
#' is 10.
#' @param k.anchor Numeric. How many neighbors (k) to use when picking anchors. Th default is set to 2.
#' @param k.filter Numeric. How many neighbors (k) to use when filtering anchors. Th default is 20.
#' @param k.score Numeric. How many neighbors (k) to use when scoring anchors. Th default is 30.
#' @param max.features Numeric. The maximum number of features to use when specifying the neighborhood search space in
#' the anchor filtering.Th default is 30.
#' @param nn.method Symbol.Method for nearest neighbor finding. Options include: "rann", "annoy". The defauly is "annoy".
#' @param n.trees Numeric. More trees gives higher precision when using annoy approximate nearest neighbor search
#' @param eps Numeric. Error bound on the neighbor finding algorithm (from RANN/Annoy)
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment
#' object or to output the result as list. By default it is set to TRUE.
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '|', 'anchor', '|', assay.name))".
#' @param prps.group Symbol. A symbol specifying the name of the PRPS group. If is 'NULL', the function will select a name
#' based on "paste0('prps|anchor|', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom Seurat VariableFeatures FindIntegrationAnchors
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom dplyr group_by top_n desc
#' @importFrom Matrix rowMeans
#' @importFrom stats setNames
#' @export

createPrPsByAnchors <- function(
        se.obj,
        assay.name,
        uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        hvg = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        anchor.features = 2000,
        scale = TRUE,
        min.ps.samples = 3,
        max.ps.samples = 5 ,
        max.prps.sets = 10,
        normalization.method = "LogNormalize",
        sct.clip.range = NULL,
        reduction = "cca",
        l2.norm = TRUE,
        dims = 1:15,
        k.anchor = 2,
        k.filter = 10,
        k.score = 30,
        max.features = 200,
        nn.method = "annoy",
        n.trees = 50,
        eps = 0,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsByAnchors function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check input #####
    if (is.list(assay.name)) {
        stop('The "assay.name" must be the name of an assay in the SummarizedExperiment object.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single name of an assay in the SummarizedExperiment object..')
    } else if (is.null(uv.variable)) {
        stop('The "uv.variable" cannot be empty. Note, unknown sources of UV can be found by the "identifyUnknownUV" function.')
    } else if (length(uv.variable) > 1) {
        stop('The "uv.variable" must be a single variable name in the SummarizedExperiment object.')
    } else if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" cannot be found in the SummarizedExperiment object')
    } else if (k.anchor == 0) {
        stop('The k.anchor cannot be 0.')
    } else if (!normalization.method %in% c("LogNormalize", "SCT")) {
        stop('The "normalization.method" must be one of the "LogNormalize" or "SCT".')
    } else if (!reduction %in% c("cca", "rpca", 'rlsi')) {
        stop('The "reduction" should be one of the "cca", "rpca" or "rlsi"')
    } else if (sum(hvg %in% row.names(se.obj)) != length(hvg)) {
        stop('All the "hvg" genes are not found in the SummarizedExperiment object.')
    }

    # Assess the SummarizedExperiment ####
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

    # Find anchors ####
    printColoredMessage(
        message = '-- Find the anchors between all the pairs of batches:',
        color = 'magenta',
        verbose = verbose
        )
    ## split the data into groups ####
    printColoredMessage(
        message = paste0(
            '* split the SummarizedExperiment object into ', length(levels(se.obj[[uv.variable]])),
            ' groups and select highly variable genes.'),
        color = 'blue',
        verbose = verbose
        )
    groups <- levels(se.obj[[uv.variable]])
    all.seurat.objects <- lapply(
        groups,
        function(x) {
            samples.index <- se.obj[[uv.variable]] == x
            seurat.obj <- SeuratObject::CreateSeuratObject(
                counts = assay(x = se.obj[, samples.index], i = assay.name),
                project = x
                )
            seurat.obj <- Seurat::NormalizeData(object = seurat.obj)
            if (!is.null(hvg)){
                Seurat::VariableFeatures(seurat.obj) <- hvg
            } else {
                seurat.obj <- Seurat::FindVariableFeatures(object = seurat.obj)
            }
            return(seurat.obj)
        })
    names(all.seurat.objects) <- groups
    all.samples.index <- c(1:ncol(se.obj))

    ## find anchors  ####
    printColoredMessage(
        message = '* find the anchors.',
        color = 'blue',
        verbose = verbose
        )
    all.anchors <- Seurat::FindIntegrationAnchors(
        object.list = all.seurat.objects,
        anchor.features = anchor.features,
        reference = NULL,
        scale = scale,
        normalization.method = normalization.method,
        sct.clip.range = sct.clip.range,
        reduction = reduction,
        l2.norm = l2.norm,
        dims = dims,
        k.anchor = k.anchor,
        k.filter = k.filter,
        k.score = k.score,
        max.features = max.features,
        nn.method = nn.method,
        n.trees = n.trees,
        eps = eps,
        verbose = verbose
        )
    score <- NULL
    all.anchors <- all.anchors@anchors
    all.anchors$score <- round(x = all.anchors$score, digits = 3)
    colnames(all.anchors)[1:2] <- c('sample1', 'sample2')
    all.anchors$dataset1.name <- groups[all.anchors$dataset1]
    all.anchors$dataset2.name <- groups[all.anchors$dataset2]
    printColoredMessage(
        message = paste0(
            nrow(all.anchors)/2, ' sample pairs (anchors) are found across all the subgroups of the "',
            uv.variable, '" variable.'),
        color = 'blue',
        verbose = verbose
        )

    ## sanity check  ####
    printColoredMessage(
        message = '* sanity check on the anchors.',
        color = 'blue',
        verbose = verbose
        )
    sanity.check <- lapply(
        1:length(groups),
        function(x) {
            max.index <- max(all.anchors$sample1[all.anchors$dataset1 == x])
            sample.size <- sum(se.obj[[uv.variable]] == groups[x])
            if (max.index > sample.size) {
                stop('There are someting wrong with the order of anchros.')
            }
            max.index <- max(all.anchors$sample2[all.anchors$dataset2 == x])
            sample.size <- sum(se.obj[[uv.variable]] == groups[x])
            if (max.index > sample.size) {
                stop('There are someting wrong with the order of anchros.')
            }
        })

    ## add overall sample index  ####
    printColoredMessage(
        message = '- Add overall sample number to the anchors.',
        color = 'blue',
        verbose = verbose
        )
    for (x in 1:length(groups)) {
        index.a <- all.anchors$sample1[all.anchors$dataset1 == x]
        all.anchors$sample.index1[all.anchors$dataset1 == x] <-
            all.samples.index[se.obj[[uv.variable]] == groups[x]][index.a]
        index.b <- all.anchors$sample2[all.anchors$dataset2 == x]
        all.anchors$sample.index2[all.anchors$dataset2 == x] <-
            all.samples.index[se.obj[[uv.variable]] == groups[x]][index.b]
    }
    rm(all.seurat.objects)
    gc()

    # Find possible sets of PRPS ####
    printColoredMessage(
        message = '- Find all possible sets of PRPS using the anchors:',
        color = 'magenta',
        verbose = verbose
        )
    rep.anchors <- nrow(all.anchors)/2
    half.anchors <- all.anchors$sample.index1[1:rep.anchors]
    second.half.anchors <- all.anchors$sample.index2[c(rep.anchors+1):nrow(all.anchors)]

    if(isTRUE(all.equal(half.anchors , second.half.anchors))){
        all.anchors <- all.anchors[1:rep.anchors, ]
    } else stop('something wrong with the anchors.')

    all.prps.sets <- split(
        x = all.anchors,
        f = all.anchors$sample.index1
        )
    groups.sets <- names(all.prps.sets)
    all.prps.sets <- lapply(
        groups.sets,
        function(x) {
            temp.anchors <- all.prps.sets[[x]]
            temp.anchors <- rbind(
                temp.anchors,
                do.call(rbind,
                    lapply(temp.anchors$sample.index2, function(j)
                        all.anchors[all.anchors$sample.index2 == j,])
                    ))
            all.datasets <- sort(unique(c(temp.anchors$dataset1, temp.anchors$dataset2)))
            temp.anchors <- temp.anchors[order(temp.anchors$score, decreasing = TRUE) , ]
            anchor.sets <- lapply(
                all.datasets,
                function(x) {
                    unique(c(
                        temp.anchors$sample.index1[temp.anchors$dataset1 == x],
                        temp.anchors$sample.index2[temp.anchors$dataset2 == x]
                    ))
                })
            average.scores <- mean(temp.anchors$score)
            anchor.datasets <-
                unlist(lapply(all.datasets, function(x) {
                    unique(c(
                        temp.anchors$dataset1.name[temp.anchors$dataset1 == x],
                        temp.anchors$dataset2.name[temp.anchors$dataset2 == x]
                    ))
                }))
            length.sets <- sapply(anchor.sets, length)
            return(
                list(average.scores = average.scores,
                    anchor.sets = setNames(anchor.sets, anchor.datasets),
                    length.sets = length.sets))
        })
    names(all.prps.sets) <- paste0('Anchor', groups.sets)

    ## filter ps set
    keep.anchor.sets <- sapply(
        names(all.prps.sets),
        function(x) sum(all.prps.sets[[x]]$length.sets >= min.ps.samples) >= 2)

    all.prps.sets <- all.prps.sets[keep.anchor.sets]
    printColoredMessage(
        message = paste0('* ', length(all.prps.sets), ' possible PRPS stes are found.'),
        color = 'blue',
        verbose = verbose
        )

    ## check initial coverage ####
    printColoredMessage(
        message = '- Check the distribution of the PRPS sets across the batches.',
        color = 'orange',
        verbose = verbose
        )
    prps.coverage <- matrix(0, nrow = length(all.prps.sets), ncol = length(groups))
    colnames(prps.coverage) <- groups
    prps.coverage <- lapply(
        seq_along(all.prps.sets),
        function(i) {
            index <- match(names(prps.coverage[i, ]), names(all.prps.sets[[i]]$anchor.sets))
            all.prps.sets[[i]]$length.sets[index]
    })
    prps.coverage <- do.call(rbind, prps.coverage)
    colnames(prps.coverage) <- groups
    prps.coverage[is.na(prps.coverage)] <- 0
    if (isTRUE(sum(colSums(prps.coverage) == 0) > 0 )) {
        printColoredMessage(
            message = paste(
                paste0(colnames(prps.coverage)[colSums(prps.coverage) == 0], collapse = ' & '),
                'are not covered by any PRPS set.'),
            color = 'red',
            verbose = verbose)
    }
    if (isTRUE(sum(rowSums(prps.coverage >= min.ps.samples) == length(groups)) > 0)) {
        printColoredMessage(
            message = paste0('* there are ', sum(rowSums(prps.coverage >= min.ps.samples) > 1),
                ' PRPS sets with at least ', min.ps.samples,
                ' samples within at least two batches across all subgroups of ', uv.variable, '.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0('* there are ', sum(rowSums(prps.coverage >= min.ps.samples) == length(groups)),
                ' PRPS sets with at least ', min.ps.samples, ' samples within each batch across all subgroups of the ',
                uv.variable, '" variable.'),
            color = 'blue',
            verbose = verbose
            )
    } else {
        # check connection
        printColoredMessage(
            message = 'There is no any PRPS sets that cover all batches. We assess the connection between differet sets:',
            color = 'blue',
            verbose = verbose
            )
        prps.connection <- list()
        for(y in 1:nrow(prps.coverage)){
            batch.names.a <- names(which(prps.coverage[y, ] >= min.ps.samples))
            con.prps <- lapply(
                c(1:nrow(prps.coverage))[-y],
                function(z) {
                    batch.names.b <- names(which(prps.coverage[z,] >= min.ps.samples))
                    inter.samples <- intersect(batch.names.a, batch.names.b)
                    if (length(inter.samples) > 0) {
                        sort(unique(c(batch.names.a, batch.names.b)), decreasing = FALSE)
                    } else
                        sort(batch.names.a, decreasing = FALSE)
                })
            all <- sort(unique(unlist(Filter(Negate(is.null), con.prps))), decreasing = FALSE)
            prps.connection[[y]] <- all
            if (sum(all %in% groups) == length(groups))
                break
        }
        covered.batches <- unique(unlist(Filter(Negate(is.null), prps.connection)))
        not.covered.batches <- groups[!groups %in% covered.batches]
        if (length(covered.batches) == length(groups)) {
            printColoredMessage(
                message = '- The connections between PRPS sets can cover all the batches.',
                color = 'white',
                verbose = verbose
                )
        } else {
            printColoredMessage(
                message = '- All batches are not covered by PRPS.',
                color = 'red',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    'The PRPS sets donot cover ',
                    paste0(not.covered.batches, collapse = ' & '),' batches.'),
                color = 'white',
                verbose = verbose
                )
        }
    }
    prps.coverage <- as.data.frame(prps.coverage)
    group <- NULL
    prps.coverage$group <- sapply(
        1:nrow(prps.coverage),
        function(x) paste0(groups[prps.coverage[x , ] >= min.ps.samples], collapse = '_')
        )
    prps.coverage$score <- unlist(lapply(
        names(all.prps.sets),
        function(x)  all.prps.sets[[x]]$average.scores)
        )
    prps.coverage$prps.sets <- names(all.prps.sets)
    index <- prps.coverage$group == paste0(groups, collapse = '_' )
    if(sum(index) > max.prps.sets){
        prps.coverage <- prps.coverage[index, ]
    }
    selected.prps <- group_by(.data = prps.coverage, group) %>%
        arrange(desc(score)) %>%
        top_n(max.prps.sets, wt = score)
    selected.prps <- selected.prps$prps.sets
    all.prps.sets <- all.prps.sets[selected.prps]

    ## filter PRPS sets based on PS sample size ####
    printColoredMessage(
        message = '- Filter the PRPS sets based on sample size.',
        color = 'blue',
        verbose = verbose
        )
    all.prps.sets.filtered <- lapply(
        names(all.prps.sets),
        function(x) {
            df <- all.prps.sets[[x]]
            for(i in names(df$anchor.sets)){
                if(length(df$anchor.sets[[i]]) > max.ps.samples){
                    df$anchor.sets[[i]] <- df$anchor.sets[[i]][1:max.ps.samples]
                }
            }
            return(df)
        })
    names(all.prps.sets.filtered) <- names(all.prps.sets)

    printColoredMessage(
        message = paste0( 'There are ', length(all.prps.sets), ' PRPS sets with at least ', min.ps.samples,
            ' samples and maximum = ', max.ps.samples, ' within at least two batches across all subgroups of ', uv.variable, '.'),
        color = 'blue',
        verbose = verbose
        )

    # final check coverage ####
    ## check initial coverage ####
    printColoredMessage(
        message = '- Check the distribution of the PRPS sets across the batches.',
        color = 'black',
        verbose = verbose
    )
    prps.coverage <- matrix(0, nrow = length(all.prps.sets), ncol = length(groups))
    colnames(prps.coverage) <- groups
    prps.coverage <- lapply(
        seq_along(all.prps.sets),
        function(i) {
            index <- match(names(prps.coverage[i, ]), names(all.prps.sets[[i]]$anchor.sets))
            all.prps.sets[[i]]$length.sets[index]
        })
    prps.coverage <- do.call(rbind, prps.coverage)
    colnames(prps.coverage) <- groups
    prps.coverage[is.na(prps.coverage)] <- 0
    if (isTRUE(sum(colSums(prps.coverage) == 0) > 0 )) {
        printColoredMessage(
            message = paste(
                paste0(colnames(prps.coverage)[colSums(prps.coverage) == 0], collapse = ' & '),
                'are not covered by any PRPS set.'),
            color = 'red',
            verbose = verbose)
    }
    if (isTRUE(sum(rowSums(prps.coverage >= min.ps.samples) == length(groups)) > 0)) {
        printColoredMessage(
            message = paste0('* there are ', sum(rowSums(prps.coverage >= min.ps.samples) > 1),
                             ' PRPS sets with at least ', min.ps.samples,
                             ' samples within at least two batches across all subgroups of ', uv.variable, '.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0('* there are ', sum(rowSums(prps.coverage >= min.ps.samples) == length(groups)),
                             ' PRPS sets with at least ', min.ps.samples, ' samples within each batch across all subgroups of the ',
                             uv.variable, '" variable.'),
            color = 'blue',
            verbose = verbose
        )
    } else {
        # check connection
        printColoredMessage(
            message = 'There is no any PRPS sets that cover all batches. We assess the connection between differet sets:',
            color = 'blue',
            verbose = verbose
        )
        prps.connection <- list()
        for(y in 1:nrow(prps.coverage)){
            batch.names.a <- names(which(prps.coverage[y, ] >= min.ps.samples))
            con.prps <- lapply(
                c(1:nrow(prps.coverage))[-y],
                function(z) {
                    batch.names.b <- names(which(prps.coverage[z,] >= min.ps.samples))
                    inter.samples <- intersect(batch.names.a, batch.names.b)
                    if (length(inter.samples) > 0) {
                        sort(unique(c(batch.names.a, batch.names.b)), decreasing = FALSE)
                    } else
                        sort(batch.names.a, decreasing = FALSE)
                })
            all <- sort(unique(unlist(Filter(Negate(is.null), con.prps))), decreasing = FALSE)
            prps.connection[[y]] <- all
            if (sum(all %in% groups) == length(groups))
                break
        }
        covered.batches <- unique(unlist(Filter(Negate(is.null), prps.connection)))
        not.covered.batches <- groups[!groups %in% covered.batches]
        if (length(covered.batches) == length(groups)) {
            printColoredMessage(
                message = '- The connections between PRPS sets can cover all the batches.',
                color = 'white',
                verbose = verbose
            )
        } else {
            printColoredMessage(
                message = '- All batches are not covered by PRPS.',
                color = 'red',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'The PRPS sets donot cover ',
                    paste0(not.covered.batches, collapse = ' & '),' batches.'),
                color = 'white',
                verbose = verbose
            )
        }
    }
    ## create PRPS data ####
    # data transformation and normalization ####
    printColoredMessage(
        message = '-- Data transformation and normalization:',
        color = 'magenta',
        verbose = verbose
        )
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0('- apply log2 the ', assay.name,' + ', pseudo.count, ' data.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0('- apply log2 on the ', assay.name,' data.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    prps.data <- lapply(
        names(all.prps.sets),
        function(x) {
            prps.sets <- all.prps.sets[[x]]$anchor.sets
            temp.data <- lapply(
                1:length(prps.sets),
                function(y) {
                    if(length(prps.sets[[y]]) >= min.ps.samples){
                        rowMeans(expr.data[, prps.sets[[y]], drop = FALSE])
                    }
                })
            temp.data <- do.call(cbind, temp.data)
            colnames(temp.data) <- rep(
                x = paste0(uv.variable, '.unsu.', x),
                ncol(temp.data)
                )
            return(temp.data)
        })
    prps.data <- do.call(cbind, prps.data)
    if (sum(table(colnames(prps.data)) == 1)) {
        stop( 'There are something wrong with the prps.data. All the column names of the prps.data are the same.')
    }
    if(sum(is.na(prps.data)) !=0){
        stop( 'There NA in the PRPS data, please check the data input and parameters.')
    }
    se.obj[[uv.variable]] <- initial.variable

    # Saving the output ####
    ## select output name ####
    if(is.null(output.name))
        output.name <- paste0(uv.variable, '|', 'anchor', '|', assay.name)
    if (is.null(prps.group))
        prps.group <- paste0('prps|anchor|', uv.variable)
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
        if (!'prps.data' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- prps.data

        printColoredMessage(message = '------------The createPrPsByAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    if (isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The createPrPsByAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(prps.data = prps.data))
    }
}



