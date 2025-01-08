#' Create PRPS for a categorical source of unwanted variation.

#' @author Ramyar Molania

#' @description
#' Distinct group of PRPS are created for each source of unwanted variation defined in the 'main.uv.variable' argument
#' within an homogeneous group of samples. The grouping of samples are created based on each biological subtype defined
#' using the 'bio.variables' using the 'bio.clustering.method' selected and it might be also be combined with
#' 'other.uv.variables' if requested. By default 'other.uv.variables' is set up to NULL, meaning that the homogeneous
#' grouping of samples created is based solely on the biological subtype(s) defined in 'bio.variables'. If
#' 'other.uv.variables' is not set to NULL, the creation of homogeneous groups of samples is based on the biological subtype
#' defined in 'bio.variables' combined with 'other.uv.variables' using the 'other.uv.clustering.method' selected. A
#' pseudo-sample will be created for each homogeneous group of samples by averaging all the samples within that group of
#' samples.For example to correct for batch effect arising related to a batch variable defined in the 'main.uv.variable'
#' argument, several group of pseudo-samples will be created, one for each homogeneous group of samples and for each batch.
#' All those pseudo-samples belonging to the same group across batches will be defined as pseudo-replicates which constitutes
#' a PRPS set.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol specifying an assay (data) name within the SummarizedExperiment object. The selected
#' assay should be the one that will be used for RUV-III normalization.
#' @param bio.variables Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
#' subtypes, tumor purity, ... within the sample annotation (colData) of the SummarizedExperiment object. This can comprise
#' a vector containing either categorical, continuous, or a combination of both variables.
#' @param main.uv.variable Symbols. A symbol representing the label of a categorical unwanted variable, such as batch effects,
#' , ... within the sample annotation (colData) of the SummarizedExperiment object. Sets of PRPS data will be generated
#' to remove the impact of the this variable form the data.
#' @param other.uv.variables  Symbol. A symbol or symbols representing the label of unwanted variable(s), such as
#' library size, ... within the sample annotation (colData) of the SummarizedExperiment object. This can comprise a vector
#' containing either categorical, continuous, or a combination of both variables. These variable(s) will be considered when
#' generating PRPS sets fo the "main.uv.variable". This can avoid potential contamination for the "other.uv.variables" in
#' the PRPS data of the "main.uv.variable". The default is set to 'NULL'.
#' @param min.sample.for.ps Numeric. A numeric value that indicates the minimum number of biologically homogeneous samples
#' to be averaged to create one pseudo-sample (PS). The default it is set to 3.
#' @param bio.clustering.method Symbol. A symbol indicating which clustering method should be used to group each continuous
#' biological variable, if there are any. The option include 'kmeans', cut' and 'quantile'. The default is to 'kmeans'.
#' We refer to the createHomogeneousBioGroups() function for more details.
#' @param nb.bio.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous biological
#' variable specified in the "bio.variables". The default it set to 3.
#' @param other.uv.clustering.method Symbol. A symbol indicating which clustering method shoudl be used to group each
#' continuous unwanted variable, if there are specified in the "other.uv.variables". The option include 'kmeans', cut'
#' and 'quantile'. The default is to 'kmeans'. We refer to the createHomogeneousUVGroups() function for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous unwanted
#' variable. The default it set to 3.
#' @param check.prps.connectedness Logical. Indicates whether to assess the connectedness between the PRPS sets or not.
#' The default is set to "TRUE". See the details for more information.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before creating PS. The default
#' is set to 'TRUE'.
#' @param pseudo.count Numeric. A numeric value as a pseudo count to be added to all measurements before log transformation.
#' The default is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. If 'TRUE', the checkSeObj()
#' function will be applied inside the function. The default is set to 'TRUE'.
#' @param remove.na Symbol. A symbol that is indicating whether to remove NA or missing values from either the 'assays',
#'  the sample.annotation', both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be
#' excluded. If sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and
#' 'uv.variables' will be excluded. The default is set to both'.
#' @param assess.variables Logical. Indicates whether to assess the association between the biological or unwanted
#' variable(s) separately. The default is set to 'FALSE'. We refer to the 'assessVariableAssociation' for more details.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each
#' pair of categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each
#' pair of 'bio.variables'. The correlation is computed by the function ContCoef from the DescTools package. If the correlation
#' of a pair of variable is higher than the cut-off, then only the variable that has the highest number of factor will be
#' kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient
#' between each pair of continuous variables. The first one is between each pair of 'uv.variables' and the second one is
#' between each pair of 'bio.variables'. If the correlation of a pair of variable is higher than the cut-off, then only
#' the variable that has the highest variance will be kept and the other one will be excluded from the remaining analysis.
#' By default they are both set to 0.7.
#' @param point.size Numeric. A numeric value of the size of points in the PRPS map. The default is set to 4.
#' @param save.se.obj Logical. Indicates whether to save the results in the metadata of the SummarizedExperiment object
#' or to output the result as list. The default by is set to 'TRUE'.
#' @param plot.prps.map Logical. Indicates whether to generate the PRPS map plot for individual sources of unwanted variation.
#' The default is 'TRUE'.
#' @param output.name Symbol. A symbol to specify a name for the PRPS data. The default is set to 'NULL'. Refer to to
#' the detail for more information.
#' @param prps.group Symbol. A symbol that is specify an name for PRPS data sets created for a specific group. The default
#'  is set to 'NULL'. Refer to to the detail for more information.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object that contains all the PRPS data and PPRS map plot in the metadata or a list
#' that contains all the results.

#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix rowMeans
#' @importFrom dplyr count n
#' @importFrom tidyr %>%
#' @import ggplot2
#' @export

createPrPsForCategoricalUV <- function(
        se.obj,
        assay.name,
        bio.variables,
        main.uv.variable,
        other.uv.variables = NULL,
        min.sample.for.ps = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        check.prps.connectedness = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'none',
        assess.variables = FALSE,
        cat.cor.coef = c(0.95, 0.95),
        cont.cor.coef = c(0.95, 0.95),
        plot.prps.map = TRUE,
        point.size = 4,
        save.se.obj = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsForCategoricalUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####
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

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, bio.variables, other.uv.variables),
            remove.na = remove.na
            )
    }

    # Data transformation ####
    printColoredMessage(
        message = '-- Applying log transformation on the data:',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count, ' (pseudo.count) on the ',
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
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(x = se.obj, i = assay.name)
    }

    # Creating PRPS data sets ####
    printColoredMessage(
        message = '-- Creating PRPS sets:',
        color = 'magenta',
        verbose = verbose
        )
    ## creating PRPS with considering other uv variables ####
    if(!is.null(other.uv.variables)){
        printColoredMessage(
            message = '- Creating PRPS sets with considering the variable(s) provided by "other.uv.variables" :',
            color = 'blue',
            verbose = verbose
            )
        ### create homogeneous sample groups with respect to biology ####
        printColoredMessage(
            message = '-- Creating homogeneous sample groups with respect to the variable(s) provided by "bio.variables":',
            color = 'magenta',
            verbose = verbose
            )
        homo.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            cont.cor.coef = cont.cor.coef,
            cat.cor.coef = cat.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        if(sum(table(homo.bio.groups) >=  2*min.sample.for.ps) == 0)
            stop(paste0(
            'All the homogeneous biological groups of samples have less than ',
            2*min.sample.for.ps,
            ' samples. Then, PRPS cannot be created.')
            )

        ## create homogeneous sample groups with respect to unwanted variables ####
        printColoredMessage(
            message = '-- Creating homogeneous sample groups with respect to the variable(s) provided "other.uv.variables":',
            color = 'blue',
            verbose = verbose
            )
        ## create all possible sample groups with respect to other.uv.variables ####
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
        ## putting all biological and unwanted groups together ####
        all.groups <- data.frame(
            homo.uv.groups = homo.uv.groups,
            homo.bio.groups = homo.bio.groups
            )
        bio.batch <- uv.group <- NULL
        all.groups$bio.batch <- paste(
            all.groups$homo.uv.groups,
            all.groups$homo.bio.groups,
            sep = '||')
        all.groups$uv.group <- se.obj[[main.uv.variable]]
        printColoredMessage(
            message = paste0('- ',
                length(unique(all.groups$bio.batch)),
                ' of ',
                c(length(unique(all.groups$homo.bio.groups)) *length(unique(all.groups$homo.uv.groups))),
                ' possible groups (',
                length(unique(all.groups$homo.bio.groups)),
                ' biological groups * ',
                length(unique(all.groups$homo.uv.groups)),
                ' unwanted variation groups) ',
                'groups are created.'),
            color = 'blue',
            verbose = verbose
            )
        ## checking the samples distribution ####
        samples.dis <- table(all.groups$bio.batch, all.groups$uv.group)
        selected.groups <- sum(rowSums(samples.dis >= min.sample.for.ps) > 1)
        if(selected.groups == 0){
            stop(paste0(
                ' None of the ',
                length(unique(all.groups$bio.batch)),
                ' homogeneous samples have at least ',
                min.sample.for.ps,
                ' samples in at least two batches of ',
                main.uv.variable, '. PRPS cannot be created.'))
        }
        printColoredMessage(
            message = paste0(
                selected.groups,
                ' sample groups of ',
                length(unique(all.groups$bio.batch)),
                ' homogeneous samples have at least ',
                min.sample.for.ps,
                ' samples in at least two batches of the ',
                main.uv.variable,
                'variable.'),
            color = 'blue',
            verbose = verbose
            )
        ## check connection between PRPS sets ####
        if(isTRUE(check.prps.connectedness)){
            printColoredMessage(
                message = paste0(
                    '-- Check the connection between possible PRPS sets across batches of the ',
                    main.uv.variable,
                    ' variable.'),
                color = 'magenta',
                verbose = verbose )
            samples.dis <- checkPRPSconnectedness(
                data.input = samples.dis,
                min.samples = min.sample.for.ps,
                batch.name = main.uv.variable
                )
        }
        ## checking samples abundance ####
        printColoredMessage(
            message = '-- Checking samples abundance within each groups before creating PRPS sets:',
            color = 'magenta',
            verbose = verbose
            )
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.ps) > 1 , ]
        mean.samples <- round(mean(samples.dis[samples.dis >= min.sample.for.ps]), digits = 0)
        printColoredMessage(
            message = paste0(
                'The average sample size of groups to select ',
                min.sample.for.ps,
                ' samples for creating pseudo-sample is ',
                mean.samples,
                '.'),
            color = 'red',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                'Please note, the high average sample size,',
                ' may result in a contamination of PRPS sets with some other variation.'),
            color = 'red',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                'To avoid this, you may consider increase the number ',
                'of clustering of continuous biological or unwanted variables or both.'),
            color = 'red',
            verbose = verbose
            )
        ## creating PRPS sets ####
        prps.sets <- lapply(
            1:nrow(samples.dis),
            function(x) {
                selected.batches <- colnames(samples.dis)[samples.dis[x , ] >= min.sample.for.ps]
                ps.matrix <- sapply(
                    selected.batches,
                    function(y) {
                        index.samples <- all.groups$bio.batch == row.names(samples.dis)[x] & all.groups$uv.group == y
                        rowMeans(expr.data[, index.samples])
                    })
                colnames(ps.matrix) <- rep(
                    paste(row.names(samples.dis)[x], main.uv.variable, sep = '||'),
                    ncol(ps.matrix))
                ps.matrix
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0(
                length(unique(colnames(prps.sets))),
                ' PRPS sets with the total number of ',
                ncol(prps.sets),
                ' pseudo-samples are created.'),
            color = 'blue',
            verbose = verbose)
    }
    ## create PRPS without considering other uv variables ####
    if(is.null(other.uv.variables)) {
        printColoredMessage(
            message = '- Creating PRPS sets without considering any other uv variables:',
            color = 'blue',
            verbose = verbose
            )
        homo.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            cont.cor.coef = cont.cor.coef,
            cat.cor.coef = cat.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        if(sum(table(homo.bio.groups) == 1) == length(unique(homo.bio.groups)))
            stop('All the homogeneous biological group of samples have only 1 sample. PRPS cannot be created.')
        all.groups <- data.frame(
            uv.group = se.obj[[main.uv.variable]],
            bio.groups = homo.bio.groups
            )
        samples.dis <- table(all.groups$bio.groups, all.groups$uv.group)

        ### check connection between PRPS sets ####
        if(isTRUE(check.prps.connectedness)){
            printColoredMessage(
                message = paste0(
                    '-- Checking the connection between possible PRPS sets across the batches of ',
                    main.uv.variable,
                    '.'),
                color = 'magenta',
                verbose = verbose
                )
            samples.dis <- checkPRPSconnectedness(
                data.input = samples.dis,
                min.samples = min.sample.for.ps,
                batch.name = main.uv.variable
                )
        }
        ## check samples abundance ####
        printColoredMessage(
            message = '-- Checking samples abundance of groups before create PRPS:',
            color = 'magenta',
            verbose = verbose
            )
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.ps) > 1 , ]
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.ps) > 1 , ]
        mean.samples <- round(mean(samples.dis[samples.dis >= min.sample.for.ps]), digits = 0)
        printColoredMessage(
            message = paste0(
                'The average sample size of groups to select ',
                min.sample.for.ps, ' samples for PRPS is ',
                mean.samples, '.'),
            color = 'red',
            verbose = verbose)
        printColoredMessage(
            message = 'Please note, the high average sample size, may result in a contamination of PRPS with some other variation.',
            color = 'red',
            verbose = verbose
            )
        printColoredMessage(
            message = 'To avoid this, you may consider increase the number of clusters of continuous biological and unwanted variables.',
            color = 'red',
            verbose = verbose
            )
        ## create PRPS sets ####
        prps.sets <- lapply(
            1:nrow(samples.dis),
            function(x) {
                selected.batches <- colnames(samples.dis)[samples.dis[x , ] >= min.sample.for.ps]
                ps.matrix <- sapply(
                    selected.batches,
                    function(y) {
                        index.samples <- all.groups$bio.groups == row.names(samples.dis)[x] &
                            all.groups$uv.group == y
                        rowMeans(expr.data[, index.samples])
                    })
                colnames(ps.matrix) <- rep(
                    paste(row.names(samples.dis)[x], main.uv.variable, sep = '||'),
                    ncol(ps.matrix))
                ps.matrix
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0(
            length(unique(colnames(prps.sets))),
            ' PRPS sets with the total number of ',
            ncol(prps.sets),
            ' pseudo-samples are created.'),
            color = 'blue',
            verbose = verbose)
    }
    ## plot PRPS map #####
    if (isTRUE(plot.prps.map)) {
        ## PRPS map plot
        printColoredMessage(
            message = '-- Plotting the PRPS map:',
            color = 'magenta',
            verbose = verbose
            )
        if(!is.null(other.uv.variables)){
            groups <- all.groups$bio.batch
        } else groups <- homo.bio.groups
        catvar <- use <- n <- NULL
        sample.info <- as.data.frame(
            x = colData(se.obj),
            stringsAsFactors = FALSE
            )
        colnames(sample.info) <- names(colData(se.obj))
        sample.info$catvar <- as.factor(paste0(sample.info[, main.uv.variable]))
        sample.info$groups <- as.factor(groups)
        sample.info <- sample.info %>% dplyr::count(catvar, groups)
        sample.info$new.name <- paste0(sample.info$catvar, sample.info$groups)
        sample.info.filtered <- sample.info %>%
            group_by(groups) %>%
            filter(n >= min.sample.for.ps) %>%
            filter(n() >= 2)
        sample.info.filtered$use <- 'Selected'
        sample.info <- left_join(
            x = sample.info,
            y = sample.info.filtered[ , c('new.name', 'use')],
            by = 'new.name'
            )
        sample.info$use[is.na(sample.info$use)] <- 'Unselected'
        prps.map.plot <- ggplot(sample.info, aes(x = catvar, y = groups)) +
            geom_point(aes(color = use), size = point.size) +
            geom_text(aes(label = n , hjust = 0.5, vjust = 0.5)) +
            xlab(main.uv.variable) +
            ylab('Homogeneous groups') +
            theme_bw() +
            theme(
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                axis.line = element_line(colour = 'black', linewidth = .85),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                axis.text.y = element_text(size = 12),
                legend.position = 'right') +
            guides(color = guide_legend(title = "PS"))
        if(isTRUE(verbose)) print(prps.map.plot)
    }
    # Save the results ####
    ## select output name ####
    if(!is.null(other.uv.variables)) {
        output.name <- paste0(
            main.uv.variable,
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(other.uv.variables, collapse = '&'),
            '|',
            assay.name)
    } else{
        output.name <- paste0(
            main.uv.variable,
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            assay.name)
    }
    if (is.null(prps.group)) {
        prps.group <- paste0('prps|supervised|', main.uv.variable)
    }
    ## save the output in the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        ## check
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## check
        if (!'supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['supervised']] <- list()
        }
        ## check
        if (!prps.group %in% names(se.obj@metadata[['PRPS']][['supervised']])) {
            se.obj@metadata[['PRPS']][['supervised']][[prps.group]] <- list()
        }
        ## check
        if (!'prps.data' %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.data']] <- list()
        }
        ## check
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.data']])) {
            se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.data']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.data']][[output.name]] <- prps.sets

        ## plot
        if(isTRUE(plot.prps.map)){
            if (!'prps.map.plot' %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]])) {
                se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']])) {
                se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']][[output.name]] <- list()
            }
            se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']][[output.name]] <- prps.map.plot
        }
        printColoredMessage(
            message = 'The PRPS data and PRPS map plot are saved to th metdata in the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The createPrPsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose
                            )
        return(se.obj)
    }
    ## save the output in as data frame ####
    if(isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose)
        if(isTRUE(plot.prps.map)){
            return(list(
                prps.sets = prps.sets,
                prps.map.plot = prps.map.plot)
                )
        } else if (isFALSE(plot.prps.map)){
            return(list(prps.sets = prps.sets))
        }

    }
}


