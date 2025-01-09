#' Create all possible homogeneous groups with respect to biological variables.--Finalized

#' @author Ramyar Molania

#' @description
#' This function generates all possible homogeneous sample groups based on the specified biological variables.

#' @details
#' The function generates all possible homogeneous sample groups based on the specified biological variables. If continuous
#' variables are provided, the function splits each into a number of clusters determined by ’nb.clusters’, using the
#' clustering method specified in ’clustering.method’. Ultimately, all combinations of all clusters are created and
#' each such combination is regarded as a homogeneous sample group concerning biological variables.

#' @param se.obj A 'SummarizedExperiment' object.
#' @param bio.variables Character. A character string or a vector of character strings specifying the column names of biological
#' variables in the sample annotation of the 'SummarizedExperiment  object. These 'bio.variables' can be either categorical
#' or continuous variables.
#' @param clustering.method Character. A character string specifying the clustering method to be applied for grouping each
#' continuous biological variable. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans' clustering.
#' @param nb.clusters Numeric. A value indicating the number of groups for continuous sources of biological variation.
#' The default is 3. This implies that each continuous variable will be split into 3 groups using the specified
#' 'clustering.method'.
#' @param assess.se.obj Logical. Whether to assess the 'SummarizedExperiment' object or not. If 'TRUE', the function
#' 'checkSeObj' will be applied. The default is set to ' TRUE' .
#' @param remove.na Character. Indicates whether to remove missing values from the specified variables. The options are
#' 'sample.annotation' or 'none'. The default is set to 'sample.annotation', meaning that missing values in the variables
#' will be removed.
#' @param save.se.obj Logical. Indicates whether to save the results to the metadata of the ' SummarizedExperiment'  object
#' or not. If ' TRUE' , all the possible homogeneous groups will be saved into 'se.obj->metadata->HomogeneousGroups->
#' BiologicalVariables ; otherwise, the results will be returned as a vector. The default is set to ' TRUE' .
#' @param verbose Logical. If ' TRUE' , displays messages for different steps of the function.

#' @return Either a ' SummarizedExperiment' object containing all possible homogeneous groups in 'metadata->HomogeneousGroups
#' ->BiologicalVariables' or a vector of all possible homogeneous sample groups.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable
#' @export

createHomogeneousBioGroups <- function(
        se.obj,
        bio.variables,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The createHomogeneousBioGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check  the inputs ####
    if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty.')
    }
    if (!is.vector(bio.variables)){
        stop('The "bio.variables" must be a vector of variable name(s).')
    }
    if (sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables)) {
        stop('All or some of the "bio.variables" cannot be found in the SummarizedExperiment object.')
    }
    if (!clustering.method %in% c('kmeans', 'cut', 'quantile')){
        stop('The "clustering.method" must be one of: "kmeans", "cut" or "quantile".')
    }
    if (!remove.na %in% c('sample.annotation', 'none.')){
        stop('The "remove.na" mist be either "sample.annotation" or "none".')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical (TRUE or FALSE).')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical (TRUE or FALSE).')
    }

    # Assessing SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        printColoredMessage(
            message = '-- Assessing the SummarizedExperiment object:',
            color = 'magenta',
            verbose = verbose
            )
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = NULL,
            variables = bio.variables,
            remove.na = remove.na,
            verbose = verbose)
        }

    # Finding the class of biological variables ####
    printColoredMessage(
        message = '-- Finding the class and levels of the variable(s):',
        color = 'magenta',
        verbose = verbose
        )
    class.bio.var <- unlist(lapply(
        bio.variables,
        function(x) class(se.obj[[x]]))
        )
    ## categorical variables ####
    categorical.bio.var <- bio.variables[class.bio.var %in% c('factor', 'character')]
    if (length(categorical.bio.var) > 0){
        printColoredMessage(
            message = paste0(
                '- ',
                length(categorical.bio.var),
                ' categorical variable(s) are provided.'),
            color = 'blue',
            verbose = verbose
            )
        categorical.bio.var <- unlist(lapply(
            categorical.bio.var,
            function(x){
                if (length(unique(colData(se.obj))[[x]]) == 1){
                    printColoredMessage(
                        message = paste0(
                            '- the "',
                            x,
                            '" variable contains 1 level. Then, this variable will be excluded.'),
                        color = 'blue',
                        verbose = verbose)
                } else if (length(unique(colData(se.obj))[[x]]) > 1) {
                    printColoredMessage(
                        message = paste0(
                            '* the "',
                            x,
                            '" variable contains ',
                            length(unique(colData(se.obj)[[x]])),
                            ' levels.'),
                        color = 'blue',
                        verbose = verbose)
                    keep.variable <- x
                }
                return(keep.variable)
            }))
    }
    ## continuous variables ####
    continuous.bio.var <- bio.variables[class.bio.var %in% c('numeric', 'integer')]
    if (length(continuous.bio.var) > 0){
        printColoredMessage(
            message = paste0(
                '- ',
                length(continuous.bio.var),
                ' continuous variable(s) are provided.'),
            color = 'blue',
            verbose = verbose
            )
        continuous.bio.var <- unlist(lapply(
            continuous.bio.var,
            function(x){
                if (var(colData(se.obj)[[x]]) == 0){
                    printColoredMessage(
                        message = paste0(
                            '- the variance of the "',
                            x,
                            '" variable is 0 . Then, this will be excluded.'),
                        color = 'blue',
                        verbose = verbose
                        )
                } else if (var(colData(se.obj)[[x]]) > 0) {
                    printColoredMessage(
                        message = paste0(
                            '* the variance of the "',
                            x,
                            '" variable is ~ ',
                            round(x = var(colData(se.obj)[[x]]), digits = 5) , '.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    keep.variable <- x
                }
                return(keep.variable)
            }))
    }

    # Clustering the continuous variable ####
    if (length(continuous.bio.var) > 0){
        ## check clustering value for continuous variable
        if (!is.numeric(nb.clusters)){
            stop('The "nb.clusters" must be a numeric value.')
        }
        if (length(nb.clusters) > 1){
            stop('The "nb.clusters" must be a numeric value.')
        }
        if (nb.clusters == 1 | nb.clusters == 0){
            stop('The value of "nb.clusters" must be bigger than 1.')
        }
        ## cluster continuous variables ####
        printColoredMessage(
            message = '-- Clustering continuous sources of biological variation:',
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '- Each source will be divided into ',
                nb.clusters,
                ' groups using the ',
                clustering.method,
                ' method.'),
            color = 'blue',
            verbose = verbose
        )
        ## kmeans ####
        if (clustering.method == 'kmeans') {
            set.seed(3456)
            continuous.bio.groups <- lapply(
                continuous.bio.var,
                function(x){
                    bio.cont.clusters <- stats::kmeans(
                        x = colData(se.obj)[[x]],
                        centers = nb.clusters,
                        iter.max = 1000)
                    paste0(x, '_group', bio.cont.clusters$cluster)
                })
            names(continuous.bio.groups) <- continuous.bio.var
            continuous.bio.groups <- as.data.frame(do.call(cbind, continuous.bio.groups))
            continuous.bio.groups <- apply(
                continuous.bio.groups,
                1,
                paste , collapse = "..")
        }
        ## cut ####
        if (clustering.method == 'cut') {
            continuous.bio.groups <- lapply(
                continuous.bio.var,
                function(x) {
                    bio.cont.clusters <- as.numeric(
                        cut(x = colData(se.obj)[[x]],
                            breaks = nb.clusters,
                            include.lowest = TRUE))
                    paste0(x, '_group', bio.cont.clusters)
                })
            names(continuous.bio.groups) <- continuous.bio.var
            continuous.bio.groups <- as.data.frame(do.call(cbind, continuous.bio.groups))
            continuous.bio.groups <- apply(
                continuous.bio.groups,
                1,
                paste , collapse = "..")
        }
        ## quantile ####
        if (clustering.method == 'quantile') {
            continuous.bio.groups <- lapply(
                continuous.bio.var,
                function(x) {
                    quantiles <- quantile(
                        x = colData(se.obj)[[x]],
                        probs = seq(0, 1, 1 / nb.clusters)
                        )
                    bio.cont.clusters <- as.numeric(
                        cut(x = colData(se.obj)[[x]],
                            breaks = quantiles,
                            include.lowest = TRUE))
                    paste0(x, '_group', bio.cont.clusters)
                })
            names(continuous.bio.groups) <- continuous.bio.var
            continuous.bio.groups <- as.data.frame(do.call(cbind, continuous.bio.groups))
            continuous.bio.groups <- apply(
                continuous.bio.groups,
                1,
                paste , collapse = "..")
        }
    }

    # Checking categorical variable ####
    if (length(categorical.bio.var) > 0){
        categorical.bio.data <- as.data.frame(colData(se.obj)[, categorical.bio.var, drop = FALSE])
    }
    # Create all possible homogeneous groups ####
    printColoredMessage(
        message = '-- Creating all possible combination of the groups:',
        color = 'magenta',
        verbose =  verbose
        )
    ## use both continuous and categorical ####
    if (length(categorical.bio.var) > 0 & length(continuous.bio.groups) > 0 ){
        printColoredMessage(
            message = '- Creating all possible combination of the groups using both continuous and categorical variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- categorical.bio.data
        all.groups$cont <- continuous.bio.groups
        all.groups <- unname(apply(
            all.groups,
            1,
            paste , collapse = ".."
        ))
        printColoredMessage(
            message = paste0(
                '* in totall, ',
                length(unique(all.groups)),
                ' homogeneous sample groups  are created:'),
            color = 'blue',
            verbose = verbose
        )
        cat("", file = stdout(), sep = "")
        if(isTRUE(verbose))
            print(x = kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## use only continuous ####
    if (length(categorical.bio.var) == 0 & length(continuous.bio.var) > 0){
        printColoredMessage(
            message = '- Creating all possible combination of the groups using continuous variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- continuous.bio.groups
        printColoredMessage(
            message = paste0(
                '* in totall, ',
                length(unique(all.groups)),
                ' homogeneous sample groups are created:'),
            color = 'blue',
            verbose = verbose
        )
        if(isTRUE(verbose))
            print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## use only categorical ####
    if (length(categorical.bio.var) > 0 & length(continuous.bio.var) == 0){
        printColoredMessage(
            message = '- Creating all possible combination of the groups using categorical variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- apply(
            categorical.bio.data,
            1,
            paste , collapse = "..")
        printColoredMessage(
            message = paste0(
                '*in totall, ',
                length(unique(all.groups)),
                ' homogeneous groups are created:'),
            color = 'blue',
            verbose = verbose)
        if(isTRUE(verbose))
            print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## error ####
    if(length(categorical.bio.var) == 0 & length(continuous.bio.var) == 0){
        stop('There is no categorical or continuous variables to create homogeneous groups.')
    }

    # Save the results ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(message = '-- Saving the results',
                        color = 'magenta',
                        verbose = verbose
                        )
    out.put.name <- paste0(
        length(unique(all.groups)),
        ' groups|',
        paste0(bio.variables, collapse = '&'),
        '|Clus:',
        clustering.method,
        '|nb:',
        nb.clusters
        )
    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Saving the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if (!'HomogeneousGroups' %in%  names(se.obj@metadata)) {
            se.obj@metadata[['HomogeneousGroups']] <- list()
        }
        if (!'BiologicalVariables' %in%  names(se.obj@metadata[['HomogeneousGroups']]) ) {
            se.obj@metadata[['HomogeneousGroups']][['BiologicalVariables']] <- list()
        }
        se.obj@metadata[['HomogeneousGroups']][['BiologicalVariables']][[out.put.name]] <- all.groups
        printColoredMessage(
            message = '- The homogeneous groups are saved to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The createHomogeneousBioGroups function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## save the output as a vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The homogeneous groups are outputed as a vector.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The createHomogeneousBioGroups function finished.',
                            color = 'white',
                            verbose = verbose
        )
        return(all.groups)
    }
}
