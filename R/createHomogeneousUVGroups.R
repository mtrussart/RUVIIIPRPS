#' Create all possible homogeneous groups with respect to unwanted variables.

#' @author Ramyar Molania

#' @description
#' This function generates all possible homogeneous sample groups based on the specified unwanted variables.

#' @details
#' The function generates all possible homogeneous sample groups based on the specified unwanted variables. If continuous
#' variables are provided, the function splits each into a number of clusters determined by ’nb.clusters’, using the
#' clustering method specified in ’clustering.method’. Ultimately, all combinations of all clusters are created and
#' each such combination is regarded as a homogeneous sample group concerning biological variables.
#' For example, if library size and time e.g. years from 2010-2014, are provided as sources of unwanted variation,
#' the createHomogeneousUVGroups() will divide the library size into n groups and then the homogeneous sample
#' groups respect to these variables will be combinations of one of 5 years * n groups based on library size. These
#' combinations are the sample groups which will be regarded as homogeneous for the known unwanted variation.
#'

#' @param se.obj A SummarizedExperiment object.
#' @param uv.variables Symbol. A symbol or a vector of symbols specifying the column names of unwanted variables in
#' the sample annotation of the SummarizedExperiment object. These 'uv.variables' can be either categorical or continuous
#' variables.
#' @param clustering.method Symbol. A symbol specifying the clustering method to be applied for grouping each continuous
#' source of unwanted variables. Options include 'kmeans', 'cut', and 'quantile'. The default is 'kmeans' clustering.
#' @param nb.clusters Numeric. A value indicating the number of groups for continuous sources of unwanted variation.
#' The default is 3. This implies that each continuous source will be split into 3 groups using the specified
#' 'clustering.method' method.
#' @param assess.variables Logical. Indicates whether to assess correlations between the unwanted variables. If 'TRUE',
#' the function 'assessVariablesCorrelation' will be applied. For more details refer to the 'assessVariablesCorrelation'
#' function.The default is 'TRUE'.
#' @param cat.cor.coef Vector. Containing two numerical values that specifies the cut-off for the correlation coefficient
#' between each pair of categorical variables. The first value applies to correlations between each pair of 'uv.variables',
#' and the second value applies to correlations between each pair of 'bio.variables'. Correlation is computed using the
#' ContCoef function from the DescTools R package. If the correlation between a pair of variables exceeds the cut-off,
#' only the variable with the highest number of factors will be retained, and the other will be excluded from further
#' analysis. By default, both values are set to 0.9.
#' @param cont.cor.coef Vector. Containing two numerical values that specifies the cut-off for the correlation coefficient
#' between each pair of continuous variables. The first value applies to correlations between each pair of continuous
#' 'uv.variables', and the second value applies to correlations between each pair of continuous 'bio.variables'.
#' Correlation is computed using the ContCoef function from the DescTools R package. If the correlation between a pair of
#' variables exceeds the cut-off, only the variable with the highest variance will be retained, and the other will be
#' excluded from further analysis. By default, both values are set to 0.9.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeobj' will be applied. The default is 'TRUE'.
#' @param remove.na Symbol. Indicates whether to remove missing values from the 'uv.variables'. The options are
#' 'sample.annotation' or 'none'. The default is 'sample.annotation', indicating the missing values from the variables
#' will be removed.
#' @param save.se.obj Logical. Indicates whether to save the results to the metadata of the SummarizedExperiment object
#' or not. If 'TRUE', all the possible homogeneous groups will be saved into "metadata$HGgroups$UVgroups", otherwise
#' the results will outputted as a vector. The default is 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object containing the all possible homogeneous groups in the "metadata$HGgroups$UVgroups"
#' or a vector of the possible homogeneous groups.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable
#' @export

createHomogeneousUVGroups <- function(
        se.obj,
        uv.variables,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        assess.variables = TRUE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The createHomogeneousUVGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check  the inputs ####
    if(is.null(uv.variables)){
        stop('The uv.variables cannot be empty.')
    } else if(!is.vector(uv.variables)){
        stop('The "uv.variables" must be a vector.')
    } else if(!clustering.method %in% c('kmeans', 'cut', 'quantile')){
        stop('The clustering.method should be one of: kmeans, cut or quantile.')
    } else if(max(cat.cor.coef) > 1){
        stop('The maximum value for the cat.cor.coef is 1.')
    } else if(max(cont.cor.coef) > 1){
        stop('The maximum value for the cont.cor.coef is 1.')
    } else if( remove.na %in% c('both', 'measuerments')){
        stop('The remove.na should be either sample.annotation or none.')
    }

    # Assess association between variables ####
    if (isTRUE(assess.variables)) {
        printColoredMessage(
            message = '-- Assess the correlation between unwanted variation variables:',
            color = 'magenta',
            verbose = verbose)
        se.obj <- assessVariablesAssociation(
            se.obj = se.obj,
            bio.variables = NULL,
            uv.variables = uv.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = TRUE,
            remove.na = remove.na,
            verbose = verbose)
        uv.variables <- se.obj$uv.variables
        se.obj <- se.obj$se.obj
    }

    # Find the class of unwanted variables ####
    printColoredMessage(
        message = '-- Find the class and levels of the variable(s):',
        color = 'magenta',
        verbose = verbose
    )
    class.uv.var <- unlist(lapply(
        uv.variables,
        function(x) class(se.obj[[x]])))
    categorical.uv <- uv.variables[class.uv.var %in% c('factor', 'character')]
    if(length(categorical.uv) > 1){
        printColoredMessage(
            message = paste0('- ', length(categorical.uv), ' categorical variables are provided.'),
            color = 'blue',
            verbose = verbose
            )
        categorical.uv <- unlist(lapply(
            categorical.uv,
            function(x){
                if(length(unique(colData(se.obj))[[x]]) == 1){
                    printColoredMessage(
                        message = paste0('- the "', x, '" variable contains 1 levels. Then, this will be excluded.'),
                        color = 'blue',
                        verbose = verbose)
                } else {
                    printColoredMessage(
                        message = paste0('- the "', x, '" variable contains ',  length(unique(colData(se.obj)[[x]])), ' levels.'),
                        color = 'blue',
                        verbose = verbose)
                    keep.variable <- x
                }
                return(keep.variable)
            }))
    }
    continuous.uv <- uv.variables[class.uv.var %in% c('numeric', 'integer')]
    if(length(continuous.uv) > 0){
        printColoredMessage(
            message = paste0('- ', length(categorical.uv), ' continuous variables are provided.'),
            color = 'blue',
            verbose = verbose
        )
        continuous.uv <- unlist(lapply(
            continuous.uv,
            function(x){
                if(var(colData(se.obj)[[x]]) == 0){
                    printColoredMessage(
                        message = paste0('- the variance of the "', x, '" variable is 0 . Then, this will be excluded.'),
                        color = 'blue',
                        verbose = verbose
                        )
                } else {
                    printColoredMessage(
                        message = paste0('- the variance of the "', x, '" variable is ~ ',  round(x = var(colData(se.obj)[[x]]), digits = 5) , '.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    keep.variable <- x
                }
                return(keep.variable)
            }))
    }

    # Cluster continuous variable ####
    if(length(continuous.uv) > 0){
        ## check clustering value for continuous variable
        if(nb.clusters == 1){
            stop('The value of "nb.clusters" must be bigger than 1.')
        } else if (is.null(nb.clusters)){
            stop('The "nb.clusters" cannot be empty.')
        } else if ( nb.clusters == 0) {
            stop('The value of "nb.clusters" cannot be 0.')
        }
        ## cluster continuous variables ####
        printColoredMessage(
            message = '-- Cluster continuous sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose
            )
            printColoredMessage(message = paste0(
                '- each source will be divided into ', nb.clusters, ' groups using the ', clustering.method, ' method.'),
                color = 'blue',
                verbose = verbose
                )
            ## kmeans ####
            if (clustering.method == 'kmeans') {
                    set.seed(3456)
                    continuous.uv.groups <- lapply(
                        continuous.uv,
                        function(x){
                            uv.cont.clusters <- stats::kmeans(
                                x = colData(se.obj)[[x]],
                                centers = nb.clusters,
                                iter.max = 1000)
                            paste0(x, '_group', uv.cont.clusters$cluster)
                        })
                    names(continuous.uv.groups) <- continuous.uv
                    continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                    continuous.uv.groups <- apply(
                        continuous.uv.groups,
                        1,
                        paste , collapse = "..")
            }
            ## cut ####
            if (clustering.method == 'cut') {
                continuous.uv.groups <- lapply(
                    continuous.uv,
                    function(x) {
                        uv.cont.clusters <- as.numeric(
                            cut(x = colData(se.obj)[[x]],
                                breaks = nb.clusters,
                                include.lowest = TRUE))
                        paste0(x, '_group', uv.cont.clusters)
                    })
                names(continuous.uv.groups) <- continuous.uv
                continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                continuous.uv.groups <- apply(
                    continuous.uv.groups,
                    1,
                    paste , collapse = "..")
            }
            ## quantile ####
            if (clustering.method == 'quantile') {
                continuous.uv.groups <- lapply(
                    continuous.uv,
                    function(x) {
                        quantiles <- quantile(x = colData(se.obj)[[x]], probs = seq(0, 1, 1 / nb.clusters))
                        uv.cont.clusters <- as.numeric(
                            cut(x = colData(se.obj)[[x]],
                                breaks = quantiles,
                                include.lowest = TRUE))
                        paste0(x, '_group', uv.cont.clusters)
                    })
                names(continuous.uv.groups) <- continuous.uv
                continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                continuous.uv.groups <- apply(
                    continuous.uv.groups,
                    1,
                    paste , collapse = "..")
            }
    }

    # Check categorical variable ####
    if(length(categorical.uv) > 0){
        categorical.uv.groups <- as.data.frame(colData(se.obj)[, categorical.uv, drop = FALSE])
    }
    # Create all possible homogeneous groups ####
    printColoredMessage(
        message = '-- Create all possible combination of the groups:',
        color = 'magenta',
        verbose =  verbose
        )
    ## use both continuous and categorical ####
    if(length(categorical.uv) > 0 & length(continuous.uv) > 0 ){
        printColoredMessage(
            message = '- create all possible combination of the groups using both continuous and categorical variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- categorical.uv.groups
        all.groups$cont <- continuous.uv.groups
        all.groups <- unname(apply(
            all.groups,
            1,
            paste , collapse = ".."
        ))
        printColoredMessage(
            message = paste0(
                '- in totall, ',
                length(unique(all.groups)),
                ' homogeneous groups  are created:'),
            color = 'blue',
            verbose = verbose
        )
        cat("", file = stdout(), sep = "")
        if(isTRUE(verbose))
            print(x = kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## use only continuous ####
    if (length(categorical.uv) == 0 & length(continuous.uv) > 0){
        printColoredMessage(
            message = '- create all possible combination of the groups using continuous variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- continuous.uv.groups
        printColoredMessage(
            message = paste0(
                '- in totall, ',
                length(unique(all.groups)),
                ' homogeneous groups are created:'),
            color = 'blue',
            verbose = verbose
            )
        if(isTRUE(verbose))
            print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## use only categorical ####
    if(length(categorical.uv) > 0 & length(continuous.uv) == 0){
        printColoredMessage(
            message = '- create all possible combination of the groups using categorical variables:',
            color = 'blue',
            verbose =  verbose
        )
        all.groups <- apply(
            categorical.uv.groups,
            1,
            paste , collapse = "..")
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous groups are created:'),
            color = 'blue',
            verbose = verbose)
        if(isTRUE(verbose))
            print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    ## error ####
    if(length(categorical.uv) == 0 & length(continuous.uv) == 0){
        stop('There is no categorical or continuous variables to create homogeneous groups.')
    }

    # Save the results ####
    ## add results to the SummarizedExperiment object ####
    out.put.name <- paste0(
        length(unique(all.groups)),
        ' groups|',
        paste0(uv.variables, collapse = '&'),
        '|Clus:',
        clustering.method,
        '|nb:',
        nb.clusters
        )
    printColoredMessage(message = '-- Save the results',
        color = 'magenta',
        verbose = verbose)
    if(isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- save the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if(!'HomogeneousGroups' %in%  names(se.obj@metadata)) {
            se.obj@metadata[['HomogeneousGroups']] <- list()
        }
        if(!'UnwantedVariables' %in%  names(se.obj@metadata[['HomogeneousGroups']]) ) {
            se.obj@metadata[['HomogeneousGroups']][['UnwantedVariables']] <- list()
        }
        se.obj@metadata[['HomogeneousGroups']][['UnwantedVariables']][[out.put.name]] <- all.groups
        printColoredMessage(
            message = '- The homogeneous groups are saved to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The createHomogeneousUVGroups function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    }
    ## save the output as a vector ####
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The homogeneous groups are outputed as a vector.',
            color = 'blue',
            verbose = verbose
            )
       printColoredMessage(message = '------------The createHomogeneousUVGroups function finished.',
                           color = 'white',
                           verbose = verbose
                           )
    return(all.groups)
    }
}
