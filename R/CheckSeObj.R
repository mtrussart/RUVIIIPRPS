#' Assess SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function assesses the structure of a SummarizedExperiment object and removes any missing values from both
#' assay(s) and sample annotation. When multiple assays are provided, if they are missing in only one of the assays,
#' the corresponding rows in other assays will be remove as well. Please note that, the current RUV-III-PRPS method does
#' not support missing values in the assay(s).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols to specify the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is "all," indicating that all assays in the SummarizedExperiment object will
#' be assessed.
#' @param variables Symbol. A symbol or a vector of symbols specifying the name(s) of the column(s) in the sample
#' annotation of the SummarizedExperiment object. By default, it is set to 'all', then all the columns will be examined.
#' We recommend specifying those column(s) that are of interest to your analysis.
#' @param remove.na Symbol. A Symbol that specifies whether to eliminate missing values from either 'assays', 'sample.annotation',
# 'both', or 'none'. When 'assays' is chosen, genes containing missing values will be omitted. If 'sample.annotation'
# is selected, samples with NA or missing values for each 'variables' will be excluded. The default is 'both'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object.

#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom stats complete.cases
#' @export


checkSeObj <- function(
        se.obj,
        assay.names = 'all',
        variables = 'all',
        remove.na = 'both',
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The checkSeObj function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check function inputs ####
    if (remove.na =='both'){
        if(is.null(assay.names) | is.null(variables))
            stop('The "assay.names" or "variables" cannot be empty when the remove.na = "both".')
    } else if (remove.na == 'assays'){
        if(is.null(assay.names))
            stop('The "assay.names" cannot be empty when the remove.na = "assays".')
    } else if (remove.na == 'sample.annotation'){
        if(is.null(variables))
        stop('The "variables" cannot be empty when the remove.na = "sample.annotation".')
    }
    if (!is.null(remove.na) & !remove.na %in% c("assays", "sample.annotation" , "both", "none")) {
        stop('The "remove.na" must be one of the "assays", "sample.annotation", "both" or "none".')
    }
    if(is.list(variables)){
        stop('The "variables" must be a vector of the column names or variables = "all".')
    }
    if (length(variables) == 1 && variables == 'all') {
        variables <- colnames(colData(se.obj))
    }

    # Check the class the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Check the class, structure and dimention of the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
        )
    if (class(se.obj)[1] %in% c('SummarizedExperiment', 'RangedSummarizedExperiment')) {
        printColoredMessage(
            message = paste0('- The "se.obj" is a SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0('- The object contains: ', length(assays(se.obj)), ' data (assays).'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0('- The object contains: ', ncol(se.obj), ' (samples) and ', nrow(se.obj),' (genes / measurements).'),
            color = 'blue',
            verbose = verbose
            )

    } else stop('The "se.obj" provided is not a class of the SummarizedExperiment object.')

    # Check the assays ####
    ## check whether the provided assays exist ####
    if(!is.null(assay.names)){
        printColoredMessage(
            message = '-- Check the assay name(s):',
            color = 'magenta',
            verbose = verbose
            )
        if(is.list(assay.names))
            stop('The "assay.names" must be a vector of the assay names or assay.names = "all".')
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names, levels = assay.names)

        if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            not.found.assay <- assay.names[!assay.names %in% names(assays(se.obj)) ]
            stop(paste0(
                'The assay(s): ',
                paste(not.found.assay, collapse = ', '),
                ' cannot be found in the SummarizedExperiment object.'))
        } else {
            printColoredMessage(
                message = paste0(
                    'All the assay(s): ',
                    paste0(assay.names, collapse = ', '),
                    ' are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
    }
    # Checking the variables exist ####
    if(!is.null(variables)){
        printColoredMessage(
            message = '-- Check the variable name(s):',
            color = 'magenta',
            verbose = verbose
            )
        if(!sum(variables %in% colnames(colData(se.obj))) == length(variables)){
            not.found.variables <- variables[!variables %in% colnames(colData(se.obj))]
            stop(paste0(
                'The variable(s) ',
                paste(not.found.variables, collapse = ', '),
                ' cannot be found in the SummarizedExperiment object.'))
        } else{
            printColoredMessage(
                message = paste0(
                    'All the variable(s): ',
                    paste0(variables, collapse = ', '),
                    ' are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
    }
    # Find missing/NA values ####
    printColoredMessage(
        message = '-- Find missing or NA values:',
        color = 'magenta',
        verbose = verbose
        )
    ## variables ####
    if(!is.null(variables)){
        printColoredMessage(
            message = '- Find missing/NA values in the variable(s):',
            color = 'orange',
            verbose = verbose)
        variables.with.na <- sapply(
            variables,
            function(x) sum(is.na(colData(se.obj)[[x]])))
        if (sum(variables.with.na > 0) > 0) {
            printColoredMessage(
                message = '* The variable(s) with missing values status.',
                color = 'blue',
                verbose = verbose)
            if(isTRUE(verbose))
                print(knitr::kable(x = variables.with.na, col.names = 'number of NA'))
            } else printColoredMessage(
                message = '* 0 missing or NA values are found in the variable(s).',
                color = 'blue',
                verbose = verbose
            )
        } else if (is.null(variables))  printColoredMessage(
            message = '* Any variable(s) are specified to find missing or NA values.',
            color = 'blue',
            verbose = verbose
            )
    ## datasets (assays) ####
    if (!is.null(assay.names)) {
        ## checking na in the assays ####
        printColoredMessage(
            message = '- Find the genes (measurements) with missing value:',
            color = 'orange',
            verbose = verbose)
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        colnames(measurements.with.na) <- assay.names
        if (sum(measurements.with.na) > 0) {
            printColoredMessage(
                message = 'The assay(s) with missing values status.',
                color = 'blue',
                verbose = verbose)
            if(isTRUE(verbose))
                knitr::kable(
                    x = colSums(measurements.with.na),
                    col.names = 'Number of genes (measurements) with NA')
            printColoredMessage(
                message = 'Note, the current RUV-III-PRPS package does not support NA values.',
                color = 'red',
                verbose = verbose)
        } else if (sum(measurements.with.na) == 0) printColoredMessage(
            message = '* 0 missing or NA values are found in the dataset(s)(assays).',
            color = 'blue',
            verbose = verbose
            )
    } else if (is.null(assay.names)) {
        printColoredMessage(
            message = '* Any dataset(s)(assays) are specified to find missing or NA values.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = 'Note, the current RUV-III method does not support missing values in the data.',
            color = 'red',
            verbose = verbose
            )
    }
    # Remove missing or NA values ####
    printColoredMessage(
        message = '-- Remove missing or NA values:',
        color = 'magenta',
        verbose = verbose
        )
    ## both assays and sample annotation ####
    if(remove.na == 'both'){
        printColoredMessage(
            message = '- Remove missing or NA values from variable(s):',
            color = 'orange',
            verbose = verbose
        )
        ### variables ####
        keep.samples <- complete.cases(as.data.frame(colData(se.obj)[, variables, drop = FALSE]))
        if(sum(keep.samples) < ncol(se.obj)){
            printColoredMessage(
                message = paste0(
                    sum(!keep.samples),
                    ' samples are removed due to the presence of missing or NA values in the variable(s).'),
                color = 'blue',
                verbose = verbose
                )
            se.obj <- se.obj[ , keep.samples]
        } else if(sum(keep.samples) == ncol(se.obj)){
            printColoredMessage(
                message = '* 0 samples are removed due to the presence of missing or NA values in the variable(s).',
                color = 'blue',
                verbose = verbose
            )
        } else if (sum(keep.samples) == 0){
            stop('All samples are removed due to the presence of missing or NA values in the variable(s)')
        }
        printColoredMessage(
            message = '- Remove missing or NA values from dataset(s):',
            color = 'orange',
            verbose = verbose
        )
        ## datasets ####
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        if(sum(measurements.with.na) > 0){
            printColoredMessage(
                message = paste0(
                    sum(rowSums(measurements.with.na)!=0),
                    ' gene(s) are removed due to the presence of missing or NA values in the dataset(s).'),
                color = 'blue',
                verbose = verbose
                )
            se.obj <- se.obj[rowSums(measurements.with.na) == 0, ]
        } else if (sum(measurements.with.na) ==  0){
            printColoredMessage(
                message =  ' 0 gene(s) are removed due to the presence of missing or NA values in the dataset(s).',
                color = 'blue',
                verbose = verbose
            )
        } else if (sum(measurements.with.na) ==  nrow(se.obj)){
            stop('All genes are removed due to the presence of missing or NA values in the variable(s)')
        }

    }
    ## datasets (assays) ####
    if (remove.na == 'assays'){
        printColoredMessage(
            message = '- Remove missing or NA values from dataset(s):',
            color = 'orange',
            verbose = verbose
        )
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        if(sum(measurements.with.na) > 0){
            printColoredMessage(
                message = paste0(
                    sum(rowSums(measurements.with.na)!=0),
                    ' gene(s) are removed due to the presence of missing or NA values in the dataset(s).'),
                color = 'blue',
                verbose = verbose
            )
            se.obj <- se.obj[rowSums(measurements.with.na) == 0, ]
        } else if (sum(measurements.with.na) ==  0){
            printColoredMessage(
                message =  ' 0 gene(s) are removed due to the presence of missing or NA values in the dataset(s).',
                color = 'blue',
                verbose = verbose
            )
        } else if (sum(measurements.with.na) ==  nrow(se.obj)){
            stop('All genes are removed due to the presence of missing or NA values in the variable(s)')
        }
    }
    ## variables ####
    if (remove.na == 'sample.annotation'){
        printColoredMessage(
            message = '- Remove missing or NA values from variable(s):',
            color = 'orange',
            verbose = verbose
        )
        keep.samples <- complete.cases(as.data.frame(colData(se.obj)[, variables, drop = FALSE]))
        if(sum(keep.samples) < ncol(se.obj)){
            printColoredMessage(
                message = paste0(
                    sum(!keep.samples),
                    ' samples are removed due to the presence of missing or NA values in the variable(s).'),
                color = 'blue',
                verbose = verbose
            )
            se.obj <- se.obj[ , keep.samples]
        } else if(sum(keep.samples) == ncol(se.obj)){
            printColoredMessage(
                message = '* 0 samples are removed due to the presence of missing or NA values in the variable(s).',
                color = 'blue',
                verbose = verbose
            )
        } else if (sum(keep.samples) == 0){
            stop('All samples are removed due to the presence of missing or NA values in the variable(s)')
        }
    }
    ## none ####
    if (remove.na == 'none'){
        printColoredMessage(
            message = 'Any missing or NA values from both assays and variable will not be removed.',
            color = 'red',
            verbose = verbose
            )
        printColoredMessage(
            message = 'Note, the current RUV-III-PRPS package does not support NA values.',
            color = 'red',
            verbose = verbose
        )
    }
    # Return the ####
    printColoredMessage(
        message = '-- Retrun final SummarizedExperiment:',
        color = 'magenta',
        verbose = verbose
    )
    printColoredMessage(
        message = paste0(
            '- The final SummarizedExperiment object contains: ',
            length(assays(se.obj)), ' data (assays).'),
        color = 'blue',
        verbose = verbose
    )
    printColoredMessage(
        message = paste0(
            '- The final SummarizedExperiment object contains: ',
            ncol(se.obj), ' (samples) and ', nrow(se.obj),' (genes / measurements).'),
        color = 'blue',
        verbose = verbose
    )
    printColoredMessage(message = '------------The checkSeObj function finished.',
                        color = 'white',
                        verbose = verbose)

    return(se.obj)
}
