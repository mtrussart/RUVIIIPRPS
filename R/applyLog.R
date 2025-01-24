#' Apply log2 on assay(s) + pseudo.count in SummarizedExperiment object.

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object. Thses assays will be log2 + pseudo.count transformed. The default is set to "all, which
#' indicates all the assays of the SummarizedExperiment object will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param assessment Symbol. Indicates the downstream analyses impacted by not applying log transformation to the data.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return logged data

applyLog <- function(
        se.obj,
        assay.names = 'all',
        pseudo.count = 1,
        assessment,
        verbose = TRUE
        ){
    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Apply logs ####
    all.assays.loged <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (!is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('- Apply log2 on the "', x, '" + ', pseudo.count, ' (pseudo.count) data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('- Apply log2 on the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                expr <- log2(assay(x = se.obj, i = x))
            }
            return(expr)
        })
    names(all.assays.loged) <- levels(assay.names)
    return(all.assays.loged)
}


