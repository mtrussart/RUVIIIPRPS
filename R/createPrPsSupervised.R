#' Create supervised pseudo-replicates of pseudo samples (PRPS).

#' @author Ramyar Molania

#' @description
#' This function creates different PRPS sets when all sources of unwanted and biological variation are known.

#' @details
#' We will create distinct group of pseudo-replicates for each source of unwanted variation defined in the 'uv.variables'
#' argument. For example to correct for batch effect if defined in the 'uv.variables' argument, several group of
#' pseudo-samples will be created by averaging the samples of the same biological subtype defined in 'bio.variables' in
#' each batch. Then those pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set. For example
#' to correct for library size if defined in the 'uv.variables' argument, several group of pseudo-samples will be created
#' by averaging the top and bottom-ranked samples by library size of the same biological subtype in each batch. Then those
#' pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set. Similarly to correct for purity if
#' defined in the uv.variables' argument, several group of pseudo-samples will be created by averaging the top and
#' bottom-ranked samples by purity of the same biological subtype in each batch. Then those pseudo-samples will be defined
#' as pseudo-replicates which constitutes a PRPS set.

#' @references
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol used to specify the assay name within the SummarizedExperiment object. The selected
#' assay should be the one that will be used for RUV-III-PRPS normalization.
#' @param bio.variables Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
#' subtypes, tumour purity, ... within the SummarizedExperiment object. This can comprise a vector containing either
#' categorical, continuous, or a combination of both variables.
#' @param uv.variables Symbol. A symbol or symbols representing the label of unwanted variable(s), such as batch effects,
#' library size, ... within the SummarizedExperiment object. This can comprise a vector containing either categorical,
#' continuous, or a combination of both variables.
#' @param apply.other.uv.variables Logical. Determines whether to include other specified unwanted variables when generating
#' PRPS sets for individual ones. The default is set to 'TRUE'.
#' @param min.sample.for.prps Numeric. Indicates the minimum number of biologically homogeneous samples to be averaged
#' to create one pseudo-sample. The default it is set to 3.
#' @param bio.clustering.method Symbol. A symbol indicating the clustering method used to group each continuous biological
#' variable. The option include 'kmeans', cut' and 'quantile'. The default is to 'kmeans'. We refer to the
#' createHomogeneousBioGroups() function for more details.
#' @param nb.bio.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous biological
#' variable. The default it set to 3.
#' @param other.uv.clustering.method Symbol. A symbol indicating the clustering method used to group each continuous unwanted
#' variable. The option include 'kmeans', cut' and 'quantile'. The default is to 'kmeans'. We refer to the
#' createHomogeneousUVGroups() function for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous biological
#' variable. The default it set to 3.
#' @param check.prps.connectedness Logical. Indicates whether to assess the connectedness between the PRPS sets or not.
#' The default is set to 'TRUE'. See the details for more information.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. The default is set to 'TRUE'.
#' @param pseudo.count Numeric. A numeric value as a pseudo count to be added to all measurements before log transformation.
#' The default is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. If 'TRUE', the checkSeObj
#' will be applied inside the function. The default is set to 'TRUE'.
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
#' @param remove.na Symbol. A symbol Indicating whether to remove NA or missing values from either the 'assays', the
#' 'sample.annotation', both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be
#' excluded. If sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and
#' 'uv.variables' will be excluded. The default is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the results in the metadata of the SummarizedExperiment object
#' or to output the result as list. The default by is set to 'TRUE'.
#' @param plot.output Logical. Indicates whether to generate the PRPS map plot for individual sources of unwanted variation.
#' The default is 'TRUE'.
#' @param prps.name Symbol. A symbol to specify the name of all PRPS sets that will be created for all specified source(s)
#' of unwanted variation. The default is set to 'NULL'. The, the function creates a name based on paste0('prps_', uv.variable).
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object that contains all the PRPS data and PPRS map plot in the metadata or a list
#' that contains all the results.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @export

createPrPsSupervised <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        apply.other.uv.variables = TRUE,
        min.sample.for.prps = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        check.prps.connectedness = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output = TRUE,
        prps.name = NULL,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The createPrPsSupervised function starts:',
                        color = 'white',
                        verbose = verbose)
    # define categorical and continuous variables ####
    uv.class <- sapply(uv.variables,
                       function(x) class(colData(se.obj)[[x]]))
    categorical.uv <- names(uv.class[which(uv.class %in% c('character', 'factor'))])
    continuous.uv <- uv.variables[!uv.variables %in% categorical.uv]

    # prps for categorical variables ####
    if (length(categorical.uv) > 0) {
        printColoredMessage(
            message = '-- Create PRPS for all categorical sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose
        )
        if (!save.se.obj) {
            categorical.uv.prps <- lapply(
                categorical.uv,
                function(x) {
                    if (apply.other.uv.variables) {
                        other.uv.variables <- uv.variables[!uv.variables %in% x]
                    } else other.uv.variables <- NULL
                    createPrPsForCategoricalUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        bio.variables = bio.variables,
                        main.uv.variable = x,
                        other.uv.variables = other.uv.variables,
                        min.sample.for.prps = min.sample.for.prps,
                        bio.clustering.method = bio.clustering.method,
                        nb.bio.clusters = nb.bio.clusters,
                        other.uv.clustering.method = other.uv.clustering.method,
                        nb.other.uv.clusters = nb.other.uv.clusters,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        check.prps.connectedness = check.prps.connectedness,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        remove.na = remove.na,
                        prps.name = prps.name,
                        verbose = verbose
                    )
                })
            names(categorical.uv.prps) <- categorical.uv
            categorical.uv.prps.all <- do.call(cbind, categorical.uv.prps)
        } else {
            for (x in categorical.uv) {
                if (apply.other.uv.variables) {
                    other.uv.variables <- uv.variables[!uv.variables %in% x]
                } else other.uv.variables <- NULL
                se.obj <- createPrPsForCategoricalUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    bio.variables = bio.variables,
                    main.uv.variable = x,
                    other.uv.variables = other.uv.variables,
                    min.sample.for.prps = min.sample.for.prps,
                    bio.clustering.method = bio.clustering.method,
                    nb.bio.clusters = nb.bio.clusters,
                    other.uv.clustering.method = other.uv.clustering.method,
                    nb.other.uv.clusters = nb.other.uv.clusters,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    cat.cor.coef = cat.cor.coef,
                    cont.cor.coef = cont.cor.coef,
                    check.prps.connectedness = check.prps.connectedness,
                    assess.se.obj = FALSE,
                    save.se.obj = save.se.obj,
                    remove.na = remove.na,
                    prps.name = prps.name,
                    verbose = verbose
                )
            }
        }
    }
    if (length(continuous.uv) > 0) {
        printColoredMessage(message = '-- Create PRPS for all continuous sources of unwanted variation:',
                            color = 'magenta',
                            verbose = verbose)
        if (!save.se.obj) {
            continuous.uv.prps <- lapply(
                continuous.uv,
                function(x) {
                    if (apply.other.uv.variables) {
                        other.uv.variables <- uv.variables[!uv.variables %in% x]
                    } else other.uv.variables <- NULL
                    createPrPsForContinuousUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        bio.variables = bio.variables,
                        main.uv.variable = x,
                        other.uv.variables = other.uv.variables,
                        min.sample.for.prps = min.sample.for.prps,
                        bio.clustering.method = bio.clustering.method,
                        nb.bio.clusters = nb.bio.clusters,
                        other.uv.clustering.method = other.uv.clustering.method,
                        nb.other.uv.clusters = nb.other.uv.clusters,
                        cat.cor.coef = cat.cor.coef,
                        cont.cor.coef = cont.cor.coef,
                        apply.log = apply.log,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        pseudo.count = pseudo.count,
                        remove.na = remove.na,
                        prps.name = prps.name,
                        verbose = verbose
                    )
                })
            names(continuous.uv.prps) <- continuous.uv
            continuous.uv.prps.all <- do.call(cbind, continuous.uv.prps)
        } else{
            for (x in continuous.uv) {
                if (apply.other.uv.variables) {
                    other.uv.variables <- uv.variables[!uv.variables %in% x]
                } else other.uv.variables <- NULL
                se.obj <- createPrPsForContinuousUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    bio.variables = bio.variables,
                    main.uv.variable = x,
                    other.uv.variables = other.uv.variables,
                    min.sample.for.prps = min.sample.for.prps,
                    bio.clustering.method = bio.clustering.method,
                    nb.bio.clusters = nb.bio.clusters,
                    other.uv.clustering.method = other.uv.clustering.method,
                    nb.other.uv.clusters = nb.other.uv.clusters,
                    cat.cor.coef = cat.cor.coef,
                    cont.cor.coef = cont.cor.coef,
                    apply.log = apply.log,
                    assess.se.obj = FALSE,
                    save.se.obj = save.se.obj,
                    pseudo.count = pseudo.count,
                    remove.na = remove.na,
                    prps.name = prps.name,
                    verbose = verbose
                )
            }
        }
    }
    # save the output ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsSupervised function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The createPrPsSupervised function finished.',
                            color = 'white',
                            verbose = verbose)
        if (length(continuous.uv) > 0 &
            length(categorical.uv) > 0) {
            return(
                list(
                    categorical.uv.prps = categorical.uv.prps,
                    continuous.uv.prps = continuous.uv.prps,
                    all.prps = cbind(categorical.uv.prps.all, continuous.uv.prps.all))
            )
        } else if (length(continuous.uv) > 0 &
                   length(categorical.uv) == 0) {
            return(
                list(continuous.uv.prps = continuous.uv.prps,
                     all.prps = continuous.uv.prps.all))
        } else if (length(continuous.uv) == 0 &
                   length(categorical.uv) > 0) {
            return(
                list(categorical.uv.prps = categorical.uv.prps,
                     all.prps = categorical.uv.prps.all))
        }
    }

}
