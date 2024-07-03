#' Assess the different sets of NCGs.

#' @author Ramyar Molania

#' @description
#' This function assesses the performance of various sets of negative control genes in terms of their ability to capture
#' unwanted variation and their lack of association with biological variation. First, for each individual set of negative
#' control genes (NCG), the function applies PCA to the specified data. Then, it calculates the association of the 'nb.pcs'
#' principal components (PCs) with both specified biological and unwanted variables. The correlation values obtained for
#' the biological variables will be subtracted from 1. Finally, the average of all these correlations is computed to
#' determine the final correlation values. A higher correlation indicates better performance.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol indicating the name of the assay in the SummarizedExperiment object. This assay
#' should be the one that will be used as an input data for the RUV-III-PRPS normalization.
#' @param variables Symbol. A symbol or a set of symbols indicating the names of the variables in the SummarizedExperiment
#' object. These will be used to assess the performance of NCGS.
#' @param bio.variables Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
#' subtypes, tumor purity, ... within the SummarizedExperiment object. This can comprise a vector containing either
#' categorical, continuous, or a combination of both variables.
#' @param uv.variables Symbol. A symbol or symbols representing the label of unwanted variable(s), such as cancer
#' batch effects, library size, ... within the SummarizedExperiment object. This can comprise a vector containing either
#' categorical, continuous, or a combination of both variables.
#' @param ncg.type Symbol. A symbol indicating the type of NCGs that are obtained either by supervised or un.supervised
#'  approaches. The options are 'supervised' or 'un.supervised'.
#' @param ncg.group Symbol. A symbol specifying the name of a group of NCGs under the supervised or un.supervised slots
#' in the SummarizedExperiment object.
#' @param ncg.set.names A symbol specifying th exact name of a NCG set under the 'ncg.group' name in the supervised or
#' un.supervised slots in the SummarizedExperiment object. The default is set to 'all'. This will select all avaiable
#' NCGs sets.
#' @param nb.pcs Numeric. A numeric value indicating the number of PCs to be computed and used for the assessment analysis.
#' The default is set to 5.
#' @param center Logical. Indicates whether to scale the data or not prior to apply SVD. If center is 'TRUE', then
#' centering is performed by subtracting the column means of the assay from their corresponding columns. The default is 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is 'TRUE', then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is 'TRUE', and the root
#' mean square otherwise. The default is FALSE.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before computing the SVD. The
#' default is 'TRUE'. The data must be in log scale before computing the SVD.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before applying log transformation.
#' The default is 1. This argument cannot be NULL or negative.
#' @param plot.output Logical. Indicates whether to plot the assessment plot or not while running the function. The default
#' is set to 'TRUE'.
#' @param output.name Symbol. A symbol indicating the name under which the assessment plot will be saved in the
#' SummarizedExperiment object. If 'NULL', the function will select the name as follow:
#' 'paste0('comparsion_', paste0(names(ncgs), collapse = '&'))'
#' @param save.se.obj Logical. Indicates whether save the assessment plot to the SummarizedExperiment object or not. The
#' default is set to 'TRUE'. The plot will be save to:
#' @param verbose Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.

#' @importFrom dplyr bind_rows summarise mutate group_by
#' @importFrom fastDummies dummy_cols
#' @importFrom gtools mixedorder
#' @importFrom tidyr pivot_longer
#' @export

assessNCGs <- function(
        se.obj ,
        assay.name,
        variables ,
        bio.variables = NULL,
        uv.variables = NULL,
        ncg.type,
        ncg.group,
        ncg.set.names = 'all',
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        plot.output = TRUE,
        output.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Selection of NCGs ####
    if(ncg.type == 'supervised'){
        if(ncg.set.names == 'all'){
            ncgs <- se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set
        } else {
            if(sum(ncg.set.names %in% names(se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set)) != length(ncg.set.names) ){
                stop('All or some of the "ncg.set.names" cannot be found in the SummarizedExperiment.')
            } else{
                ncgs <- lapply(
                    ncg.set.names,
                    function(x){
                        se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set[[x]]
                    })
                names(ncgs) <- ncg.set.names
            }
        }
    }
    if(ncg.type == 'un.supervised'){
        if(ncg.set.names == 'all'){
            ncgs <- se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set
        } else {
            if(sum(ncg.set.names %in% names(se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set)) != length(ncg.set.names) ){
                stop('All or some of the "ncg.set.names" cannot be found in the SummarizedExperiment.')
            } else{
                ncgs <- lapply(
                    ncg.set.names,
                    function(x){
                        se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set[[x]]
                    })
                names(ncgs) <- ncg.set.names
            }
        }
    }

    # Perform PCA ####
    all.pcs.pn.ncgs <- lapply(
        names(ncgs),
        function(x){
            pcs <- computePCA(
                se.obj = se.obj[ncgs[[x]] , ],
                assay.names = assay.name,
                fast.pca = TRUE,
                nb.pcs = nb.pcs,
                center = center,
                scale = scale,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                svd.bsparam = bsparam(),
                assess.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE
            )
            pcs[[assay.name]]$svd$u
        })
    names(all.pcs.pn.ncgs) <- names(ncgs)

    ## variables class ####
    all.vars <- c(bio.variables, uv.variables)
    class.all.vars <- sapply(
        all.vars,
        function(x) class(colData(x = se.obj)[[x]])
    )
    cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
    cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])

    ## linear regression ####
    cont.vars.r.squareds <- lapply(
        names(ncgs),
        function(x) {
            pcs <- all.pcs.pn.ncgs[[x]]
            r.squareds <- sapply(
                cont.vars,
                function(y) {
                    lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ pcs))
                    lm.reg$r.squared
                })
            names(r.squareds) <- cont.vars
            r.squareds
        })
    names(cont.vars.r.squareds) <- names(ncgs)
    cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)

    ## vector correlation ####
    cat.vars.vec.corr <- lapply(
        names(ncgs),
        function(x) {
            pcs <- all.pcs.pn.ncgs[[x]]
            vec.corr <- sapply(
                cat.vars,
                function(y) {
                    catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca <- cancor(x = pcs, y = catvar.dummies)
                    1 - prod(1 - cca$cor ^ 2)
                })
            names(vec.corr) <- cat.vars
            vec.corr
        })
    names(cat.vars.vec.corr) <- names(ncgs)
    cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)

    ### put all together ####
    groups <- corr <- data <- assess <- variable <- NULL
    all <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
        round(digits = 3)
    all <- mutate(.data = all, var = row.names(all)) %>%
        pivot_longer(cols = -var, values_to = 'corr', names_to = 'data') %>%
        mutate(groups = 'unwanted') %>%
        mutate(groups = ifelse(var %in% bio.variables, 'wanted', groups)) %>%
        mutate(corr = ifelse(groups == 'wanted', 1 - corr, corr)) %>%
        group_by(data, groups) %>%
        summarise(corr = mean(corr)) %>%
        summarise(assess = corr[groups == 'wanted']/2 + corr[groups == 'unwanted']/2) %>%
        arrange(desc(assess)) %>%
        mutate(data = factor(data, levels = data))
    p.assess.ncgs <- ggplot(all, aes(x = data, y = assess)) +
        geom_bar(stat = 'identity', fill = 'grey') +
        ggtitle('Comparison of the NCG sets') +
        xlab('') +
        ylab('Correlations') +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            plot.title = element_text(size = 15),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
            axis.text.y = element_text(size = 12)
            )
    if(isTRUE(plot.output)) print(p.assess.ncgs)

    # Save the results
    printColoredMessage(
        message = '-- Save the results:',
        color = 'magenta',
        verbose = verbose
        )
    if(is.null(output.name)){
        output.name <- paste0('comparsion_', paste0(names(ncgs), collapse = '&'))
    }
    if(isTRUE(save.se.obj)){
        if (length(se.obj@metadata) == 0) {
            se.obj@metadata[['NCG']] <- list()
        }
        if(!ncg.type %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][[ncg.type]] <- list()
        }

        if(!ncg.group %in% names(se.obj@metadata[['NCG']][[ncg.type]])){
            se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]] <- list()
        }
        if(!'assessment.plot' %in% names(se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]])){
            se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]][['assessment.plot']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]][['assessment.plot']] )){
            se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]][['assessment.plot']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][[ncg.type]][[ncg.group]][['assessment.plot']][[output.name]]<- p.assess.ncgs
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(list(p.assess.ncgs = p.assess.ncgs))
    }
}


