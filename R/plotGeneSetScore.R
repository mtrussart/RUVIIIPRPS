#' Plot gene set score of a gene set.

#' @author Ramyar Molania

#' @description
#' This function computes all possible gene-gene pairwise ordinary and partial correlation of the assays in the
#' SummarizedExperiment object.

#' @details
#' Partial correlation is used to estimate correlation between two variables while controlling for third
#' variables

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol or a vector of symbols specifying the name(s) of the assay(s) in the SummarizedExperiment
#' object to calculate sample-wise scores. The default is set to "all," indicating all assays in the SummarizedExperiment
#' object will be selected.
#' @param reference.score Symbol.
#' @param gene.set.name Symbol. A symbol indicating the name to be used to save the score in the SummarizedExperiment object.
#' if 'NULL', the function will select a name as follow :
#' gene.set.name <- paste0('singscore|',length(c(upset.genes, downset.genes)),'genes')
#' @param plot.output Logical. Whether to plot the gene set score or not. The default is set to 'TRUE'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to 'TRUE'. Refer to the 'checkSeObj' function for more details.
#' @param save.se.obj Indicates whether to save the score results in the metadata of the SummarizedExperiment object
#' or to output the result as list. The default is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.
#' @export

plotGeneSetScore <- function(
        se.obj,
        assay.names = 'all',
        reference.score = NULL,
        gene.set.name,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGeneSetScore function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check inputs ####

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    # Scatter plot ####
    ref <- NULL
    if(!is.null(reference.score)){
        all.assay.score <- lapply(
            levels(assay.names),
            function(x){
                se.obj@metadata$Metrics[[x]]$global.level$GeneSetScore$singscore[[gene.set.name]]$score
            })
        colnames(all.assay.score) <- levels(assay.names)
        all.assay.score <- as.data.frame(all.assay.score)
        if(!reference.score %in% levels(assay.names)){
            all.assay.score$ref <- se.obj[[reference.score]]
        } else if (reference.score %in% levels(assay.names)){
            index <- colnames(all.assay.score) %in% reference.score
            colnames(all.assay.score)[index] <- 'ref'
        }
        scores <- NULL
        all.assay.score <- pivot_longer(data = all.assay.score, -ref, names_to = 'datasets', values_to = 'scores')
        all.assay.score$datasets <- factor(x = all.assay.score$datasets, levels = assay.names)
        p.all.scores.plot <- ggplot(data = all.assay.score, aes(x = scores, y = ref)) +
            geom_point(color = 'grey') +
            ggtitle('Comapre the gene set scores to the reference score.') +
            xlab('Scores') +
            ylab('Reference scores')  +
            facet_wrap(~datasets) +
            geom_smooth(formula = y ~ x, method = 'lm', colour = "darkgreen") +
            ggpubr::stat_cor(
                aes(label = after_stat(r.label)),
                color = "navy") +
            theme(panel.background = element_blank(),
                  axis.line = element_line(colour = 'black', linewidth = 1),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 9)
            )
        if (isTRUE(plot.output)) print(p.all.scores.plot)
    }
    if(is.null(reference.score)){
        all.assay.score <- lapply(
            levels(assay.names),
            function(x){
                se.obj@metadata$Metrics[[x]]$global.level$GeneSetScore$singscore[[gene.set.name]]$score
            })
        all.assay.score <- do.call(cbind, all.assay.score)
        colnames(all.assay.score) <- levels(assay.names)
        all.assay.score <- as.data.frame(all.assay.score)
        p.all.scores.plot <- ggpairs(data = all.assay.score) +
            theme(axis.line = element_line(colour = 'black', linewidth = 1),
                  axis.text.x = element_text(size = 6),
                  axis.text.y = element_text(size = 6))
        if (isTRUE(plot.output)) print(p.all.scores.plot)
    }

    # Save the results ####
    if(isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Save all the score plots into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        se.obj <- addOverallPlotToSeObj(
            se.obj = se.obj,
            slot = 'Plots',
            assessment.type = 'global.level',
            assessment = 'GeneSetSocore',
            method = 'singscore',
            variables = gene.set.name,
            file.name = 'general',
            plot.data = p.all.scores.plot
        )
        printColoredMessage(
            message = '------------The plotGeneSetScore function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    printColoredMessage(message = '------------The plotGeneSetScore function finished.',
                        color = 'white',
                        verbose = verbose)
    ## save the plots as list ####
    if(isFALSE(save.se.obj)){
        return(list(p.all.scores.plot = p.all.scores.plot))
    }
}
