#' Create PRPS for a continuous source of unwanted variation.

#' @author Ramyar Molania

#' @description
#' The function creates PRPS sets for a continuous source of unwanted variation defined in the 'main.uv.variable' argument.
#' For example to correct for library size if defined in the 'uv.variables' argument, several group of pseudo-samples
#' will be created by averaging the top and bottom-ranked samples by library size of the same biological subtype in each batch.
#' Then those pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set
#'

#' @details
#' Distinct group of PRPS are created for each source of unwanted variation defined in the 'main.uv.variable' argument
#' within an homogeneous group of samples. The grouping of samples are created based on each biological subtype defined
#' using the 'bio.variables' using the 'bio.clustering.method' selected and it might be also be combined with
#' 'other.uv.variables' if requested. By default 'other.uv.variables' is set up to NULL, meaning that the homogeneous
#' grouping of samples created is based solely on the biological subtype(s) defined in 'bio.variables'. If
#' 'other.uv.variables' is not set to NULL, the creation of homogeneous groups of samples is based on the biological subtype
#' defined in 'bio.variables' combined with 'other.uv.variables' using the 'other.uv.clustering.method' selected. A
#' pseudo-sample will be created for each homogeneous group of samples by averaging all the samples within that group of
#' samples. For example to correct for library size if defined in the 'main.uv.variable' argument, a pair of pseudo-samples
#' will be created for each homogeneous group of samples, each pair by averaging the top and the bottom-ranked samples by
#' library size.Each pair of pseudo-samples belonging to the same group will be defined as pseudo-replicates which
#' constitutes a PRPS set.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol specifying an assay (data) name within the SummarizedExperiment object. The selected
#' assay should be the one that will be used for RUV-III normalization.
#' @param bio.variables Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
#' subtypes, tumor purity, ... within the sample annotation (colData) of the SummarizedExperiment object. This can comprise
#' a vector containing of either categorical, continuous, or a combination of both variables.
#' @param main.uv.variable Symbol. A symbol representing the label of a continuous unwanted variable, such as library
#' size, ... within the sample annotation (colData) of the SummarizedExperiment object. Sets of PRPS will be created for
#' removing the impact of this variable from the data.
#' @param other.uv.variables Symbol. A symbol or symbols representing the label of unwanted variable(s), such as batch
#' effects, ... within the sample annotation (colData) of the SummarizedExperiment object. This can comprise a vector
#' containing either categorical, continuous, or a combination of both variables. These variable(s) will be considered
#' when generating PRPS data fo the "main.uv.variable". The default is set to 'NULL'. Refere to the details of the function
#' for more information.
#' @param min.sample.for.ps Numeric. A numeric value that indicates the minimum number of biologically homogeneous samples
#' to be averaged to create one pseudo-sample (PS). The default it is set to 3.
#' @param bio.clustering.method Symbol. A symbol indicating which clustering method should be used to group each continuous
#' biological variable, if they are specified. The options include 'kmeans', cut' and 'quantile'. The default is to 'kmeans'.
#' We refer to the createHomogeneousBioGroups() function for more details.
#' @param nb.bio.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous biological
#' variable. The default it set to 2. Then individual continuous biological variable will be divided into 2 groups using
#' the "bio.clustering.method" .
#' @param other.uv.clustering.method Symbol. A symbol specifying which clustering method should be used to group each
#' continuous unwanted variable specified in the "other.uv.variables". The options include 'kmeans', cut' and 'quantile'.
#' The default is to 'kmeans'. We refer to the createHomogeneousUVGroups() function for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous
#' unwanted variable specified in the "other.uv.variables". The default it set to 2.
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
#' @param plot.prps.map Logical. Indicates whether to generate the PRPS map plot for individual sources of unwanted variation.
#' The default is 'TRUE'.
#' @param output.name Symbol. A symbol to specify the name of all PRPS sets that will be created for all specified
#' source (s) of unwanted variation. The default is set to 'NULL'. Then, the function creates a name based on
#' paste0('prps_', uv.variable).
#' @param prps.group Symbol.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return SummarizedExperiment A SummarizedExperiment object containing the PRPS data or just PRPS data.

#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr group_by arrange slice desc add_count filter mutate
#' @importFrom tidyr %>%
#' @import ggplot2
#' @export

createPrPsForContinuousUV <- function(
        se.obj,
        assay.name,
        bio.variables,
        main.uv.variable,
        other.uv.variables,
        min.sample.for.ps = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        remove.na = 'both',
        plot.prps.map = TRUE,
        output.name = NULL,
        prps.group = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The createPrPsForContinuousUV function starts.',
                        color = 'white',
                        verbose = verbose)
    # Check the function inputs ####
    if (length(assay.name) > 1) {
        stop('The function can only take a single assay name.')
    } else if (length(main.uv.variable) > 1) {
        stop('The main.uv.variable should be a single continuous source of unwanted variation.')
    } else if (min.sample.for.ps <= 1) {
        stop('The minimum value for the min.sample.for.ps is 2.')
    } else if (var(se.obj[[main.uv.variable]]) == 0) {
        stop(paste(
                'The variable ',
                main.uv.variable,
                ' has no variation. Then, there is no need to create PRPS.')
        )
    } else if (!class(se.obj[[main.uv.variable]]) %in% c('numeric', 'integer')) {
        stop(paste(
            'The variable ',
            main.uv.variable,
            ' should be a continuous source of unwanted variation.'
        ))
    } else if (main.uv.variable %in% other.uv.variables){
        stop('The other.uv.variables should not contain the main.uv.variable.')
    } else if ( sum(c(main.uv.variable, other.uv.variables) %in% bio.variables) > 1){
        stop('The variable should be biological or unwantd variation.')
    } else if (pseudo.count < 0){
        stop('The value for pseudo.count can not be negative.')
    } else if (max(cat.cor.coef) > 1 | max(cont.cor.coef) > 1){
        stop('The maximum value for cat.cor.coef or cont.cor.coef cannot be more than 1.')
    }
    # check the SummarizedExperiment ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, other.uv.variables, bio.variables),
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # data transformation ####
    printColoredMessage(message = '-- Applying data transformation:',
                        color = 'magenta',
                        verbose = verbose
                        )
    if (isTRUE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count,' (pseudo.count',') on the ',
                assay.name,
                ' data before creating PRPS.'),
            color = 'blue',
            verbose = verbose)
        expre.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without log for creating PRPS.'),
            color = 'blue',
            verbose = verbose)
        expre.data <- assay(x = se.obj, i = assay.name)
    }
    # Create homogeneous biological groups of samples ####
    printColoredMessage(
        message = '-- Creating all possible homogeneous bioloical populations:',
        color = 'magenta',
        verbose = verbose
        )
    homo.bio.groups <- createHomogeneousBioGroups(
        se.obj = se.obj,
        bio.variables = bio.variables,
        nb.clusters = nb.bio.clusters,
        clustering.method = bio.clustering.method,
        assess.se.obj = FALSE,
        assess.variables = assess.variables,
        cat.cor.coef = cat.cor.coef,
        cont.cor.coef = cont.cor.coef,
        save.se.obj = FALSE,
        remove.na = 'none',
        verbose = verbose
        )
    # Create homogeneous groups of samples with respect to unwanted variation ####
    if (!is.null(other.uv.variables)) {
        printColoredMessage(
            message = paste0('-- Creating all possible homogeneous sample groups with respect to the',
                paste0(other.uv.variables, collapse = ' & '), ':'),
            color = 'magenta',
            verbose = verbose
            )
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.uv.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = assess.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        annot.data <- data.frame(
            uv.variable = se.obj[[main.uv.variable]],
            homo.uv.groups = homo.uv.groups,
            homo.bio.groups = homo.bio.groups
            )
        printColoredMessage(
            message = '-- Creating all possible homogeneous samples groups with respect to both bioloical and unwanted variables:',
            color = 'magenta',
            verbose = verbose
            )
        bio.batch <- NULL
        annot.data$bio.batch <- paste(
            annot.data$homo.uv.groups,
            annot.data$homo.bio.groups,
            sep = '||'
            )
        printColoredMessage(
            message = paste0(
                length(unique(annot.data$bio.batch)),
                ' ( maximun possible groups',
                length(unique(annot.data$homo.bio.groups)),
                ' biological groups * ',
                length(unique(annot.data$homo.uv.groups)),
                ' unwanted variation groups) ',
                'groups are created.'),
            color = 'blue',
            verbose = verbose
            )
        bio.batch.groups <- findRepeatingPatterns(
            vec = annot.data$bio.batch,
            n.repeat = c(2 * min.sample.for.ps)
            )
        if(length(bio.batch.groups) == 0){
            stop(paste0(
                    'There is not any sample groups with enough samples',
                    ' to create PRPS for the contenious source of unwanted variation.'))
        }
        cover <- c(length(bio.batch.groups)/length(unique(annot.data$bio.batch))) * 100
        printColoredMessage(
            message = paste0(length(bio.batch.groups),' (', round(cover, digits = 0),
                '%) groups have enough samples (2 * min.sample.for.ps) for creating PRPS.'),
            color = 'blue',
            verbose = verbose
            )
        # create PRPS data for each group of batch * bio ####
        printColoredMessage(
            message = '-- Creating a PRPS set with two pseudo samples for each individual homogeneous group:',
            color = 'magenta',
            verbose = verbose
            )
        annot.data <- annot.data %>%
            mutate(sOrder = c(1:ncol(se.obj))) %>%
            add_count(bio.batch) %>%
            filter(n >= c(min.sample.for.ps * 2))
        top <- annot.data %>%
            arrange(desc(uv.variable)) %>%
            group_by(bio.batch) %>%
            slice(1:min.sample.for.ps)
        bot <- annot.data %>%
            arrange(uv.variable) %>%
            group_by(bio.batch) %>%
            slice(1:min.sample.for.ps)
        selected.groups <- unique(bot$bio.batch)
        prps.sets <- lapply(
            selected.groups,
            function(x) {
                index.top <- top$sOrder[top$bio.batch == x]
                index.bot <- bot$sOrder[bot$bio.batch == x]
                ps.all <- cbind(rowMeans(expre.data[, index.top]),
                          rowMeans(expre.data[, index.bot]))
                colnames(ps.all) <- rep(paste0(main.uv.variable, '||', x), ncol(ps.all))
                ps.all
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0('In totall ', length(unique(colnames(prps.sets))), ' PRPS sets with ',
                ncol(prps.sets),' total number of pseudo-samples are created for the removal of the ',
                main.uv.variable,' effects.'),
            color = 'blue',
            verbose = verbose
            )
        # plot prps sets ####
        if(isTRUE(plot.prps.map)){
            prps.sets.plot <- sapply(
                selected.groups,
                function(x) {
                    index.top <- top$sOrder[top$bio.batch == x]
                    index.bot <- bot$sOrder[bot$bio.batch == x]
                    top.samples <- se.obj[[main.uv.variable]][index.top]
                    bot.samples <- se.obj[[main.uv.variable]][index.bot]
                    c(top.samples, bot.samples) }) %>%
                data.frame(.) %>%
                pivot_longer(everything(), values_to = 'uv.var', names_to = 'bio.batch') %>%
                arrange(bio.batch, uv.var) %>%
                mutate(groups = rep(rep(c('bottom', 'top'), each = min.sample.for.ps), length(selected.groups)))
            prps.map.plot <- data.frame(
                bio.batch = rep(main.uv.variable, ncol(se.obj)),
                uv.var = se.obj[[main.uv.variable]],
                groups = main.uv.variable) %>%
                rbind(., prps.sets.plot) %>%
                mutate(groups = factor(groups, levels= c(main.uv.variable, unique(prps.sets.plot$groups)))) %>%
                mutate(new.g = ifelse(groups == main.uv.variable, 'UV', 'PRPS sets')) %>%
                mutate(new.g = factor(new.g, levels = c('UV', 'PRPS sets'))) %>%
                ggplot(aes(x = bio.batch, y = uv.var, color = groups)) +
                geom_boxplot() +
                geom_point() +
                scale_color_manual(values = c('darkgreen', 'tomato', 'navy')) +
                facet_grid(.~new.g, scales = 'free_x', space = 'free') +
                scale_x_discrete(expand = c(0, 0.9)) +
                ylab(main.uv.variable) +
                xlab('Homogeneous groups') +
                ylim(c(min(se.obj[[main.uv.variable]]), max(se.obj[[main.uv.variable]]))) +
                geom_hline(yintercept = c(min(se.obj[[main.uv.variable]]), max(se.obj[[main.uv.variable]])), color = 'gray70') +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.x = element_text(size = 15))
            print(prps.map.plot)
        }
    } else if (is.null(other.uv.variables)){
        printColoredMessage(
            message = paste0(
                '-- Finding homogeneous biological groups with at least ',
                2 * min.sample.for.ps,
                ' (2* min.sample.for.ps) samples.'),
            color = 'magenta',
            verbose = verbose)
        bio.groups <- findRepeatingPatterns(
            vec = homo.bio.groups,
            n.repeat = c(2 * min.sample.for.ps)
            )
        if(length(bio.groups) == 0){
            stop(paste0(
                'There is not any homogeneous biological groups with enough samples',
                ' to create PRPS for the contenious source of unwanted variation.'))
        }
        cover <- c(length(bio.groups)/length(unique(homo.bio.groups))) * 100
        printColoredMessage(
            message = paste0(
                length(bio.groups),
                ' (',
                round(cover, digits = 0),
                '%) homogeneous biological groups have enough samples (2 * min.sample.for.ps) for creating PRPS.'),
            color = 'blue',
            verbose = verbose
            )
        # creating PRPS within individual homogeneous biological group ####
        printColoredMessage(
            message = '-- Creating a PRPS set with two pseudo samples for each individual homogeneous biological groups:',
            color = 'magenta',
            verbose = verbose
            )
        uv.variable <- n <- everything <- NULL
        annot.data <- data.frame(
            uv.variable = colData(se.obj)[[main.uv.variable]],
            bio.variable = homo.bio.groups)
        annot.data <- annot.data %>%
            mutate(sOrder = c(1:ncol(se.obj))) %>%
            add_count(bio.variable) %>%
            filter(n >= c(min.sample.for.ps * 2))
        top <- annot.data %>%
            arrange(desc(uv.variable)) %>%
            group_by(bio.variable) %>%
            slice(1:min.sample.for.ps)
        bot <- annot.data %>%
            arrange(uv.variable) %>%
            group_by(bio.variable) %>%
            slice(1:min.sample.for.ps)
        selected.groups <- unique(bot$bio.variable)
        prps.sets <- lapply(
            selected.groups,
            function(x) {
                index.top <- top$sOrder[top$bio.variable == x]
                index.bot <- bot$sOrder[bot$bio.variable == x]
                ps.all <- cbind(rowMeans(expre.data[, index.top]),
                          rowMeans(expre.data[, index.bot]))
                colnames(ps.all) <- rep(paste0(main.uv.variable, '||', x), ncol(ps.all))
                ps.all
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0(
                'In totall ',
                length(unique(colnames(prps.sets))),
                ' PRPS sets with ',
                ncol(prps.sets),
                ' total number of pseudo-samples are created for the removal of the ',
                main.uv.variable,
                ' effects.'),
            color = 'blue',
            verbose = verbose
        )
        # plot prps sets ####
        uv.var <- groups <- new.g <- bio.variable <- NULL
        if(isTRUE(plot.prps.map)){
            prps.sets.plot <- sapply(
                selected.groups,
                function(x) {
                    index.top <- top$sOrder[top$bio.variable == x]
                    index.bot <- bot$sOrder[bot$bio.variable == x]
                    top.samples <- se.obj[[main.uv.variable]][index.top]
                    bot.samples <- se.obj[[main.uv.variable]][index.bot]
                    c(top.samples, bot.samples)
                }) %>%
                data.frame(.) %>%
                pivot_longer(everything(), values_to = 'uv.var', names_to = 'bio.variable') %>%
                arrange(bio.variable, uv.var) %>%
                mutate(groups = rep(rep(c('bottom', 'top'), each = min.sample.for.ps), length(selected.groups)))
            prps.map.plot <- data.frame(
                bio.variable = rep(main.uv.variable, ncol(se.obj)),
                uv.var = se.obj[[main.uv.variable]],
                groups = main.uv.variable) %>%
                rbind(., prps.sets.plot) %>%
                mutate(groups = factor(groups, levels= c(main.uv.variable, unique(prps.sets.plot$groups)))) %>%
                mutate(new.g = ifelse(groups == main.uv.variable, 'UV', 'PRPS sets')) %>%
                mutate(new.g = factor(new.g, levels = c('UV', 'PRPS sets')))
            prps.map.plot <- ggplot(data = prps.map.plot, aes(x = uv.var , y = bio.variable , color = groups)) +
                geom_boxplot() +
                geom_point() +
                scale_color_manual(values = c('darkgreen', 'tomato', 'navy')) +
                facet_grid(new.g~., scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                xlab(main.uv.variable) +
                ylab('Homogeneous groups') +
                xlim(c(
                    min(se.obj[[main.uv.variable]]),
                    max(se.obj[[main.uv.variable]])
                )) +
                geom_hline(yintercept = c(
                    min(se.obj[[main.uv.variable]]),
                    max(se.obj[[main.uv.variable]])), color = 'gray70') +
                theme_bw() +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.y = element_text(size = 15)
                    )
            if(isTRUE(plot.prps.map)) print(prps.map.plot)
        }
    }
    # Saving the output ####
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
    } else if (is.null(other.uv.variables)) {
        output.name <- paste0(
            main.uv.variable,
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            assay.name)
    }
    ## select prps group name ####
    if (is.null(prps.group)) {
        prps.group <- paste0('prps|supervised|', main.uv.variable)
    }
    if (isTRUE(save.se.obj)) {
        ## check if PRPS exists
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
            ## check
            if (!'prps.map.plot' %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]])) {
                se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']])) {
                se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']][[output.name]] <- list()
            }
            se.obj@metadata[['PRPS']][['supervised']][[prps.group]][['prps.map.plot']][[output.name]] <- prps.map.plot
        }
        printColoredMessage(
            message = paste0('The PRPS data and PRPS map plots are saved to the metadata of the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The createPrPsForContinuousUV function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        printColoredMessage(
            message = paste0('The PRPS data and PRPS map plots are outputed as a list.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The createPrPsForContinuousUV function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(prps.sets = prps.sets, prps.map.plot = prps.map.plot))
    }
}
