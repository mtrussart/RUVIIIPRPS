#' Assess the W matrix of RUV-III

#' @author Ramyar Molania

#' @description
#' This functions assesses the correlation between the W of RUV-III normalized data with known variables. Tt calculates
#' the association of the columns of the W matrix with both specified biological and unwanted variables. The correlation
#' values obtained for the biological variables will be subtracted from 1. Finally, the average of all these correlations
#' is computed to determine the final correlation values. A higher correlation indicates better performance.

#' @param se.obj A summarized experiment object.
#' @param variables Symbol. A symbol or symbols representing the label of variable(s), such as cancer subtypes, tumour
#' purity, librar size ... within the SummarizedExperiment object. This can comprise a vector containing either categorical,
#' continuous, or a combination of both variables.
#' @param bio.variables Symbol. A symbol or symbols representing the label of biological variable(s), such as cancer
#' subtypes, tumour purity, ... within the SummarizedExperiment object. This can comprise a vector containing either
#' categorical, continuous, or a combination of both variables.
#' @param uv.variables Symbol. A symbol or symbols representing the label of unwanted variable(s), such as cancer
#' batch effects, library size, ... within the SummarizedExperiment object. This can comprise a vector containing either
#' categorical, continuous, or a combination of both variables.
#' @param compare.w Logical.  Specifies whether to compare different W matrices in their ability to capture unwanted
#' variation and lack correlation with biological variables. See the function's details for more information.
#' @param plot.output Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.
#' @param save.se.obj Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.
#' @param verbose Logical. Indicates whether to plot the output or not. The default is set to 'TRUE'.

#' @importFrom dplyr bind_rows summarise mutate group_by
#' @importFrom fastDummies dummy_cols
#' @importFrom gtools mixedorder
#' @importFrom tidyr pivot_longer
#' @export

assessW <- function(
        se.obj ,
        variables ,
        bio.variables = NULL,
        uv.variables = NULL,
        compare.w = FALSE,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){

    # Check inputs ####
    if(isTRUE(compare.w)){
        if(is.null(bio.variables)){
            stop('To compare different k values, the "bio.variables" must be specified.')
        }
        if(is.null(uv.variables)){
            stop('To compare different k values, the "uv.variables" must be specified.')
        }
    }
    if(isFALSE(compare.w)){
        if(is.null(variables)){
            stop('The "variables" cannot be "NULL".')
        }
    }

    #  Compare different W values ####
    if(isTRUE(compare.w)){
        printColoredMessage(
            message = '-- Compare different W values',
            color = 'magenta',
            verbose = verbose
            )
        data.names <- names(se.obj@metadata$RUVIII$W)
        names.order <- mixedorder(data.names)
        data.names <- data.names[names.order]
        all.w <- lapply(
            data.names,
            function(x) se.obj@metadata$RUVIII[['W']][[x]]
        )
        data.names <- unlist(lapply(
            data.names,
            function(x){
                split.name <- strsplit(x, split = '_')
                char.len <- length(split.name[[1]])
                paste(split.name[[1]][char.len-1], split.name[[1]][char.len], sep = ':')
            }))
        names(all.w) <- data.names

        ## variables class
        all.vars <- c(bio.variables, uv.variables)
        class.all.vars <- sapply(
            all.vars,
            function(x) class(colData(x = se.obj)[[x]])
            )
        cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
        cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])

        ## linear regression class
        cont.vars.r.squareds <- lapply(
            names(all.w),
               function(x) {
                   w <- all.w[[x]]
                   r.squareds <- sapply(
                       cont.vars,
                       function(y) {
                           lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ w))
                           lm.reg$r.squared
                       })
                   names(r.squareds) <- cont.vars
                   r.squareds
               })
        names(cont.vars.r.squareds) <- names(all.w)
        cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)

        ## vector correlation
        cat.vars.vec.corr <- lapply(
            names(all.w),
            function(x) {
                w <- all.w[[x]]
                vec.corr <- sapply(
                    cat.vars,
                    function(y) {
                        catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                        catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                        cca <- cancor(x = w, y = catvar.dummies)
                        1 - prod(1 - cca$cor ^ 2)
                        })
                names(vec.corr) <- cat.vars
                vec.corr
            })
        names(cat.vars.vec.corr) <- names(all.w)
        cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)

        ### put all together
        groups <- corr <- data <- assess <- variable <- NULL
        all <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
            round(digits = 3)
        all.a <- mutate(.data = all, var = row.names(all)) %>%
            pivot_longer(cols = -var, values_to = 'corr', names_to = 'data') %>%
            mutate(data = factor(data, levels = gsub(':', '.', names(all.w))) )
        p.w.1 <- ggplot(data = all.a, aes(x = data, y = corr, color = var)) +
            geom_point(size = 3) +
            xlab('W') +
            ylab('Correlations') +
            ylim(c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12)
                ) +
            guides(color = guide_legend(title = "Variables"))
        all <- mutate(.data = all, var = row.names(all)) %>%
            pivot_longer(cols = -var, values_to = 'corr', names_to = 'data') %>%
            mutate(groups = 'unwanted') %>%
            mutate(groups = ifelse(var %in% bio.variables, 'wanted', groups)) %>%
            mutate(corr = ifelse(groups == 'wanted', 1 - corr, corr)) %>%
            group_by(data, groups) %>%
            summarise(corr = mean(corr)) %>%
            summarise(assess = corr[groups == 'wanted']/2 + corr[groups == 'unwanted']/2) %>%
            mutate(data = factor(data, levels = gsub(':', '.', names(all.w))) )
        p.w.2 <- ggplot(all, aes(x = data, y = assess)) +
            geom_bar(stat = 'identity', fill = 'grey') +
            xlab('W') +
            ylab('Correlations') +
            ylim(c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12)
                )
        p.w <- ggarrange(p.w.1, p.w.2)
        if(isTRUE(plot.output)) print(p.w)
    }

    #  Assess W ####
    if(isFALSE(compare.w)){
        printColoredMessage(
            message = '-- Assess different W values',
            color = 'magenta',
            verbose = verbose)
        data.names <- gsub('\\.', '_', names(se.obj@metadata$RUVIII$W))
        data.names <- mixedorder(data.names)
        all.w <- lapply(
            data.names,
            function(x) se.obj@metadata$RUVIII[['W']][[x]]
        )
        names(all.w) <- gsub('\\.', '_', names(se.obj@metadata$RUVIII$W))

        ## variables class
        class.all.vars <- sapply(
            variables,
            function(x) class(colData(x = se.obj)[[x]])
        )
        cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
        cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])

        ## linear regression class
        cont.vars.r.squareds <- lapply(
            names(all.w),
            function(x) {
                w <- all.w[[x]]
                r.squareds <- sapply(
                    cont.vars,
                    function(y) {
                        lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ w))
                        lm.reg$r.squared
                    })
                names(r.squareds) <- cont.vars
                r.squareds
            })
        names(cont.vars.r.squareds) <- names(all.w)
        cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)

        ## vector correlation
        cat.vars.vec.corr <- lapply(
            names(all.w),
            function(x) {
                w <- all.w[[x]]
                vec.corr <- sapply(
                    cat.vars,
                    function(y) {
                        catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                        catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                        cca <- cancor(x = w, y = catvar.dummies)
                        1 - prod(1 - cca$cor ^ 2)
                    })
                names(vec.corr) <- cat.vars
                vec.corr
            })
        names(cat.vars.vec.corr) <- names(all.w)
        cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)

        ### put all together
        all.corrs <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
            round(digits = 3)
        all.corrs <- mutate(variable = row.names(all.corrs))
        all.corrs <- pivot_longer(
            data = all.corrs,
            cols = -variable,
            values_to = 'corr',
            names_to = 'data') %>%
            mutate(data = factor(data, levels = gsub('\\.', '_', names(se.obj@metadata$RUVIII$W)) ))
        all.corrs <- data.frame(all.corrs)
        p.w <- ggplot(all.corrs, aes(x = data, y = corr, group = variable)) +
            geom_line(aes(color = variable), size = 1) +
            geom_point(aes(color = variable), size = 3) +
            ylab('Correlations') +
            xlab('W') +
            theme(
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12)
                )
        if(isTRUE(plot.output)) print(p.w)
    }
    # Save the results
    printColoredMessage(
        message = '-- Save the results:',
        color = 'magenta',
        verbose = verbose)
    if(isTRUE(save.se.obj)){
        if (length(se.obj@metadata) == 0) {
            se.obj@metadata[['RUVIII']] <- list()
        }
        # check if RUVIII already exists in the metadata
        if (!'RUVIII' %in% names(se.obj@metadata)) {
            se.obj@metadata[['RUVIII']] <- list()
        }
        ## check if W already exists in the RUVIII
        if (!'CompareW' %in% names(se.obj@metadata[['RUVIII']])) {
            se.obj@metadata[['RUVIII']][['CompareW']] <- list()
        }
        se.obj@metadata[['RUVIII']][['CompareW']] <- p.w
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(p.w = p.w)
    }
}


m <- matrix(rnorm(n = 3000), 50, 60)
colnames(m) <- paste0('sample1', 1:60)
colnames(m)[1:2] <- paste0('sampleD', 'sampleD')
colnames(m)[10:11] <- paste0('sampleA', 'sampleA')
colnames(m)[40:41] <- paste0('sampleC', 'sampleC')
colnames(m)[20:21] <- paste0('sampleS', 'sampleS')
colnames(m)[18:19] <- paste0('sampleP', 'sampleP')
colnames(m)[30:31] <- paste0('sampleM', 'sampleM')

M <- replicate.matrix(colnames(m))


ff <- function (Y, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE,
          average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE)
{
    if (is.data.frame(Y))
        Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    M = replicate.matrix(M)
    ctl = tological(ctl, n)
    if (inputcheck) {
        if (m > n)
            warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
        if (sum(is.na(Y)) > 0)
            warning("Y contains missing values.  This is not supported.")
        if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
            0)
            warning("Y contains infinities.  This is not supported.")
    }
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (ncol(M) >= m)
        newY = Y
    else if (is.null(k)) {
        ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
        newY = (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*%
                                                              ycyctinv)) %*% Y
        fullalpha = NULL
    }
    else if (k == 0) {
        newY = Y
        fullalpha = NULL
    }
    else {
        if (is.null(fullalpha)) {
            Y0 = residop(Y, M)
            fullalpha = t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M),
                                                      sum(ctl)), drop = FALSE]) %*% Y
        }
        alpha = fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
        ac = alpha[, ctl, drop = FALSE]
        W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY = Y - W %*% alpha
    }
    if (average)
        newY = ((1/apply(M, 2, sum)) * t(M)) %*% newY
    if (!return.info)
        return(newY)
    else return(list(newY = newY, M = M, W = W, fullalpha = fullalpha))
}
library(ruv)

tological <- function (ctl, n)
{
    ctl2 = rep(FALSE, n)
    ctl2[ctl] = TRUE
    return(ctl2)
}



