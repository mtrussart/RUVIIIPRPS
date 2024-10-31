#
#
#
#
#
#
#
# getMaximumK <- function(
#         se.obj,
#         assay.name,
#         prps.type = 'supervised',
#         prps.group,
#         prps.set.names = NULL,
#         ncg.type = 'supervised',
#         ncg.group,
#         ncg.set.names = NULL,
#         technical.replicates = NULL,
#         k = 1,
#         apply.log = TRUE,
#         data.to.log = 'assay',
#         pseudo.count = 1,
#         eta = NULL,
#         include.intercept = TRUE,
#         return.info = FALSE,
#         assess.se.obj = TRUE,
#         remove.na = 'none',
#         output.name = NULL,
#         save.se.obj = TRUE,
#         verbose = TRUE
# ) {
#     printColoredMessage(message = '------------The RUVIII.PRPS function starts:',
#                         color = 'white',
#                         verbose = verbose)
#     # Check inputs ####
#     if(!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name) | assay.name == 'all'){
#         stop('The "assay.name" must be the name of an assay in the SummarizedExperiment object.')
#     }
#     if(!assay.name %in% names(assays(se.obj))){
#         stop('The assay name cannot be found in the SummarizedExperiment object.')
#     }
#     if(!prps.type %in% c('supervised', 'un.supervised', 'both', 'none')){
#         stop('The "prps.type" must be of the "supervised", "un.supervised", "both" or "none".')
#     }
#     if(prps.type != 'none'){
#         if(length(se.obj@metadata) == 0 | !'PRPS' %in% names(se.obj@metadata)){
#             stop('Any PRPS data cannot be found in the SummarizedExperiment object.')
#         }
#         if(is.null(prps.type) | is.logical(prps.type)){
#             stop('The "prps.type" cannot be NULL or logical.')
#         }
#         if(prps.type == 'supervised'){
#             if(!'supervised' %in% names(se.obj@metadata$PRPS)){
#                 stop('The "supervised" object does not exist in the PRPS slot of the metada of the SummarizedExperiment object.')
#             }
#             if(length(se.obj@metadata$PRPS$supervised) == 0){
#                 stop('The "supervised" object in the PRPS slot of the metada of the SummarizedExperiment object does not have any data.')
#             }
#             if(sum(prps.group %in% names(se.obj@metadata$PRPS$supervised)) != length(prps.group)){
#                 stop('The "prps.group" cannot be found in the "supervised" part of th PRPS slot in the the metada.')
#             }
#             if(!is.null(prps.set.names)){
#                 if(sum(prps.set.names %in% names(se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data)) != length(prps.set.names)){
#                     stop('All or some of the "prps.set.names" cannot be found in the "supervised" part of th PRPS slot in the the metada.')
#                 }
#             }
#         }
#         if(prps.type == 'un.supervised'){
#             if(!'un.supervised' %in% names(se.obj@metadata$PRPS)){
#                 stop('The "un.supervised" object does not exist in the PRPS slot of the metada of the SummarizedExperiment object.')
#             }
#             if(length(se.obj@metadata$PRPS$un.supervised) == 0){
#                 stop('The "un.supervised" object in the PRPS slot of the metada of the SummarizedExperiment object does not have any data.')
#             }
#             if(sum(prps.group %in% names(se.obj@metadata$PRPS$un.supervised)) != length(prps.group) ){
#                 stop('The "prps.group" cannot be found in the "un.supervised" part of th PRPS slot in the the metada.')
#             }
#             if(!is.null(prps.set.names)){
#                 if(sum(prps.set.names %in% names(se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data) ) != length(prps.set.names) ){
#                     stop('All or some of the "prps.set.names" cannot be found in the "un.supervised" part of th PRPS slot in the the metada.')
#                 }
#             }
#         }
#         if (prps.type == 'both'){
#             if(!'supervised' %in% names(se.obj@metadata$PRPS)){
#                 stop('The "supervised" object does not exist in the PRPS slot of the metada of the SummarizedExperiment object.')
#             }
#             if(length(se.obj@metadata$PRPS$supervised) == 0){
#                 stop('The "supervised" object in the PRPS slot of the metada of the SummarizedExperiment object does not have any data.')
#             }
#             if(!'un.supervised' %in% names(se.obj@metadata$PRPS)){
#                 stop('The "un.supervised" object does not exist in the PRPS slot of the metada of the SummarizedExperiment object.')
#             }
#             if(length(se.obj@metadata$PRPS$un.supervised) == 0){
#                 stop('The "un.supervised" object in the PRPS slot of the metada of the SummarizedExperiment object does not have any data.')
#             }
#             all.prps.group <- c(names(se.obj@metadata$PRPS$supervised), names(se.obj@metadata$PRPS$un.supervised))
#             if(sum(prps.group %in% all.prps.group) != length(prps.group)){
#                 stop('The "prps.group" cannot be found in either "supervised" or "un.supervised" parts of th PRPS slot in the the metada.')
#             }
#             if(!is.null(prps.set.names)){
#                 supervised.prps.group <- all.prps.group[all.prps.group %in% names(se.obj@metadata$PRPS$supervised)]
#                 un.supervised.prps.group <- all.prps.group[all.prps.group %in% names(se.obj@metadata$PRPS$un.supervised)]
#                 all.prps.set.names <- c(
#                     names(se.obj@metadata$PRPS$supervised[[supervised.prps.group]]$prps.data),
#                     names(se.obj@metadata$PRPS$un.supervised[[un.supervised.prps.group]]$prps.data)
#                 )
#                 if(sum(prps.set.names %in% all.prps.set.names) != length(prps.set.names)){
#                     stop('All or some of the "prps.set.names" cannot be found in the "supervised" or "un.supervised" parts of th PRPS slot in the the metada.')
#                 }
#             }
#         }
#     }
#     if(!ncg.type %in% c('supervised', 'un.supervised', 'pre.selected', 'all.genes')){
#         stop('The "ncg.type" must be of the "supervised", "un.supervised", "pre.selected" or "all.genes".')
#     }
#     if(ncg.type != 'all.genes'){
#         if(!'NCG' %in% names(se.obj@metadata)){
#             stop('The "NCG" object does not exist in the metada of the SummarizedExperiment object.')
#         }
#         if(ncg.type == 'supervised'){
#             if(!'supervised' %in% names(se.obj@metadata$NCG)){
#                 stop('The supervised "NCG" object does not exist in the NCG slot of the metada of the SummarizedExperiment object.')
#             }
#             if(!ncg.group %in% names(se.obj@metadata$NCG$supervised)) {
#                 stop('The ncg.group object does not exist in the NCG slot of the metada of the SummarizedExperiment object.')
#             }
#             if(!is.null(ncg.set.names)){
#                 if(sum(ncg.set.names %in% names(se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set)) != length(ncg.set.names)){
#                     stop('All or some of "ncg.set.names" cannot be found in the SummarizedExperiment object.')
#                 }
#             }
#         }
#         if(ncg.type == 'un.supervised'){
#             if(!'un.supervised' %in% names(se.obj@metadata$NCG)){
#                 stop('The un.supervised "NCG" object does not exist in the NCG slot of the metada of the SummarizedExperiment object.')
#             }
#             if(!ncg.group %in% names(se.obj@metadata$NCG$un.supervised)){
#                 stop('The ncg.group object does not exist in the NCG slot of the metada of the SummarizedExperiment object.')
#             }
#             if(!is.null(ncg.set.names)){
#                 if(sum(ncg.set.names %in% names(se.obj@metadata$NCG$un.supervised[[ncg.group]])) != length(ncg.set.names)){
#                     stop('All or some of "ncg.set.names" cannot be found in the SummarizedExperiment object.')
#                 }
#             }
#         }
#         if(ncg.type == 'pre.selected'){
#             if(!'pre.selected' %in% names(se.obj@metadata$NCG)){
#                 stop('"pre.selected" NCG data canot be found in the metadata of the SummarizedExperiment object.')
#             }
#             if(sum(ncg.group %in% names(se.obj@metadata$NCG$pre.selected)) != length(ncg.group)){
#                 stop('All or some of "ncg.group" cannot be found in the SummarizedExperiment object.')
#             }
#         }
#     }
#
#     if(prps.type == 'none' & is.null(technical.replicates)){
#         stop('To run RUVIII either the "prps.type" or "technical.replicates" or both must be specified.')
#     }
#     if(!is.null(technical.replicates)){
#         if(length(technical.replicates) != ncol(se.obj)){
#             stop('The length of the "technical.replicate" must be the same as the length of the column name of the the SummarizedExperiment object.')
#         }
#     }
#     if (min(k) <= 0){
#         stop('k cannot be 0 or negative values.')
#     }
#     if(!data.to.log %in% c("assay", "prps", "both")){
#         stop('The "data.to.log" must be one of the "assay", "prps" or "both".')
#     }
#
#     # Assess the SummarizedExperiment object ####
#     if (isTRUE(assess.se.obj)) {
#         se.obj <- checkSeObj(
#             se.obj = se.obj,
#             assay.names = assay.name,
#             variables = NULL,
#             remove.na = remove.na,
#             verbose = verbose)
#     }
#
#     # Obtain control assays ####
#     ## technical replicates ####
#     if(!is.null(technical.replicates)){
#         printColoredMessage(
#             message = '-- Obtain technical replicates data:',
#             color = 'magenta',
#             verbose = verbose
#         )
#         rep.samples <- findRepeatingPatterns(
#             vec = technical.replicates,
#             n.repeat = 2
#         )
#         if(length(rep.samples) == 0){
#             stop('All the sample names of the "technical.replicates" are unique. Individual technical replicate sets must have the same names.')
#         } else {
#             printColoredMessage(
#                 message = paste0(
#                     'There is(are) ', length(unique(rep.samples)), ' sets of technical replicates in the data.'),
#                 color = 'blue',
#                 verbose = verbose
#             )
#         }
#         colnames(se.obj) <- technical.replicates
#     }
#     ## prps ####
#     if(prps.type != 'none'){
#         ### obtain PRPS data ####
#         printColoredMessage(
#             message = '-- Obtain the PRPS data in the SummarizedExperiment object:',
#             color = 'magenta',
#             verbose = verbose
#         )
#         ### supervised ####
#         if(prps.type == 'supervised'){
#             printColoredMessage(
#                 message = paste0('- The prps group "', prps.group, '" has ',
#                                  length(se.obj@metadata$PRPS$supervised[[prps.group]]$prps.data),
#                                  ' PRPS sub_group(s).'),
#                 color = 'blue',
#                 verbose = verbose
#             )
#             if(!is.null(prps.set.names)){
#                 prps.data <- lapply(
#                     prps.set.names,
#                     function(x) {
#                         se.obj@metadata$PRPS$supervised[[prps.group]]$prps.data[[x]]
#                     })
#                 names(prps.data) <- prps.set.names
#             }
#             if(is.null(prps.set.names)){
#                 prps.set.names <- names(se.obj@metadata$PRPS$supervised[[prps.group]]$prps.data)
#                 prps.data <- lapply(
#                     prps.set.names,
#                     function(x) {
#                         se.obj@metadata$PRPS$supervised[[prps.group]]$prps.data[[x]]
#                     })
#                 names(prps.data) <- prps.set.names
#             }
#         }
#         ### un supervised ####
#         if (prps.type == 'un.supervised'){
#             printColoredMessage(
#                 message = paste0('- The prps group "', prps.group, '" has ',
#                                  length(se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data),
#                                  ' PRPS sub_group(s).'),
#                 color = 'blue',
#                 verbose = verbose
#             )
#             if(!is.null(prps.set.names)){
#                 prps.data <- lapply(
#                     prps.set.names,
#                     function(x) {
#                         se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data[[x]]
#                     })
#                 names(prps.data) <- prps.set.names
#             }
#             if(is.null(prps.set.names)){
#                 prps.set.names <- names(se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data)
#                 prps.data <- lapply(
#                     prps.set.names,
#                     function(x) {
#                         se.obj@metadata$PRPS$un.supervised[[prps.group]]$prps.data[[x]]
#                     })
#                 names(prps.data) <- prps.set.names
#             }
#         }
#         ### both ####
#         if (prps.type == 'both'){
#             if(is.null(prps.group)){
#                 prps.data <- c(
#                     se.obj@metadata$PRPS$supervised,
#                     se.obj@metadata$PRPS$un.supervised$prps.data
#                 )
#             } else{
#                 prps.data.sup <- lapply(
#                     prps.group,
#                     function(x) se.obj@metadata$PRPS$supervised[[x]])
#                 names(prps.data.sup) <- prps.group
#                 prps.data.unsup <- lapply(
#                     prps.group,
#                     function(x) se.obj@metadata$PRPS$un.supervised[[x]])
#                 names(prps.data.unsup) <- prps.group
#                 prps.data <- c(prps.data.sup, prps.data.unsup)
#                 prps.data <- Filter(Negate(is.null), prps.data)
#             }
#             printColoredMessage(
#                 message = paste0(length(prps.data), ' supervised/unsupervised PRPS set(s) are found in the SummarizedExperiment object.'),
#                 color = 'blue',
#                 verbose = verbose
#             )
#         }
#         ### check dimension and gene names in the prps data ####
#         printColoredMessage(
#             message = '- Check row names and number of pseudo samples in the PRPS data:',
#             color = 'orange',
#             verbose = verbose
#         )
#         check.prps.data <- lapply(
#             1:length(prps.data),
#             function(x){
#                 if (nrow(se.obj) != nrow(prps.data[[x]]) ) {
#                     stop('The number of genes between the SummarizedExperiment object and the prps data is not the same.')
#                 }
#                 if (isFALSE(all.equal(row.names(se.obj), row.names(prps.data[[x]]))) ) {
#                     stop('The order of genes are not the same between the SummarizedExperiment object and the PRPS data.')
#                 }
#                 rep.samples <- findRepeatingPatterns(
#                     vec = colnames(prps.data[[x]]),
#                     n.repeat = 2
#                 )
#                 if (length(rep.samples) == 0) {
#                     stop('The names of the columns in the "prps.data" are all unique.')
#                 } else if (isTRUE(length(rep.samples) != length(unique(colnames(prps.data[[x]]))))) {
#                     stop('Some names of the columns in the "prps.data" are all unique.')
#                 }
#             })
#         printColoredMessage(
#             message = '- The PRPS data has no issues.',
#             color = 'blue',
#             verbose = verbose
#         )
#         ### merge all the PRPS data ####
#         printColoredMessage(
#             message = '- Merge all the PPRS data:',
#             color = 'orange',
#             verbose = verbose
#         )
#         prps.data <- do.call(cbind, prps.data)
#         printColoredMessage(
#             message = '- The final PRPS data has:',
#             color = 'blue',
#             verbose = verbose
#         )
#         printColoredMessage(
#             message = paste0('* ', ncol(prps.data), ' pseudo-samples'),
#             color = 'blue',
#             verbose = verbose
#         )
#         printColoredMessage(
#             message = paste0('* ', length(findRepeatingPatterns(vec = colnames(prps.data), n.repeat = 2)), ' PRPS sets.'),
#             color = 'blue',
#             verbose = verbose
#         )
#     } else prps.data <- NULL
#
#     # Control genes ####
#     if(isTRUE(ncg.type != 'all.genes')){
#         ## obtain NCG ####
#         printColoredMessage(
#             message = '-- Obtain the selected NCGs in the SummarizedExperiment object',
#             color = 'magenta',
#             verbose = verbose
#         )
#         ## supervised ####
#         if(ncg.type == 'supervised'){
#             if(is.null(ncg.set.names)){
#                 printColoredMessage(
#                     message = paste0('- The supervised NCG group "', ncg.group, '" has ',
#                                      length(se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set),
#                                      ' NCG sub_group(s).'),
#                     color = 'blue',
#                     verbose = verbose
#                 )
#                 if(length(se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set) > 1){
#                     printColoredMessage(
#                         message = '- The union of the NCGs sets will be used.',
#                         color = 'blue',
#                         verbose = verbose
#                     )
#                 }
#                 ncg <- se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set
#             }
#             if(!is.null(ncg.set.names)){
#                 ncg <- lapply(
#                     ncg.set.names,
#                     function(x) se.obj@metadata$NCG$supervised[[ncg.group]]$ncg.set[[x]]
#                 )
#                 if(length(ncg.set.names) > 1){
#                     printColoredMessage(
#                         message = '- The union of the NCGs sets will be used.',
#                         color = 'blue',
#                         verbose = verbose
#                     )
#                 }
#                 names(ncg) <- ncg.set.names
#             }
#         }
#
#         ## un supervised ####
#         if (ncg.type == 'un.supervised'){
#             if(is.null(ncg.set.names)){
#                 printColoredMessage(
#                     message = paste0('- The un.supervised NCG group "', prps.group, '" has ',
#                                      length(se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set),
#                                      ' PRPS sub_group(s).'),
#                     color = 'blue',
#                     verbose = verbose
#                 )
#                 if(length(se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set) > 1){
#                     printColoredMessage(
#                         message = '- The union of the NCGs sets will be used.',
#                         color = 'blue',
#                         verbose = verbose
#                     )
#                 }
#                 ncg <- se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set
#             }
#             if(!is.null(ncg.set.names)){
#                 ncg <- lapply(
#                     ncg.set.names,
#                     function(x) se.obj@metadata$NCG$un.supervised[[ncg.group]]$ncg.set[[x]])
#                 if(length(ncg.set.names) > 1){
#                     printColoredMessage(
#                         message = '- The union of the NCGs sets will be used.',
#                         color = 'blue',
#                         verbose = verbose
#                     )
#                 }
#                 names(ncg) <- ncg.set.names
#             }
#         }
#         ## pre selected ####
#         if (ncg.type == 'pre.selected'){
#             if(is.null(ncg.group)){
#                 ncg <- se.obj@metadata$NCG$pre.selected
#             } else{
#                 ncg <- lapply(
#                     ncg.group,
#                     function(x) se.obj@metadata$NCG$pre.selected[[x]])
#                 names(ncg) <- ncg.group
#             }
#             printColoredMessage(
#                 message = paste0(length(ncg), ' supervised NCG set(s) are found in the SummarizedExperiment object.'),
#                 color = 'blue',
#                 verbose = verbose)
#         }
#
#     } else if (ncg.type == 'all.genes') {
#         ncg <- list(all.genes = row.names(se.obj) %in% row.names(se.obj))
#     }
#     ncg <- do.call(cbind, ncg)
#     ncg <- rowSums(ncg) != 0
#     printColoredMessage(
#         message = paste0('- The final set of NCGs contains ', sum(ncg),  ' genes.'),
#         color = 'blue',
#         verbose = verbose
#     )
#
#     # Check max k for each normalization ####
#     printColoredMessage(
#         message = '-- Check max k for each normalization:',
#         color = 'magenta',
#         verbose = verbose
#     )
#     control.assay <- findRepeatingPatterns(
#         vec = c(colnames(prps.data), colnames(se.obj)),
#         n.repeat = 2
#     )
#     max.k.values <- min(length(control.assay), sum(ncg))
#     printColoredMessage(
#         message = paste0('- The maximum k for the curret control assays and NCGs is ', max.k.values, '.'),
#         color = 'blue',
#         verbose = verbose
#     )
# }
#
