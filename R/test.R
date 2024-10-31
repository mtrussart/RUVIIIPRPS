# raw <- log2(assay(brca.se.obj, 'RawCount') + 1)
# rle <- raw - rowMedians(raw)
# rle.med <- colMedians(rle)
# plot(rle.med)
# corr.med.ls <- Rfast::correls(y = rle.med, x = t(raw), type = 'spearman')
# corr.pur <- Rfast::correls(y = -brca.se.obj$tumour.purity.estimate, x = t(raw), type = 'spearman')
# corr.ls <- Rfast::correls(y = log2(brca.se.obj$library.size), x = t(raw), type = 'spearman')
# plot(corr.med.ls[,1], corr.pur[ , 1])
# plot(corr.med.ls[,1], corr.ls[ , 1])
#
# hist(corr.med.ls[,1], breaks = 100)
#
#
# cor.test(rle.med, brca.se.obj$tumour.purity.estimate)
# index <- corr.med.ls[,1] < .2
#
#
# raw <- lm(t(raw) ~ rle.med)
# raw <- t(raw$residuals)
#
# raw <- log2(assay(brca.se.obj, 'RawCount') + 1)
# raw <- raw[ index, ]
# rle <- raw - rowMedians(raw)
# rle.med <- colMedians(rle)
# plot(rle.med)
# corr.med.ls <- Rfast::correls(y = rle.med, x = t(raw), type = 'spearman')
# hist(corr.med.ls[,1], breaks = 100)
# plot(rle.med, brca.se.obj$tumour.purity.estimate)
