# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active rcavalcante/r3.4.4_bioc3.6:msig_v0.4.4
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')

library(bsseq)
library(methylSig)

pdata = read.table('samplesheet_errbs.txt', sep = '\t', header = T, as.is = T)
pdata$group = relevel(factor(pdata$group), ref = 'DS')
files = list.files(path = 'covs', full.names = TRUE)

bs = methylSigReadData(
    fileList = files,
    sample.ids = pdata$sample,
    header = FALSE,
    assembly = NA,
    context = 'CpG',
    treatment = c(1,0,1,0,1,0),
    destranded = FALSE,
    maxCount = 500, minCount = 10,
    filterSNPs = FALSE,
    num.cores = 1)

old_0.4.4_diff = methylSigCalc(
    meth = bs,
    groups = c(Treatment = 1, Control = 0),
    dispersion= 'both',
    local.disp = FALSE,
    local.meth = FALSE,
    min.per.group = c(2,2),
    T.approx = TRUE,
    num.cores = 1
)

old_0.4.4_gr = GRanges(
    seqnames = old_0.4.4_diff@data.chr,
    ranges = IRanges(
        start = old_0.4.4_diff@data.start,
        end = old_0.4.4_diff@data.start)
)

old_0.4.4_df = data.frame(
    meth_case = old_0.4.4_diff@results[, 'mu1'],
    meth_control = old_0.4.4_diff@results[, 'mu0'],
    meth_diff = old_0.4.4_diff@results[, 'meth.diff'],
    pvalue = old_0.4.4_diff@results[, 'pvalue'],
    fdr = old_0.4.4_diff@results[, 'qvalue'],
    disp_est = old_0.4.4_diff@results[, 'theta'],
    log_lik_ratio = old_0.4.4_diff@results[, 'logLikRatio'],
    df = old_0.4.4_diff@results[, 'df'],
    stringsAsFactors = F
)

mcols(old_0.4.4_gr) = old_0.4.4_df

old_0.4.4_local_diff = methylSigCalc(
    meth = bs,
    groups = c(Treatment = 1, Control = 0),
    dispersion= 'both',
    local.disp = TRUE,
    winsize.disp = 200,
    local.meth = FALSE,
    winsize.meth = 200,
    min.per.group = c(2,2),
    T.approx = TRUE,
    num.cores = 1
)

old_0.4.4_local_gr = GRanges(
    seqnames = old_0.4.4_local_diff@data.chr,
    ranges = IRanges(
        start = old_0.4.4_local_diff@data.start,
        end = old_0.4.4_local_diff@data.start)
)

old_0.4.4_local_df = data.frame(
    meth_case = old_0.4.4_local_diff@results[, 'mu1'],
    meth_control = old_0.4.4_local_diff@results[, 'mu0'],
    meth_diff = old_0.4.4_local_diff@results[, 'meth.diff'],
    pvalue = old_0.4.4_local_diff@results[, 'pvalue'],
    fdr = old_0.4.4_local_diff@results[, 'qvalue'],
    disp_est = old_0.4.4_local_diff@results[, 'theta'],
    log_lik_ratio = old_0.4.4_local_diff@results[, 'logLikRatio'],
    df = old_0.4.4_local_diff@results[, 'df'],
    stringsAsFactors = F
)

mcols(old_0.4.4_local_gr) = old_0.4.4_local_df

save(list = c('old_0.4.4_gr', 'old_0.4.4_local_gr'), file = 'rda/old_0.4.4_msig.rda')
