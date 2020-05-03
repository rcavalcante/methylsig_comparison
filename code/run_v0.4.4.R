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

v0.4.4_diff = methylSigCalc(
    meth = bs,
    groups = c(Treatment = 1, Control = 0),
    dispersion= 'both',
    local.disp = FALSE,
    local.meth = FALSE,
    min.per.group = c(2,2),
    T.approx = TRUE,
    num.cores = 1
)

v0.4.4_gr = GRanges(
    seqnames = v0.4.4_diff@data.chr,
    ranges = IRanges(
        start = v0.4.4_diff@data.start,
        end = v0.4.4_diff@data.start)
)

v0.4.4_df = data.frame(
    meth_case = v0.4.4_diff@results[, 'mu1'],
    meth_control = v0.4.4_diff@results[, 'mu0'],
    meth_diff = v0.4.4_diff@results[, 'meth.diff'],
    pvalue = v0.4.4_diff@results[, 'pvalue'],
    fdr = v0.4.4_diff@results[, 'qvalue'],
    disp_est = v0.4.4_diff@results[, 'theta'],
    log_lik_ratio = v0.4.4_diff@results[, 'logLikRatio'],
    df = v0.4.4_diff@results[, 'df'],
    stringsAsFactors = F
)

mcols(v0.4.4_gr) = v0.4.4_df

v0.4.4_local_diff = methylSigCalc(
    meth = bs,
    groups = c(Treatment = 1, Control = 0),
    dispersion= 'both',
    local.disp = TRUE,
    winsize.disp = 200,
    local.meth = TRUE,
    winsize.meth = 200,
    min.per.group = c(2,2),
    T.approx = TRUE,
    num.cores = 1
)

v0.4.4_local_gr = GRanges(
    seqnames = v0.4.4_local_diff@data.chr,
    ranges = IRanges(
        start = v0.4.4_local_diff@data.start,
        end = v0.4.4_local_diff@data.start)
)

v0.4.4_local_df = data.frame(
    meth_case = v0.4.4_local_diff@results[, 'mu1'],
    meth_control = v0.4.4_local_diff@results[, 'mu0'],
    meth_diff = v0.4.4_local_diff@results[, 'meth.diff'],
    pvalue = v0.4.4_local_diff@results[, 'pvalue'],
    fdr = v0.4.4_local_diff@results[, 'qvalue'],
    disp_est = v0.4.4_local_diff@results[, 'theta'],
    log_lik_ratio = v0.4.4_local_diff@results[, 'logLikRatio'],
    df = v0.4.4_local_diff@results[, 'df'],
    stringsAsFactors = F
)

mcols(v0.4.4_local_gr) = v0.4.4_local_df

save(list = c('v0.4.4_gr', 'v0.4.4_local_gr'), file = 'rda/v0.4.4_msig.rda')
