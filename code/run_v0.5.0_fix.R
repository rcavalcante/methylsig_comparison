# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active umichbfxcore/r_344:latest
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')

library(bsseq)
devtools::load_all('/nfs/turbo/epicore-active/rcavalca/methylSig_v0.5.0')
library(tidyverse)

pdata = read_tsv('samplesheet_errbs.txt')
pdata$group = relevel(factor(pdata$group), ref = 'DS')
files = list.files(path = 'bedgraphs', full.names = TRUE)

bs = methylSigReadData(
    fileList = files,
    pData = pdata,
    assembly = NA,
    destranded = FALSE,
    maxCount = 500, minCount = 10,
    filterSNPs = FALSE,
    num.cores = 1,
    fileType = c("cov"),
    verbose = TRUE)
rownames(pData(bs)) = pData(bs)$sample

v0.5.0_fix_gr = methylSigCalc(
    meth = bs,
    comparison = 'group',
    dispersion= 'both',
    local.info = FALSE, local.winsize = 200,
    min.per.group = c(2,2),
    weightFunc=methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1
)

colnames(mcols(v0.5.0_fix_gr)) = c(
    'disp_est',
    'log_lik_ratio',
    'df',
    'meth_case',
    'meth_control',
    'meth_diff',
    'direction',
    'pvalue',
    'fdr')

v0.5.0_fix_local_gr = methylSigCalc(
    meth = bs,
    comparison = 'group',
    dispersion= 'both',
    local.info = TRUE, local.winsize = 200,
    min.per.group = c(2,2),
    weightFunc=methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1
)

colnames(mcols(v0.5.0_fix_local_gr)) = c(
    'disp_est',
    'log_lik_ratio',
    'df',
    'meth_case',
    'meth_control',
    'meth_diff',
    'direction',
    'pvalue',
    'fdr')

save(list = c('v0.5.0_fix_gr', 'v0.5.0_fix_local_gr'), file = 'rda/v0.5.0_fix_msig.rda')
