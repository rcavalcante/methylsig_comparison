# module load singularity/3.5.2
# singularityshell r_344
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')

library(bsseq)
library(methylSig)
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

old_gr = methylSigCalc(
    meth = bs,
    comparison = 'group',
    dispersion= 'both',
    local.info = FALSE, local.winsize = 200,
    min.per.group = c(2,2),
    weightFunc=methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1
)

old_local_gr = methylSigCalc(
    meth = bs,
    comparison = 'group',
    dispersion= 'both',
    local.info = TRUE, local.winsize = 200,
    min.per.group = c(2,2),
    weightFunc=methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1
)

save(list = c('old_gr', 'old_local_gr') file = 'rda/old_msig.rda')
