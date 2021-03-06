# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active umichbfxcore/orochi_r:v1_manual
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')

library(bsseq)
library(methylSig)
library(tidyverse)

pdata = read_tsv('samplesheet_errbs.txt')
files = list.files(path = 'bedgraphs', full.names = TRUE)

bs = bsseq::read.bismark(
    files = files,
    colData = pdata,
    rmZeroCov = FALSE,
    strandCollapse = FALSE)
rownames(pData(bs)) = pData(bs)$sample

bs = filter_loci_by_coverage(bs = bs, min_count = 10, max_count = 500)

bs = filter_loci_by_group_coverage(bs = bs, group_column = 'group', min_samples_per_group = c(DR = 2, DS = 2))

v0.99.4_gr = diff_methylsig(
    bs = bs,
    group_column = 'group',
    comparison_groups = c(case = 'DR', control = 'DS'),
    disp_groups = c(case = TRUE, control = TRUE),
    local_window_size = 0,
    t_approx = TRUE,
    n_cores = 1)

v0.99.4_local_gr = diff_methylsig(
    bs = bs,
    group_column = 'group',
    comparison_groups = c(case = 'DR', control = 'DS'),
    disp_groups = c(case = TRUE, control = TRUE),
    local_window_size = 200,
    t_approx = TRUE,
    n_cores = 1)

save(list = c('v0.99.4_gr', 'v0.99.4_local_gr'), file = 'rda/v0.99.4_msig.rda')
