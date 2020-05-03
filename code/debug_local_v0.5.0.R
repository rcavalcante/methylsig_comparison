# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active umichbfxcore/r_344:latest
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

########################################

meth = bs
comparison = 'group'
dispersion= 'both'
local.info = TRUE
local.winsize = 200
min.per.group = c(2,2)
weightFunc=methylSig_weightFunc
T.approx = TRUE
num.cores = 1

########################################

if(!local.info) {
    local.winsize = 0
}

if(length(min.per.group) == 1) {
    min.per.group = c(min.per.group,min.per.group)
}

min.disp = 1e-6
min.InvDisp = 0.001
max.InvDisp = max(1/max(min.disp, 1e-6), min.InvDisp)

minMu = 0
maxMu = 1

#####################################
# Get the group labels, NOTE: THIS ASSUMES CORRECT REFERENCE LEVEL SET
pdata = bsseq::pData(meth)
group2 = levels(pdata[, comparison])[2]
group1 = levels(pdata[, comparison])[1]

# Determine which rows of pData belong to which group
# / which columns of Cov and M matrices belong to which group
group2_idx = which(pdata[,comparison] == group2)
group1_idx = which(pdata[,comparison] == group1)

# Determine which sample column indexes to use for dispersion calculation
if(dispersion == 'both') {
    disp_groups_idx = c(group2_idx, group1_idx)
} else if (dispersion == group2) {
    disp_groups_idx = group2_idx
} else if (dispersion == group1) {
    disp_groups_idx = group1_idx
} else {
    stop('"dispersion" should be one of "both", the name of group2, or the name of group1')
}

#####################################
# Determine which sites are valid to test according to min.per.group
all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))

# Estimate mu per locus within each group. The same value is used for all samples within the same group.
muEst = matrix(0, ncol = ncol(meth), nrow = nrow(meth))
muEst[, group2_idx] = base::rowSums(all_meth[, group2_idx]) / (base::rowSums(all_cov[, group2_idx]) + 1e-100)
muEst[, group1_idx] = base::rowSums(all_meth[, group1_idx]) / (base::rowSums(all_cov[, group1_idx]) + 1e-100)

# Determine which loci satisfy min.per.group
valid_idx = which(
    base::rowSums(all_cov[, group2_idx] > 0) >= min.per.group[2] & base::rowSums(all_cov[, group1_idx] > 0) >= min.per.group[1]
)

#####################################
# Resize all_cov, all_meth, and muEst to valid_idx
# Extract the granges of the meth BSseq object
# These are all used within the mclapply below
# all_cov = all_cov[valid_idx,]
# all_meth = all_meth[valid_idx,]
# muEst = muEst[valid_idx,]
# meth_gr = GenomicRanges::granges(meth)[valid_idx,]

num_loci = length(valid_idx)

########################################

#locus_idx = seq_along(valid_idx)[1]
locus_idx = valid_idx[1]
locus_idx_2 = valid_idx[1]

# NOTE: Here, you're getting local_loci_idx on the basis of meth_gr, which is granges(meth)
# but filtered for valid loci

# Get the indices which are within the local.winsize, but also limit to 5 CpGs on either side
local_loci_idx = intersect(
    which(abs(BiocGenerics::start(meth_gr)[locus_idx] - BiocGenerics::start(meth_gr)) < local.winsize),
    max(1, locus_idx - 5):min(num_loci, locus_idx + 5))

# What happens if you change meth_gr to granges(meth) above

local_loci_idx = intersect(
    which(abs(BiocGenerics::start(granges(meth))[locus_idx] - BiocGenerics::start(granges(meth))) < local.winsize),
    max(1, locus_idx - 5):min(num_loci, locus_idx + 5))

########################################


local_loci_norm = (BiocGenerics::start(granges(meth))[local_loci_idx] - BiocGenerics::start(granges(meth))[locus_idx]) / (local.winsize + 1)

# Calculate the weights
# Each is a vector of values of the weight function (range)
local_weights = weightFunc(local_loci_norm)

# Collect Cov and M matrices for all the loci in the window
# Rows are loci and columns are samples
local_cov = all_cov[local_loci_idx, ]
local_meth = all_meth[local_loci_idx, ]

# Collect the correct rows of muEst
local_muEst = muEst[local_loci_idx, ]
