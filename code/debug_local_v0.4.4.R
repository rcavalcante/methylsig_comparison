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

########################################

meth = bs
groups = c(Treatment = 1, Control = 0)
dispersion= 'both'
local.disp = TRUE
winsize.disp = 200
local.meth = TRUE
winsize.meth = 200
min.per.group = c(2,2)
T.approx = TRUE
weightFunc=methylSig_weightFunc
num.cores = 1

########################################

min.disp=1e-6


if(meth@resolution == "tfbs") {
    local.disp=FALSE
    local.meth=FALSE
}

if(local.meth == FALSE) winsize.meth = 0
if(local.disp == FALSE) winsize.disp = 0

treatment = slot(meth,"treatment")

group1 = which(treatment == groups[1])
group2 = which(treatment == groups[2])

if(length(group1) == 0 || length(group2) == 0) {
    stop("Groups do not match your treatment in the input data")
}

if(length(min.per.group) == 1) {
    min.per.group = c(min.per.group,min.per.group)
}

methSigObject=new.env(parent=globalenv())
class(methSigObject)='pointer'
methSigObject$groups <- list(group1 = group1,
           group2 = group2,
           group3 = c(group1,group2)
          )

orderMethStart = order(meth@data.ids)

methSigObject$treads   = t(meth@data.numTs[orderMethStart,])
methSigObject$creads   = t(meth@data.numCs[orderMethStart,])
methSigObject$treads[is.na(methSigObject$treads)] = 0
methSigObject$creads[is.na(methSigObject$creads)] = 0

methSigObject$uniqueLoc = meth@data.ids[orderMethStart]
methSigObject$stepSize = ifelse(meth@destranded, 2, 1)
methSigObject$wMeth = winsize.meth
methSigObject$wMethIndex = winsize.meth/methSigObject$stepSize
methSigObject$wMethNorm = 1/(methSigObject$stepSize*(methSigObject$wMethIndex+1))
methSigObject$wDispersion = winsize.disp
methSigObject$wDispIndex = winsize.disp/methSigObject$stepSize
methSigObject$wDispNorm  = 1/(methSigObject$stepSize*(methSigObject$wDispIndex+1))
methSigObject$min.InvDisp = 0.001
methSigObject$max.InvDisp = max(1/max(min.disp,1e-6), methSigObject$min.InvDisp)
methSigObject$numValidMu = min.per.group
methSigObject$weightFunc = weightFunc

if(dispersion == "both") {
   methSigObject$dispersionGroups = c(group1,group2)
   methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) - 1, 0) + pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group1,]) - 1, 0)
}  else if(dispersion == groups[1] || (length(names(groups)[1])>0 && dispersion == names(groups)[1])) {
    methSigObject$dispersionGroups = group1
    methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group1,]) - 1, 0)
} else if(dispersion == names(groups)[2] || dispersion == groups[2]) {
    methSigObject$dispersionGroups = group2
    methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) - 1, 0)
} else {
    cat("Dispersion should be \"", names(groups)[1], "\", \"", names(groups)[2], "\" or \"both\".\n", sep="")
    return(NULL)
}

nLoci = NCOL(methSigObject$creads)
methSigObject$muEst <- matrix(0, ncol=nLoci, nrow=NROW(methSigObject$creads))

muList1 <- colSums(methSigObject$creads[group1,])/(colSums((methSigObject$creads+methSigObject$treads)[group1,])+1e-100)
muList2 <- colSums(methSigObject$creads[group2,])/(colSums((methSigObject$creads+methSigObject$treads)[group2,])+1e-100)

for(g in group1) {
    methSigObject$muEst[g,] = muList1
}
for(g in group2) {
    methSigObject$muEst[g,] = muList2
}

validLoci = ((colSums((methSigObject$creads + methSigObject$treads> 0)[group1,]) >= min.per.group[1])
           & (colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) >= min.per.group[2]))

methSigObject$whichOrd = cumsum(validLoci)

########################################

vLociIdx = which(validLoci)

loc = vLociIdx[1]
obj = methSigObject

minMu = 0
maxMu = 1
group1=obj$groups[[1]]
group2=obj$groups[[2]]

allGroupsIndex = 3

locSize = NCOL(obj$creads)

validMuList  <- max(1,loc-obj$wMethIndex):min(locSize,loc+obj$wMethIndex)
validPhiList <- max(1,loc-obj$wDispIndex):min(locSize,loc+obj$wDispIndex)

validMuList  = validMuList [which(abs(obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc]) <= obj$wMeth)]
validPhiList = validPhiList[which(abs(obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc]) <= obj$wDispersion)]

whichUseful = which(obj$validForPhiCalculate[validPhiList] > 0)
if(length(whichUseful) == 0)  return(c(loc,NA,NA,NA,NA,NA,NA))

validPhiList = validPhiList[whichUseful]

if(length(validPhiList) > 5) validPhiList = validPhiList[order(abs(obj$uniqueLoc[validPhiList]  - obj$uniqueLoc[loc]))[1:5]]

weightPhi <- obj$weightFunc((obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc])*obj$wDispNorm)
weightMu  <- obj$weightFunc((obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc])*obj$wMethNorm)
