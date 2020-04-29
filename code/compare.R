# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active umichbfxcore/orochi_r:v1_manual
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')
library(GenomicRanges)

load('rda/new_msig.rda')
load('rda/old_msig.rda')

png('images/msig_new_old_meth_case.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_case, old_gr$meth.DR, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_no_meth_case.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_case, new_local_gr$meth_case, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_old_local_no_meth_case.png', height = 6, width = 6, units = 'in', res = 300)
    plot(old_gr$meth.DR, old_local_gr$meth.DR, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_old_local_meth_case.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_local_gr$meth_case, old_local_gr$meth.DR, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()



png('images/msig_new_old_meth_control.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_control, old_gr$meth.DS, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_no_meth_control.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_control, new_local_gr$meth_control, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_old_local_no_meth_control.png', height = 6, width = 6, units = 'in', res = 300)
    plot(old_gr$meth.DS, old_local_gr$meth.DS, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_old_local_meth_control.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_local_gr$meth_control, old_local_gr$meth.DS, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()


png('images/msig_new_old_meth_diff.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_diff, old_gr$meth.diff, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_no_meth_diff.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$meth_diff, new_local_gr$meth_diff, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_old_local_no_meth_diff.png', height = 6, width = 6, units = 'in', res = 300)
    plot(old_gr$meth.diff, old_local_gr$meth.diff, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_old_local_meth_diff.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_local_gr$meth_diff, old_local_gr$meth.diff, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()


png('images/msig_new_old_pvalue.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$pvalue, old_gr$pvalue, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_no_pvalue.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_gr$pvalue, new_local_gr$pvalue, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_old_local_no_pvalue.png', height = 6, width = 6, units = 'in', res = 300)
    plot(old_gr$pvalue, old_local_gr$pvalue, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()

png('images/msig_new_local_old_local_pvalue.png', height = 6, width = 6, units = 'in', res = 300)
    plot(new_local_gr$pvalue, old_local_gr$pvalue, pch = 20)
    abline(a = 0, b = 1, col = 'red')
dev.off()
