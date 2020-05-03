# docker run -ti -v /nfs/turbo/epicore-active:/nfs/turbo/epicore-active umichbfxcore/orochi_r:v1_manual
# R

setwd('/nfs/turbo/epicore-active/rcavalca/methylsig_comparison')
library(GenomicRanges)

load('rda/v0.99.4_msig.rda')
load('rda/v0.5.0_msig.rda')
load('rda/v0.4.4_msig.rda')

# Are all the starts the same?
# If they are, don't need to do findOverlaps()
all(start(v0.99.4_gr) == start(v0.5.0_gr))
all(start(v0.5.0_gr) == start(v0.4.4_gr))

all(start(v0.99.4_local_gr) == start(v0.5.0_local_gr))
all(start(v0.5.0_local_gr) == start(v0.4.4_local_gr))

results = list(
    v0.99.4 = list(
        no_local = v0.99.4_gr,
        local = v0.99.4_local_gr
    ),
    v0.5.0 = list(
        no_local = v0.5.0_gr,
        local = v0.5.0_local_gr
    ),
    v0.4.4 = list(
        no_local = v0.4.4_gr,
        local = v0.4.4_local_gr
    )
)

comparisons = expand.grid(
    version1 = c('v0.99.4','v0.5.0','v0.4.4'),
    version2 = c('v0.99.4','v0.5.0','v0.4.4'),
    test1 = c('no_local', 'local'),
    test2 = c('no_local', 'local'),
    stringsAsFactors = F)

within_versions = subset(comparisons, version1 == version2 & test1 != test2)[1:3, ]
within_tests = subset(comparisons, version1 != version2 & test1 == test2)[c(1,2,4,7,8,10), ]
comparisons = rbind(within_versions, within_tests)

for(i in seq(nrow(comparisons))) {
    version1 = comparisons[i, 'version1']
    version2 = comparisons[i, 'version2']
    test1 = comparisons[i, 'test1']
    test2 = comparisons[i, 'test2']

    xlab_str = sprintf('%s %s', version1, test1)
    ylab_str = sprintf('%s %s', version2, test2)

    title = sprintf('%s %s versus %s %s', version1, test1, version2, test2)
    file = gsub(' ', '_', title)

    result1 = results[[version1]][[test1]]
    result2 = results[[version2]][[test2]]

    png(sprintf('images/%s_meth_case.png', file), height = 6, width = 6, units = 'in', res = 300)
        plot(result1$meth_case, result2$meth_case,
            pch = 20,
            main = sprintf('%s meth_case', title),
            xlab = sprintf('%s meth_case', xlab_str),
            ylab = sprintf('%s meth_case', ylab_str))
        abline(a = 0, b = 1, col = 'red')
    dev.off()

    png(sprintf('images/%s_meth_control.png', file), height = 6, width = 6, units = 'in', res = 300)
        plot(result1$meth_control, result2$meth_control,
            pch = 20,
            main = sprintf('%s meth_control', title),
            xlab = sprintf('%s meth_control', xlab_str),
            ylab = sprintf('%s meth_control', ylab_str))
        abline(a = 0, b = 1, col = 'red')
    dev.off()

    png(sprintf('images/%s_meth_diff.png', file), height = 6, width = 6, units = 'in', res = 300)
        plot(result1$meth_diff, result2$meth_diff,
            pch = 20,
            main = sprintf('%s meth_diff', title),
            xlab = sprintf('%s meth_diff', xlab_str),
            ylab = sprintf('%s meth_diff', ylab_str))
        abline(a = 0, b = 1, col = 'red')
    dev.off()

    png(sprintf('images/%s_pvalue.png', file), height = 6, width = 6, units = 'in', res = 300)
        plot(result1$pvalue, result2$pvalue,
            pch = 20,
            main = sprintf('%s pvalue', title),
            xlab = sprintf('%s pvalue', xlab_str),
            ylab = sprintf('%s pvalue', ylab_str))
        abline(a = 0, b = 1, col = 'red')
    dev.off()

}
