cd /nfs/turbo/epicore-active/rcavalca/methylsig_comparison

mkdir -p covs

for file in `find ./bedgraphs -name "*bedGraph"`
do
    sample=`basename $file '_CpG.bedGraph'`
    awk -v OFS='\t' '{print $1"."$2, $1, $2, "*", $5 + $6, ($5/($5 + $6)) * 100, ($6/($5 + $6))*100}' $file > covs/${sample}.cov
done
