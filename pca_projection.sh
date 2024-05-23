#PCA
## make loadings of the ref
plink2 --bfile $OUT/${refname}_commonVars --freq counts --pca allele-wts --out $OUT/${refname}_pca

##----------project studies on the loadings
for study in $studies
do
	infile=$OUT/${study}_editedvars.clean.pruned_commonVars
	plink2 --bfile $infile \
	--read-freq $OUT/${refname}_pca.acount \
	--score $OUT/${refname}_pca.eigenvec.allele 2 5 header-read no-mean-imputation \
	variance-standardize \
	list-variants \
	--score-col-nums 6-15 \
	--out $OUT/${study}_${refname}_projected
done
