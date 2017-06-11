cd $DATA_FOLDER

# rearranges the files so that data are sorted by chromosome
plink --bfile merged_clean --make-bed --out sorted

# TL;DR filter and select only variants with minor allele freq (maf) > 0.05
# A new set of plink files is generated from the data by just 
# selecting the variants that show a minor allele frequency or q
# bigger than 0.05. This means that only variants where at least 5 %
# of the samples bear the minor allele are kept and the others are discarded

# This is required because when the maf is very low and the sample size is not big enough,
# the estimation of the actual value for q is very uncertain
# What is the frequency of a very weird event? If I want to estimate this with
# low standard error, I need a huge sample size with lots of occurences of this event

#plink --make-bed --bfile sorted --geno 0.1 --maf 0.05 --out maf0.05
plink --make-bed --bfile sorted --maf 0.05 --out maf0.05

# prunes the snps and removes linkage desiquilibrium (LD)
plink --bfile maf0.05 --indep 50 5 2 --out maf0.05

## Analyze Identity By Descent individual wise
plink --bfile maf0.05 --extract maf0.05.prune.in --genome --min 0.1 --out maf0.05
grep -Fwvf <(tail -n +2 maf0.05.genome | awk '{print $3}' | sort | uniq) <(tail -n +2 maf0.05.genome | awk '{print $1}' | sort | uniq)> maf0.05.remove.fam
#plink --make-bed --bfile sorted --maf 0.05 --geno 0.1 --remove-fam maf0.05.remove.fam --out maf0.05
plink --make-bed --bfile sorted --maf 0.05 --remove-fam maf0.05.remove.fam --out maf0.05

# prunes the snps and removes linkage desiquilibrium (LD)
plink --bfile maf0.05 --indep 50 5 2 --out maf0.05


# perform pca on the snps that are not exluded
# (the ones in the maf0.05.prune.in file)
plink --bfile maf0.05 --extract maf0.05.prune.in --pca 200 --out maf0.05

## Output freq statistics
plink --bfile maf0.05 --extract maf0.05.prune.in --freq --out maf0.05

## Output heterozygosity statistics
plink --bfile maf0.05 --extract maf0.05.prune.in --het --out maf0.05


