cd $DATA_FOLDER
ADMIXTURE=../code/admixture
#for i in 4 5 6;
#do echo "K=$i"; $ADMIXTURE --cv merged_clean.bed $i; done >> cvoutput
#grep -i 'CV error' cvoutput


awk '$2 = "ITU | $2 = "CHB" | $2 = "CEU" | \
 $2 = "PEL" | $2 = "ACB" {print $1}' panel > representative_subset
plink --bfile merged_clean.bed --keep representative_subset \
  --out representative_subset

$ADMIXTURE representative_subset.bed 5

