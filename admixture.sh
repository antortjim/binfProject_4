cd $DATA_FOLDER
ADMIXTURE=../code/admixture
#for i in 4 5 6;
#do echo "K=$i"; $ADMIXTURE --cv merged_clean.bed $i; done >> cvoutput
#grep -i 'CV error' cvoutput


awk '$2 == "ITU" || $2 == "CHB" || $2 == "CEU" || \
 $2 == "PEL" || $2 == "ACB" {print $1}' \
 <(tail -n +2 panel) \
 > representative_subset

paste representative_subset representative_subset > temp
mv temp representative_subset

cat $DATA_FOLDER/merged_clean.fam | grep ^I | cut -f1,2 \
 >> representative_subset



plink --bfile merged_clean \
      --keep representative_subset \
      --make-bed --out representative_subset

$DATA_FOLDER/$ADMIXTURE $DATA_FOLDER/representative_subset.bed 5

