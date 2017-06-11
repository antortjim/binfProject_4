cd $DATA_FOLDER
ADMIXTURE=../code/admixture
#for i in 4 5 6;
#do echo "K=$i"; $ADMIXTURE --cv merged_clean.bed $i; done >> cvoutput
#grep -i 'CV error' cvoutput

$ADMIXTURE merged_clean.bed 5

