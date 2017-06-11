for i in 4 5 6;
do echo "K=$i"; admixture --cv ../merged_clean.bed $i; done >> cvoutput
grep -i 'CV error' cvoutput


admixture ../merged_clean.bed 5
