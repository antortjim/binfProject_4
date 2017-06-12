# Merge the whole thing taking care of snps that
# are annotated using the rev strand

cd $DATA_FOLDER
tput setaf 2; echo "Merging 1000Genomes and AthGene"

# Generate bed file from tped generated in R
plink --tfile athgene --make-bed --out athgene

# Flip alleles in the reverse strand!!!!!
#rm athgene_flipped.bim
#awk '$4 == "-" {print $14}' metadata_subset.txt > negative_snps
#awk -F$'\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' negative_snps athgene.bim > mysubset.bim
#
##While loop to read line by line
#while read -r LINE
#do
#    FLIPPED=$(echo $LINE | tr A "|" | tr C "@" | tr G "#" | tr T "~" | \
#               tr "|" T | tr "@" G  | tr "#" C | tr "~" A)
#
#    echo $FLIPPED >> athgene_flipped.bim
#done < mysubset.bim
#
#awk '$4 == "+" {print $14}' metadata_subset.txt > positive_snps
#awk -F$'\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' positive_snps athgene.bim > athgene_notflipped.bim
#
#mv athgene.bim old.bim
#cat athgene_notflipped.bim athgene_flipped.bim | sort -g -k1,4 | sed "s/ \+/\t/g" > athgene.bim

cd $DATA_FOLDER/$GENOMES

# Try to merge this file with the bed files parsed from the 1000Genomes gzvcfs
plink --bfile $DATA_FOLDER/athgene --merge-list merge_list --out $DATA_FOLDER/merged

# Problematic snps will be discarded
cd $DATA_FOLDER
mv merged.missnp old.missnp

plink --bfile athgene --exclude old.missnp --make-bed \
  --out athgene_clean

cd $DATA_FOLDER/$GENOMES
plink --bfile $DATA_FOLDER/athgene_clean --merge-list merge_list --out $DATA_FOLDER/merged_clean
plink --bfile $DATA_FOLDER/merged_clean --exclude $DATA_FOLDER/old.missnp --make-bed --out $DATA_FOLDER/merged_clean


cp $DATA_FOLDER/individuals.txt $DATA_FOLDER/individuals_athgene.txt

for BATCH in {1..11}
do
  for IND in {1..48}
  do
    RUN=$(((BATCH-1)*48 + IND))
    printf "I$RUN\tAthGene\tAthGene\t0\t$BATCH\n" >> $DATA_FOLDER/individuals_athgene.txt
  done
done

cd $DATA_FOLDER
