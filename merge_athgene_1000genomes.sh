# Merge the whole thing taking care of snps that
# are annotated using the rev strand

cd $DATA_FOLDER
tput setaf 2; echo "Merging 1000Genomes and AthGene"
plink --tfile athgene --make-bed --out athgene
plink --bfile $DATA_FOLDER/athgene --merge-list $GENOMES/merge_list --out $DATA_FOLDER/merged
mv merged.missnp old.missnp
plink --bfile athgene --exclude old.missnp --make-bed \
  --out athgene_clean

cd $GENOMES
plink --bfile $DATA_FOLDER/athgene_clean --merge-list merge_list --out $DATA_FOLDER/merged_clean


for BATCH in {1..11}
do
  for IND in {1..48}
  do
    RUN=$(((BATCH-1)*48 + IND))
    printf "I$RUN\tAthGene\tAthGene\t0\t$BATCH\n" >> $DATA_FOLDER/individuals.txt
  done
done

cd $DATA_FOLDER
