# Merge the whole thing taking care of snps that
# are annotated using the rev strand

cd $DATA_FOLDER
tput setaf 2; echo "Merging 1000Genomes and AthGene"

# Generate bed file from tbed generated in R
plink --tfile athgene --make-bed --out athgene
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


cp individuals.txt individuals_athgene.txt

for BATCH in {1..11}
do
  for IND in {1..48}
  do
    RUN=$(((BATCH-1)*48 + IND))
    printf "I$RUN\tAthGene\tAthGene\t0\t$BATCH\n" >> individuals_athgene.txt
  done
done

cd $DATA_FOLDER
