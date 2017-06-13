# Merge the whole thing taking care of snps that
# are annotated using the rev strand
FLIP=$1

cd $DATA_FOLDER
tput setaf 2; echo "Merging 1000Genomes and AthGene"

# Generate bed file from tped generated in R
plink --tfile athgene --make-bed --out athgene


# Flip alleles in the reverse strand!!!!!
if [ "$FLIP" = flip ]
then
  tput setaf 2; echo "Flipping negative strand snps"
  # Using illumina info
  awk '$4 == "-" {print $14}' metadata_purged.txt > negative_snps
  # Using http://www.well.ox.ac.uk/~wrayner/strand/
  #grep -Fwf \
  #  <(awk '$5 == "-" {print $1}' filtered.strand) <(cut -f1,14 metadata_purged.txt) \
  #  | cut -f2 | sort | uniq > negative_snps
  plink --bfile athgene --flip negative_snps --make-bed --out athgene_flipped

  cd $DATA_FOLDER/$GENOMES
  
  # Try to merge this file with the bed files parsed from the 1000Genomes gzvcfs
  plink --bfile $DATA_FOLDER/athgene_flipped \
      --merge-list merge_list \
      --out $DATA_FOLDER/merged

elif [ "$FLIP" = remove ]
then
  tput setaf 2; echo "Removing ambiguous snps"

  awk '$5 == "A" && $6 == "T" ||
       $5 == "T" && $6 == "A" ||
       $5 == "C" && $6 == "G" ||
       $5 == "G" && $6 == "C" {print $2}' \
       < athgene.bim \
       > negative_snps

       plink --bfile athgene --exclude negative_snps --make-bed \
             --out athgene_unambiguous

       cd $DATA_FOLDER/$GENOMES

       plink --bfile $DATA_FOLDER/athgene_unambiguous \
             --merge-list merge_list \
             --out $DATA_FOLDER/merged

       #### FIRST TRIAL
       cd $DATA_FOLDER
       mv merged.missnp old.missnp
       plink --bfile $DATA_FOLDER/athgene_unambiguous \
             --flip old.missnp \
             --make-bed \
             --out $DATA_FOLDER/merged
       
       cd $DATA_FOLDER/$GENOMES
       plink --bfile $DATA_FOLDER/merged \
           --merge-list merge_list \
           --out $DATA_FOLDER/merged_clean

       #### SECOND TRIAL 
       cd $DATA_FOLDER
        plink --bfile $DATA_FOLDER/athgene_unambiguous \
              --flip old.missnp \
              --exclude $DATA_FOLDER/merged_clean.missnp \
              --make-bed \
              --out $DATA_FOLDER/merged
       
       cd $DATA_FOLDER/$GENOMES
       plink --bfile $DATA_FOLDER/merged \
           --merge-list merge_list \
           --out $DATA_FOLDER/merged_clean
       

       cp $DATA_FOLDER/individuals.txt $DATA_FOLDER/individuals_athgene.txt
       
       for BATCH in {1..11}
       do
         for IND in {1..48}
         do
           RUN=$(((BATCH-1)*48 + IND))
           printf "I$RUN\tAthGene\tAthGene\t0\t$BATCH\n" >> $DATA_FOLDER/individuals_athgene.txt
         done
       done

       exit 0

else
  cd $DATA_FOLDER/$GENOMES
  # Try to merge this file with the bed files parsed from the 1000Genomes gzvcfs
  plink --bfile $DATA_FOLDER/athgene \
      --merge-list merge_list \
      --out $DATA_FOLDER/merged
fi


# Problematic snps will be discarded
cd $DATA_FOLDER
mv merged.missnp old.missnp

plink --bfile athgene --exclude old.missnp --make-bed \
  --out athgene_clean

cd $DATA_FOLDER/$GENOMES
plink --bfile $DATA_FOLDER/athgene_clean \
    --merge-list merge_list \
    --out $DATA_FOLDER/merged_clean

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
