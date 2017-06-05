# Merge the whole thing taking care of snps that
# are annotated using the rev strand

tput setaf 2; echo "Merging 1000Genomes and AthGene"
plink --tfile athgene --make-bed --out athgene
plink --bfile athgene --merge-list merge_list --out merged
mv merged.missnp old.missnp
plink --bfile athgene --exclude old.missnp --make-bed \
  --out athgene_clean

plink --bfile athgene_clean --merge-list merge_list --out merged_clean


for BATCH in {1..11}
do
  for IND in {1..48}
  do
    RUN=$(((BATCH-1)*48 + IND))
    printf "I$RUN\tAthGene\tAthGene\t0\t$BATCH\n" >> individuals.txt
  done
done
