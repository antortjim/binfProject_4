CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

tput setaf 2; echo "Filtering individuals in 1000Genomes"
# Generate file with individuals ID. So far let's select GBR
PANEL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Select all individuals
curl $PANEL | awk '$3 != "super_pop" {print}' > individuals.txt

# Reset the list of chromosomes that have been added
truncate -s 0 merge_list 


tput setaf 2; echo "Parsing VCF files"
# Extract genotypical data for all individuals and for SNPs run by AthGene
# Repeat for all chromosomes
for CHR in "${CHROMOSOMES[@]}"
do
  tput setaf 4; echo "Chromosome $CHR"
  # Find the SNPs in our data that belong to the current chromosome
  tail -n +2 metadata.txt | awk -F"\t" '$3 == '$CHR' {print $2}' | \
  grep -v NA | tr ',' '\n' > rs_chr${CHR}.txt

  # if rscode is missing, use the genomic position 
  #tail -n +2 metadata.txt | awk -F"\t" '$3 == '$CHR' && \
  #  $2 == "NA" {printf ("%s\t%s\n", $3, $4)}' > pos_chr$CHR.txt

  # Define the route to the vcfgz + vcfgz.tbi files
  VCFGZ_FILE="gzvcf/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

  # Open this vcfgz file, look at the snps listed in the file created above
  # keep the individuals in the individuals.txt list (can be all of them)
  # Outut a tped (text pedigree file)
  vcftools --gzvcf $VCFGZ_FILE \
    --snps rs_chr$CHR.txt \
    --keep <(cut -f1 individuals.txt) \
    --plink-tped --chr $CHR --out chr${CHR}_result

   # Transform the tped file in a binary pedigree file ready for plink
  plink --tfile chr${CHR}_result --make-bed --out chr${CHR}
  echo "chr${CHR}" >> merge_list
done
