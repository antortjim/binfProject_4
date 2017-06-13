CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
#CHROMOSOMES=(22)
GENOMES=$DATA_FOLDER/1000_genomes
cd $DATA_FOLDER
mkdir $GENOMES

tput setaf 2; echo "Filtering individuals in 1000Genomes"
# Generate file with individuals ID. So far let's select GBR
PANEL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

wget $PANEL -O $DATA_FOLDER/panel
# Select the individuals
## CEU (EUR) PEL (AMR) CHB (EAS) ITU (EAS) ACB (AFR, African Caribbean in Barbados)
#cat $DATA_FOLDER/panel | awk '$2 == "CEU" || $2 == "PEL" || $2 == "CHB" \
#                || $2 == "ITU" || $2 == "ACB" {print}' > $GENOMES/individuals.txt

cp $DATA_FOLDER/panel $DATA_FOLDER/individuals.txt

# Reset the list of chromosomes that have been added
truncate -s 0 $GENOMES/merge_list 


tput setaf 2; echo "Parsing VCF files"
# Extract genotypical data for all individuals and for SNPs run by AthGene
# Repeat for all chromosomes
for CHR in "${CHROMOSOMES[@]}"
do
  tput setaf 4; echo "Chromosome $CHR"
  # Find the SNPs in our data that belong to the current chromosome
  cat $DATA_FOLDER/metadata_purged.txt | awk -F"\t" '$2 == '$CHR' {print $14}' | \
  grep -v NA | tr ',' '\n' > $GENOMES/rs_chr${CHR}.txt

  # Define the route to the vcfgz + vcfgz.tbi files
  VCFGZ_FILE="$GENOMES/gzvcf/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

  # Open this vcfgz file, look at the snps listed in the file created above
  # keep the individuals in the individuals.txt list (can be all of them)
  # Outut a tped (text pedigree file)
  vcftools --gzvcf $VCFGZ_FILE \
    --snps $GENOMES/rs_chr$CHR.txt \
    --keep <(cut -f1 $DATA_FOLDER/individuals.txt) \
    --plink-tped --chr $CHR --out $GENOMES/chr${CHR}_result

   # Transform the tped file in a binary pedigree file ready for plink
  plink --tfile $GENOMES/chr${CHR}_result --make-bed --out $GENOMES/chr${CHR}
  echo "chr${CHR}" >> $GENOMES/merge_list
done
