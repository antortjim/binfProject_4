#! /bin/bash

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

cd $DATA_FOLDER/$GENOME
mkdir gzvcf
cd gzvcf

for CHR in ${CHROMOSOMES[@]}
do
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
done
