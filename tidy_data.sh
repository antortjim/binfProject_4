#! /bin/bash

GTREPORT=$DATA_FOLDER/gtReport.txt
DATA=$DATA_FOLDER/data
MARKER_IDS=$DATA_FOLDER/marker_ids
RSIDS=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_b144_rsids.txt
ANNOTATED=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.annotated.txt
LOCUS_REPORT=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_LocusReport.txt
MANIFEST=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.csv
STRAND=$DATA_FOLDER/Exome_24/InfiniumCoreExome-24v1-1_A-b37.strand
PROBLEM=$DATA_FOLDER/problematic_marker_ids

# Extract marker_ids
tput setaf 2; echo "Extracting marker ids"
cut -f 1 $GTREPORT | tail -n +2  > $MARKER_IDS
cut -f 2- $GTREPORT | tail -n +2  > $DATA.tsv
cp $MARKER_IDS $DATA_FOLDER/all_marker_ids
cut -f 2- $DATA.csv > $DATA.tsv

tput setaf 2; echo "Filtering support files with available marker_ids"
####################################################
####################################################
# Subset the support files with the data available

# Get exm-number, chr, coord and strand from MANIFEST
tput setaf 2; echo "Manifest"
awk -F',' 'NR==FNR{c[$1]++;next};c[$2] > 0' \
   $MARKER_IDS \
   <(tail -n +9 $MANIFEST | head -n -24) \
   | cut -f2,11,10,21 -d',' | tr "," "\t" | sort -k1 > \
   $DATA_FOLDER/strand.txt


# NF = 4
# Join with annotated, locus_report and rsids
# Get genomic information (transcript, mutation...)
tput setaf 2; echo "annotated"
join -t$'\t' -e . $DATA_FOLDER/strand.txt \
    <(tail -n +2 $ANNOTATED | cut -f 1,5-8 | sort -k1) \
    > $DATA_FOLDER/joined
perl -p -i -e "s/\r//g" $DATA_FOLDER/joined

# NF = 4 + 4
# Get allele frequencies
tput setaf 2; echo "maf"
join -t$'\t' -e . $DATA_FOLDER/joined \
    <(tail -n +2 $LOCUS_REPORT | cut -f 2,5-9 | sort -k1) \
    > $DATA_FOLDER/joined_2

# NF = 8 + 5 = 13
# Get rscode
tput setaf 2; echo "rscode"
join -t$'\t' -e . $DATA_FOLDER/joined_2 <(tail -n +2 $RSIDS | sort -k1) \
    > $DATA_FOLDER/metadata.txt
perl -p -i -e "s/\r//g" $DATA_FOLDER/metadata.txt
# NF = 13 + 1 = 14


#rm $DATA_FOLDER/joined $DATA_FOLDER/joined_2
####################################################


tput setaf 2; echo "Filtering problematic markers"
####################################################
####################################################
# remove the marker ids that lack rscode
# remove the marker ids in chr 0
# remove INDELS
# remove snps with identical rscode but different exm-number

awk '$14 == "." || $2 == 0 || $1 ~ "IND" {print $1}' $DATA_FOLDER/metadata.txt \
    > $PROBLEM

awk '$14 != "." && $2 != 0 && $1 !~ "IND" {print}' $DATA_FOLDER/metadata.txt \
    > $DATA_FOLDER/metadata_clean.txt

# remove non first entries of snps with identical rscode but different exm
sort -u -k14 $DATA_FOLDER/metadata_clean.txt > $DATA_FOLDER/temp
mv $DATA_FOLDER/temp $DATA_FOLDER/metadata_clean.txt

grep "IND" $MARKER_IDS >> $PROBLEM
sort $PROBLEM | uniq > $DATA_FOLDER/temp
mv $DATA_FOLDER/temp $PROBLEM

grep -Fwvf $PROBLEM $MARKER_IDS > $DATA_FOLDER/out
cp $DATA_FOLDER/out $MARKER_IDS

# Index them starting at 1 and sort them alphabetically
awk '{printf("%d %s\n", NR, $0)}' $MARKER_IDS | sort -k 2 \
    > $DATA_FOLDER/marker_ids.indexed


tput setaf 2; echo "Filtering maf < 0.05"
# we are also filtering snps containing indels (there might be snps
# with ids that don't suggest it's an indel
awk '$13 > 0.05' $DATA_FOLDER/metadata_clean.txt \
    | sort -k1 \
    > $DATA_FOLDER/metadata_purged.txt

#rm $DATA_FOLDER/metadata.txt

tput setaf 2; echo "Filtering dataset"
######################################
## Select the marker_ids in the locus subset where maf > 0.05
cut -f1 $DATA_FOLDER/metadata_purged.txt > ${MARKER_IDS}_purged
cut -f14 $DATA_FOLDER/metadata_purged.txt > $DATA_FOLDER/rscodes

# 1
paste $DATA_FOLDER/all_marker_ids \
    <(tail -n +2 $DATA.tsv) \
    | sort -k 1 -t$'\t' \
    > $DATA_FOLDER/data_rownames.tsv

# 2
grep -Fwf ${MARKER_IDS}_purged $DATA_FOLDER/data_rownames.tsv \
    | sort -k 1 > $DATA_FOLDER/data_purged.tsv

# Generate tped
######################################

# 3
cut -f2- $DATA_FOLDER/data_purged.tsv | expand -t1 | tr "-" "0" \
    | sed -e 's/ //g' | sed -e 's/\(.\)/\1\t/g' \
    > $DATA_FOLDER/data.tped


cd $DATA_FOLDER
NMARKS=$(cat metadata_purged.txt | wc -l)
yes 0 | head -n $NMARKS > cM

# 4
paste <(cut -f2 metadata_purged.txt) \
      <(cut -f14 metadata_purged.txt) cM \
      <(cut -f3 metadata_purged.txt) \
      > left.tped

paste left.tped data.tped > athgene.tped
#rm cM left.tped data.tped

truncate -s 0 athgene_individuals.txt iid father phenotype athgene.tfam
# Generate tfam
######################################
NINDIV=528
for i in $(seq $NINDIV)
do
  echo "I$i" >> athgene_individuals.txt
  echo "1" >> iid
  echo "0" >> father
  echo "-9" >> phenotype
done

paste athgene_individuals.txt iid father father father phenotype \
    > athgene.tfam

rm iid father phenotype
