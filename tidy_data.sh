#! /bin/bash

DATA_FOLDER="../data"
GTREPORT=$DATA_FOLDER/data.csv
DATA=$DATA_FOLDER/data.tsv
MARKER_IDS=$DATA_FOLDER/marker_ids
RSIDS=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_b144_rsids.txt
ANNOTATED=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.annotated.txt
LOCUS_REPORT=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_LocusReport.txt
MANIFEST=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.csv
STRAND=$DATA_FOLDER/Exome_24/InfiniumCoreExome-24v1-1_A-b37.strand
PROBLEM=$DATA_FOLDER/problematic_marker_ids

# Extract marker_ids
tput setaf 2; echo "Extracting marker ids"
cp $MARKER_IDS $DATA_FOLDER/all_marker_ids
cut -f 2- $GTREPORT > $DATA


tput setaf 2; echo "Filtering problematic markers"
# remove the marker ids that lack rscode
awk '$2 == "." {print $1}' $RSIDS > $PROBLEM
# remove the marker ids in chr 0
awk '$2 == 0 {print $1}' $ANNOTATED >> $PROBLEM

# remove INDELS
grep IND $MARKER_IDS >> $DATA_FOLDER/problematic_marker_ids

grep -Fwvf $DATA_FOLDER/problematic_marker_ids $MARKER_IDS > $DATA_FOLDER/out
mv $DATA_FOLDER/out $MARKER_IDS

# Index them starting at 1 and sort them alphabetically
awk '{printf("%d %s\n", NR, $0)}' $MARKER_IDS | sort -k 2 > $DATA_FOLDER/marker_ids.indexed

tput setaf 2; echo "Filtering support files with available marker_ids"
# Subset the support files with the data available

# Get exm-number, chr, coord and strand
awk -F',' 'NR==FNR{c[$1]++;next};c[$2] > 0' \
   <(cut -f 2 -d' ' $DATA_FOLDER/marker_ids.indexed) \
   <(tail -n +9 $MANIFEST | head -n -24) | cut -f2,11,10,21 -d',' | tr "," "\t" \
   > $DATA_FOLDER/strand.txt

# NF = 4

# Join sowyer with annotated, locus_report and rsids
# Get genomic information (transcript, mutation...)
join -t$'\t' -e . <(sort $DATA_FOLDER/strand.txt) <(tail -n +2 $ANNOTATED | cut -f 1,5-8 | sort -k 1) > $DATA_FOLDER/joined
perl -p -i -e "s/\r//g" $DATA_FOLDER/joined 

# NF = 4 + 4
# Get allele frequencies
join -t$'\t' -e . $DATA_FOLDER/joined <(tail -n +2 $LOCUS_REPORT | cut -f 2,5-9 | sort -k 1) > $DATA_FOLDER/joined_2

# NF = 8 + 5 = 13
# Get rscode
join -t$'\t' -e . $DATA_FOLDER/joined_2 <(tail -n +2 $RSIDS | sort -k1) > $DATA_FOLDER/joined_3
perl -p -i -e "s/\r//g" $DATA_FOLDER/joined_3 

# NF = 13 + 1 = 14
rm $DATA_FOLDER/joined $DATA_FOLDER/joined_2

tput setaf 2; echo "Filtering maf < 0.05"
# we are also filtering snps containing indels (there might be snps
# with ids that don't suggest it's an indel
awk '$13 > 0.05' $DATA_FOLDER/joined_3 | tail -n +2 > $DATA_FOLDER/metadata_subset.txt
rm $DATA_FOLDER/joined_3


## Select the marker_ids in the locus subset where maf > 0.05
awk '{print $1}' $DATA_FOLDER/metadata_subset.txt > $MARKER_IDS.frequent

#diff <(cat marker_ids.shared | cut -f 2 -d' ') <(cat locus_subset.txt | cut -f 2)

## Filter the marker_ids.indexed list using the marker_ids.shared list
## and keep the row id
## These are the rows we want to filter out in the main dataset
## becuase they correspond to snps that
#
#  # appear in the support files
#  # their maf in the support file is > 0.05
#
## Rows shared is a list of row numbers that contain data for markers available
## in the support files and with MAF > 0.05
## The numbers are sorted so that their marker ids are sorted alphabetically
grep -Fwf $MARKER_IDS.frequent $MARKER_IDS.indexed | sort -k 2 | cut -f 1 -d' ' > $DATA_FOLDER/selected_rows

#
## Select id of good markers
## Index 0 to account for header
#tput setaf 2; echo "Filtering dataset"
paste $DATA_FOLDER/all_marker_ids $DATA | sort -k 1 -t$'\t' > $DATA_FOLDER/data_rownames.tsv
## sort the rows of the data according to the marker ids
grep -Fwf $MARKER_IDS.frequent $DATA_FOLDER/data_rownames.tsv | sort -k 1 > $DATA_FOLDER/data_purged.tsv


### Add rscode to variants that lack it in the support files by taking it from SeattleSNPs
## This command looks up the chr and the coord of the snps missing rscode and uses that to query 
## the seattlesnp and returns the positions queried along with the found rscode
#
##grep -Fwf <(grep -Fwf <(awk -F"\t" '$2 !~ "rs"' nomenclature_subset.txt | cut -f 1) annotation_subset.txt | awk '$2 != 0 {printf ("%s\t%s\n", $2, $3)}') /home/antortjim/MEGA/AthGene/data/seatle_snps/SeattleSeqAnnotation138.seattle_custom.395898255493.txt | cut -f 2,3,11 | uniq | sort -k1 -k2 > found_variants
##
##paste found_variants <(grep -Fwf <(awk '{printf "%s-%s\n", $1, $2}' found_variants) <(sort -k2 -k3 -t$'\t' annotation_subset.txt | awk '{printf "%s-%s\t%s\n", $2, $3, $1}') | awk '{print $2}') | cut -f 3,4 > out
##grep -vE "^0" out | sed 's/^/rs/' | awk '{printf "%s\t%s\n", $2, $1}' > found_variants
##rm out

# transpose data
tput setaf 2; echo "Transposing table"
bash transpose.sh $DATA_FOLDER/data_purged.tsv > $DATA_FOLDER/data_purged_transposed.tsv
