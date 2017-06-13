#! /bin/bash
export DATA_FOLDER="/home/antortjim/MEGA/AthGene/data"
export GTREPORT=$DATA_FOLDER/data.csv
export DATA=$DATA_FOLDER/data.tsv
export MARKER_IDS=$DATA_FOLDER/marker_ids
export RSIDS=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_b144_rsids.txt
export ANNOTATED=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.annotated.txt
export LOCUS_REPORT=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1_LocusReport.txt
export MANIFEST=$DATA_FOLDER/Exome_24/InfiniumExome-24v1-1_A1.csv
export STRAND=$DATA_FOLDER/Exome_24/InfiniumCoreExome-24v1-1_A-b37.strand
export PROBLEM=$DATA_FOLDER/problematic_marker_ids
export FLIP="remove"

# Uncomment if sql database available,
# otherwise stick to the supplied rs_code_fitness.txt in data
#python extract_rscodes_fitness.py $DATA_FOLDER && \

# Tidy data and prepare to merge with 1000Genomes
bash tidy_data.sh && \

bash query_1000genomes.sh && \

# Merge with 1000Genomes
bash merge_athgene_1000genomes.sh $FLIP

# Run plink on the merged data
#bash plink_analysis.sh && \

## Visualize
#Rscript visualization.R $DATA_FOLDER && \

#bash admixture.sh && \
#Rscript admixture.R && \
#bash mail_test.sh
