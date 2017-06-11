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

bash tidy_data.sh && \
Rscript read_snps.R $DATA_FOLDER && \
bash query_1000genomes.sh && \
bash merge_athgene_1000genomes.sh && \
bash plink_analysis.sh && \
Rscript visualization.R $DATA_FOLDER && \
bash admixture.sh && \
Rscript admixture.R && \
bash mail_test.sh
