rm(list = ls())

pca <- read.table("/home/antortjim/MEGA/AthGene/last/pca_ancestry.txt",
                  col.names = c("sample", "super"))
admixture <- read.table("/home/antortjim/MEGA/AthGene/last/admixture_ancestry.txt",
                        col.names = c("sample", "super"))

data <- merge(pca, admixture, by = "sample")
agreement <- data[,2] == data[,3]
table(agreement)

data[!agreement,]
