arguments <- commandArgs(trailingOnly = T)
data.dir <- arguments[[1]]
print(getwd())
suppressMessages(library("pegas"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("plyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))
suppressMessages(library("beepr"))

flipper <- c("A", "C", "G", "T")
names(flipper) <- rev(flipper)

flip_strand <- function(myRow) {

    # transforms metadata of markers in reverse strand (BOT) so
    # that they are displayed as if they were in the forward (TOP)
    if(myRow["strand"] == "-") {
      myRow["allele_1"] <- flipper[myRow["allele_1"]]
      myRow["allele_2"] <- flipper[myRow["allele_2"]]
    }
  return(myRow)
}


message("Reading metadata")
# read genes metadata
fl <- paste(data.dir, "metadata_subset.txt", sep = "/")
metadata <- read.table(file = fl, na.strings = c("."),
                       stringsAsFactors = F,
                       col.names = c("illumina_code", "chr", "coord",
                                     "strand", "transcript",
                                     "gene", "in-exon", "mutation", "freqAA",
                                     "freqAB", "freqBB", "call-freq", "MAF", "rs_code"))



metadata$cM <- 0

metadata %>% arrange(chr) %>% write.table(file = "metadata.txt",
                                          quote = F, row.names = F, sep = "\t")

message("Reading data")
# read data from AthGene
data.fl <- paste(data.dir, "data_purged_transposed.tsv", sep = "/")
id.fl <- paste(data.dir, "marker_ids.frequent", sep = "/")
index.fl <- paste(data.dir, "selected_rows", sep = "/")

id <- read.table(file = id.fl, stringsAsFactors = F)[,1]
index <- read.table(file = index.fl)[,1]

column_names <- metadata %>% apply(MARGIN = 1, function(x) {
  x[c("illumina_code", "chr", "rs_code", "cM", "coord")] %>% as.character %>% paste(collapse = ":") }) %>%
  gsub(pattern = " ", replacement = "")
n_batches <- 11
batches <- rep((1:n_batches), each = 48) %>% as.factor

df <- read.table(file = data.fl, sep = " ", header = T,
                 na.strings =  "--", stringsAsFactors = F, check.names = F,
                 col.names = column_names, nrows = 48 * n_batches)

#df <- df[,1:2]
message("Generating genind object")
genind_df <- df2genind(df, sep = "")
# takes a while

message("Generating loci object")
loci_df <- as.loci(genind_df)

# saveRDS(df, "df.rds")
# saveRDS(genind_df, "genind_df.rds")
# saveRDS(loci_df, "loci_df.rds")
# 
# loci_df <- readRDS("loci_df.rds")

message("Exporting tped and tfam files")
tped.fl <- paste(data.dir, "athgene", sep = "/")
tfam.fl <- paste(data.dir, "athgene.tfam", sep = "/")
transpose.script <- paste(data.dir, "../code/transpose.sh", sep = "/")
write.loci(x = loci_df, file = tped.fl, loci.sep = "\t", allele.sep = ":",
           quote = F, row.names = F, na = "0:0")
system(paste('bash ', transpose.script, " ", tped.fl, ' | tr ": " "\t" | cut -f 2- > ', tped.fl, '.tped', sep = ""))

tfam <- data.frame(fid = paste("I", 1:(48 * n_batches), sep = ""),
                   iid = 1,
                   father = 0,
                   mother = 0,
                   sex = 0,
                   phenotype = -9)

write.table(tfam, file = tfam.fl, quote = F, col.names = F, row.names = F)
beep()
