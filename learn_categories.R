data.dir  <- "~/MEGA/AthGene/data"
#data.dir <- commandArgs(trailingOnly = T)[[1]]
setwd(data.dir)
library("readxl")
library("plyr")
library("dplyr")

excel.fl <- paste(data.dir, "category.xlsx", sep = "/")
data <- read_xlsx(path = excel.fl)
old <- data
colnames(data)[11:13] <- c("maj/maj",	"maj/min", "min/min")

toMatch <- c("sc\\.", "impact", "Rel.", "Cit", "art", "X__4")
           
remove_columns <- grep(pattern = paste(toMatch,collapse="|"),
                       ignore.case = T,
                       x = colnames(data))

data <- data[-1, -remove_columns]

colnames(data)[c(2, 3, 4, 8, 9, 14, 18)] <- c("gene", "gene.name", "protein.type", "chip.genotype",
                                              "ensembl", "confidence", "eur.freq")

categories <- c("Alcohol", "Appetite", "Bone density", "Caffeine", "Calcium",
                "Carbohydrates", "Cravings/Emotional eating/Binge eating",
                "Endurance", "Endurance vs power", "Fat", "Gluten intolerance",
                "Heart rate trainability", "Heart Health/Cardiovascular Disease",
                "High altitude training", "Inflammation", "Injury CT", "Injury muscles",
                "Iron", "Irritable bowel syndrome", "Lactose intolerance", "Memory",
                "Migraine", "Motivation to train", "Muscle fatigue", "Muscle fiber composition",
                "Muscle growth", "Omega 3/6", "Power", "Protein metabolism - Homocysteine levels",
                "Protein metabolism - Thyrosinemia", "Protein Uptake", "Response to diet", "Salt",
                "Skin (Aging)", "Skin (Cellulite)", "Skin (Scarring)", "Skin (Sun)", "Sleep Quality",
                "Sweating/body odour", "Sweet tooth", "Testosterone (male only)", "Thrill seeking", 
                "Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E", "VO2max")

# Split dataframe by categories
start <- which(data[[1]] %in% categories)
end <- c(start[-1] - 1, nrow(data))
discard <- c()

data_frame_list <- list()
for (i in 1:length(categories)) {
  print(i)
   # select the corresponding rows
   result <- data[(start[i]+1):end[i],-1]
   # remove empty rows
   result <- result[rowSums(is.na(result)) != ncol(result),]
   # add colum saying the category
   if(nrow(result) == 0) {
     discard <- c(discard, categories[i])
   } else {
     
   result <- cbind(result, category = categories[i])
   
   # add it to the list
   data_frame_list[[categories[i]]] <- result
   }
  }

tidy_data <- do.call(rbind, data_frame_list)
rownames(tidy_data) <- NULL



rs_code.fl <- paste(data.dir, "rs_code_fitness.txt", sep = "/")
snps <- read.table(file = rs_code.fl, stringsAsFactors = F)[,1]

# remove snps that are not in the sql database
# why do we have to do this? Because some of the snps in the excel are not used anymore
tidy_data <- tidy_data[(tidy_data$`exm-number` %in% snps),]

# remove snps without score annotation
tidy_data <- tidy_data[!(is.na(tidy_data[["maj/maj"]]) |
                        is.na(tidy_data[["maj/min"]]) |
                        is.na(tidy_data[["min/min"]])),]

snps_categories <- dlply(tidy_data, .(category), .fun = function(x) x[["exm-number"]])



myregex <- "(.)?([ACGT])(.)?([ACGT])(.*)?"
tidy_data$allele_1 <- gsub(pattern = myregex, "\\2", x = tidy_data$chip.genotype)
tidy_data$allele_2 <- gsub(pattern = myregex, "\\4", x = tidy_data$chip.genotype)
tidy_data <- select(tidy_data, -chip.genotype)


write(x = tidy_data$`exm-number`,
      file = "scoring_markers_illumina.txt",
      ncolumns = 1)

snps[!snps %in% tidy_data[["exm-number"]]]

# Filter tidy_data and the customer data so that only snps available in the SQL database
# the tsv file and the excel file are used
system("grep -Fwf scoring_markers_illumina.txt data_purged.tsv > the_subset.tsv")

the_subset <- read.table(file = "data.tsv", sep = "\t",
                 na.strings =  "--", stringsAsFactors = F, check.names = F,
                 header = T)

markers <- read.table("all_marker_ids", header = F, stringsAsFactors = F)[,1]

sample_size <- ncol(the_subset)

# The 9 extra rows in tidy_data compared to the subset account to 9 snps used in two different categories
tidy_data <- tidy_data[!duplicated(tidy_data[["exm-number"]]),]
tidy_data <- tidy_data[tidy_data[["exm-number"]] %in% markers,]
the_subset <- the_subset[markers %in% tidy_data[["exm-number"]],]
markers <- markers[markers %in% tidy_data[["exm-number"]]]

# Order data frame with same order as in the snps data file
# tidy_data <- tidy_data[match(tidy_data[["exm-number"]], rownames(the_subset)),]


# customer_alleles es una lista en la que cada elemento es una tabla con cuentas de genotipos
customer_alleles <- apply(X = the_subset, MARGIN = 1, FUN = function(x) na.omit(x) %>% table)
#customer_alleles <- apply(X = the_subset, MARGIN = 1, FUN = function(x) x)

#scores <- rep(x = 0, nrow(tidy_data))
scores <- matrix(data = rep(0, nrow(tidy_data) * sample_size),
                 nrow = sample_size)
colnames(scores) <- tidy_data[["exm-number"]]


for(i in 1:nrow(tidy_data)) {
  current_row <- tidy_data[i,]
  current_exm <- current_row[["exm-number"]]
  scoring <- c(current_row[["maj/maj"]], current_row[["maj/min"]], current_row[["min/min"]]) %>% as.numeric
  genotypes <- the_subset[markers == current_exm, ] %>% as.character
  
  
  major <- paste(current_row$allele_1, current_row$allele_1, sep = "")
  minor <- paste(current_row$allele_2, current_row$allele_2, sep = "")
  hetero1 <- paste(current_row$allele_1, current_row$allele_2, sep = "")
  hetero2 <- paste(current_row$allele_2, current_row$allele_1, sep = "")
  
  
  
  major_id <- which(names(table(genotypes)) == major)
  minor_id <- which(names(table(genotypes)) == minor)
  hetero_id <- which(names(table(genotypes)) %in% c(hetero1, hetero2))
  # check we are matching the right genotypes
  if(length(major_id) + length(minor_id) == 0) {
    print(i)
    #print(scores)
    scores[i] <- NA
    } else {
      print(i)
      genotypes[genotypes != major & genotypes != minor] <- "hetero"
      
      names(scoring) <- c(major, "hetero", minor)
      
      current_scores <- scoring[genotypes]
      names(current_scores) <- NULL
      
      scores[,i] <- current_scores
      # 
      # total_score <- sum(genotypes[major_id] * scoring[1],
      #                    genotypes[hetero_id] * scoring[2],
      #                    genotypes[minor_id] * scoring[3])
      # mean_score <- total_score / sample_size
      # scores[i] <- mean_score
      # names(scores)[i] <- current_exm
    }
  }

#the_subset[rownames(the_subset) == "exm-rs1366594",]

#tidy_data[which(is.na(scores)),] %>% View

scores.fl <- paste(data.dir, "scores_99.txt", sep = "/")
tidy_scores.fl <- paste(data.dir, "tidy_data_scores.txt", sep = "/")
the_subset.fl <- paste(data.dir, "the_subset.rds", sep = "/")
snps_categories.fl <- paste(data.dir, "snps_categories.rds", sep = "/")

write.table(scores, scores.fl, quote = F, col.names = T)
write.table(tidy_data, tidy_scores.fl, quote = F, sep = "\t", na = "NA", row.names = F)
saveRDS(object = the_subset, file = the_subset.fl)
saveRDS(object = snps_categories, file = snps_categories.fl)

ambiguous_snps <- tidy_data %>%
  filter(allele_1 == "C" & allele_2 == "G" |
         allele_1 == "G" & allele_2 == "C" |
         allele_1 == "A" & allele_2 == "T" |
         allele_1 == "T" & allele_2 == "A") %>%
  .[["exm.number"]]