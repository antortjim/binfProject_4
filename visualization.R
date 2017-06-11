rm(list = ls())
library("ggplot2")
theme_set(theme_bw())
library("ggrepel")
library("dplyr")
library("GGally")
data.dir <- commandArgs(trailingOnly = T)
setwd(data.dir)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
athgene <- paste("I", 1:528, sep = "")

fl <- paste(data.dir, "maf0.05.eigenvec", sep = "/")
df <- read.table(file = fl, header = F)
NF <- system(paste("awk '{print NF; exit}' ", fl, sep = ""), intern = T) %>% as.numeric
if(NF < 202) {
  colnames(df) <- c("ind", "X", paste("C", 1:(NF-2), sep = ""))
  interval <- paste("C1:C", NF-2, sep = "")
  NF <- NF - 2
  } else {
  colnames(df) <- c("ind", "X", paste("C", 1:200, sep = ""))
  NF <- 200
  interval <- paste("C1:C", NF, sep = "")
}

density_computer <- function(df, col1, col2) {
  result <- fields::interp.surface(MASS::kde2d(df[[col1]], df[[col2]]), df[,c(col1, col2)])
}

dot_density <- density_computer(df, "C1", "C2")
df$density <- dot_density
#df$density <- 0.2

fl <- paste(data.dir, "maf0.05.eigenval", sep = "/")
eigenval <- read.table(file = fl, header = F)$V1
individuals.fl <- paste(data.dir, "1000_genomes/individuals.txt", sep = "/")
people <- read.table(file = individuals.fl,
                     col.names = c("ind", "pop", "super", "gender", "batch"),
                     fill = T)
people$batch <- factor(people$batch)
df <- merge(df, people, by = "ind")

components <- select_(df, interval)
groups <- select_(df, paste("-(", interval, ")", sep = ""))
# p <- scatmat(data.frame(components[,1:5], groups), color = "super") +
#   scale_colour_manual(values=cbPalette)
# ggsave(p, dpi = 360, filename = "scatmat.png", height = 20, width = 30)

ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
       mapping = aes(x = C1, y = C2, col = super)) +
  geom_point() +  scale_colour_manual(values=cbPalette) +
  coord_fixed()


temp <- filter(df, pop == "AthGene")
temp$ind <- factor(temp$ind %>% as.character, levels = athgene)
temp <- arrange(temp, ind)

query <- select(temp, C1, C2, C3, C4) %>% as.matrix
ref <- select(filter(df, pop != "AthGene"), C1, C2, C3, C4) %>% as.matrix

remove_fam <- read.table(file = "maf0.05.remove.fam")[,1]
athgene <- athgene[!(athgene %in% remove_fam)]
super_pop <- filter(df, pop != "AthGene") %>% .$super

#k_nearest_neighbors <-  5
#distances <- apply(X = query, MARGIN = 1, FUN = function(x) sqrt(rowSums((sweep(ref, 2, x))^2)))
#winners <- distances %>% apply(MARGIN = 2,
#                                      FUN = function(x) super_pop[match(sort(x), x)[1:k_nearest_neighbors]] %>%
#                                        table %>% which.max
#                                        )
#winners <- levels(super_pop)[winners]
#
#ancestry <- data.frame(
#  sample = athgene,
#  super = winners)
#
#write.table(ancestry, file = "pca_ancestry.txt", quote = F, row.names = F, col.names = F)
  
  
ggsave("../plots/c1c2.png", limitsize = F, dpi = 720)
ggsave("../plots/c1c2.svg", limitsize = F, dpi = 720)

ggplot(data = filter(df, pop == "AthGene"),
       mapping = aes(x = C1, y = C2, col = batch)) +
  geom_text(mapping = aes(label = batch), size = 3)

ggsave("../plots/batches.png")

# scatmat(filter(data.frame(components[,1:5], groups), pop == "AthGene"),
#         color = "batch")

#filter(df, pop == "AthGene") %>% select(C1, C2, batch) %>% group_by(batch) %>% summarise(n())


vals_df <- data.frame(value = eigenval, component = 1:NF)

ggplot(data = vals_df,
       mapping = aes(x = component, y = cumsum(value) / sum(value))) +
  geom_bar(stat = "identity", width = 1)

ggplot(data = vals_df,
       mapping = aes(x = component, y = cumsum(value)/sum(value))) +
  geom_line()


#cumsum(eigenvals) / sum(eigenvals)

genome.fl <- paste(data.dir, "maf0.05.genome", sep = "/")
genome <- read.table(file = fl, header = T)

#genome %>% filter(Z2 == 1) %>% View
ggplot(data = filter(genome, PI_HAT > 0.2),
       mapping = aes(x = Z1, y = Z2)) + geom_point() +
  geom_text_repel(aes(label = paste(FID1, FID2))) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))
