  rm(list = ls())
  library("ggplot2")
  theme_set(theme_bw())
  library("ggrepel")
  library("dplyr")
  library("GGally")
  #data.dir <- commandArgs(trailingOnly = T)
  data.dir <- "~/MEGA/AthGene/data"
  setwd(data.dir)
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  athgene <- paste("I", 1:528, sep = "")
  
  fl <- paste(data.dir, "maf0.05.eigenvec", sep = "/")
  df <- read.table(file = fl, header = F)
  
  
  NF <- system(paste("awk '{print NF; exit}' ", fl, sep = ""), intern = T) %>% as.numeric
  if(NF < 202) {
    colnames(df) <- c("sample", "X", paste("C", 1:(NF-2), sep = ""))
    interval <- paste("C1:C", NF-2, sep = "")
    NF <- NF - 2
    } else {
    colnames(df) <- c("sample", "X", paste("C", 1:200, sep = ""))
    NF <- 200
    interval <- paste("C1:C", NF, sep = "")
    }
  
  fl <- paste(data.dir, "maf0.05.eigenval", sep = "/")
  eigenval <- read.table(file = fl, header = F)$V1
  vals_df <- data.frame(value = eigenval, component = 1:NF)
  
  
  density_computer <- function(df, col1, col2) {
    result <- fields::interp.surface(MASS::kde2d(df[[col1]], df[[col2]]), df[,c(col1, col2)])
  }
  
  dot_density <- density_computer(df, "C1", "C2")
  df$density <- dot_density
  #df$density <- 0.2
  
  individuals.fl <- paste(data.dir, "individuals_athgene.txt", sep = "/")
  people <- read.table(file = individuals.fl,
                       col.names = c("sample", "pop", "super", "gender", "batch"),
                       fill = T, header = T)
  people$batch <- factor(people$batch)
  df <- merge(df, people, by = "sample")
  
  components <- select_(df, interval)
  groups <- select_(df, paste("-(", interval, ")", sep = ""))
  # p <- scatmat(data.frame(components[,1:5], groups), color = "super") +
  #   scale_colour_manual(values=cbPalette)
  # ggsave(p, dpi = 360, filename = "scatmat.png", height = 20, width = 30)
  
  df$density <- fields::interp.surface(
    MASS::kde2d(df$C1, df$C2), df[,c("C1", "C2")])
  
  df[df$super == "AthGene", "density"] <- 1
  
  
  #ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
  #df$alpha <- ifelse(df$super == "AthGene", 0.8, 0.2)
  ggplot(data = df,
         mapping = aes(x = C1, y = C2,
                       col = super, alpha = 1/density)) +
    geom_point() +  scale_colour_manual(values = cbPalette,
                                        name = "Population") +
    guides(alpha = F) +
    labs(x = paste("PC1 ", round(vals_df$value[1] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""),
         y = paste("PC2 ", round(vals_df$value[2] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""))
  
  

  ggsave("../plots/c1c2.png", limitsize = F, dpi = 720)
ggsave("../plots/c1c2.svg", limitsize = F, dpi = 720)

ggplot(data = filter(df, pop == "AthGene"),
       mapping = aes(x = C1, y = C2, col = batch)) +
  geom_text(mapping = aes(label = batch), size = 3, show.legend = F)

ggsave("../plots/batches.png")


# temp <- filter(df, pop == "AthGene")
# temp$sample <- factor(temp$sample %>% as.character, levels = athgene)
# temp <- arrange(temp, sample)
# 
# query <- select(temp, C1, C2, C3, C4) %>% as.matrix
# ref <- select(filter(df, pop != "AthGene"), C1, C2, C3, C4) %>% as.matrix
# 
# remove_fam <- read.table(file = "maf0.05.remove.fam")[,1]
# athgene <- athgene[!(athgene %in% remove_fam)]
# super_pop <- filter(df, pop != "AthGene") %>% .$super

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

# scatmat(filter(data.frame(components[,1:5], groups), pop == "AthGene"),
#         color = "batch")

#filter(df, pop == "AthGene") %>% select(C1, C2, batch) %>% group_by(batch) %>% summarise(n())


ggplot(data = vals_df,
       mapping = aes(x = component, y = cumsum(value) / sum(value))) +
  geom_bar(stat = "identity", width = 1)

ggplot(data = vals_df,
       mapping = aes(x = component, y = cumsum(value)/sum(value))) +
  geom_line() + geom_point(size = 0.5) +
  labs(x = "Component", y = "Cumulative variance")

ggsave("../plots/cumulative_variance.png")


genome.fl <- paste(data.dir, "maf0.05.genome", sep = "/")

genome <- read.table(file = genome.fl, header = T)

genome <- merge(genome, select(people, FID1 = sample, super), by = "FID1")
group_by(genome, super) %>% summarise(count = n())
filter(genome, PI_HAT > 0.1) %>% nrow
genome %>% nrow


#genome %>% filter(Z2 == 1) %>% View
ggplot(data = genome,
       mapping = aes(x = Z1, y = Z2, col = super, alpha = PI_HAT)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = cbPalette) +
  guides(alpha = FALSE)

ggsave("../plots/IBD.png")
