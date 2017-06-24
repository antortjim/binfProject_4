rm(list = ls())
library("ggplot2")
theme_set(theme_bw())
library("ggrepel")
library("dplyr")
library("GGally")
#data.dir <- commandArgs(trailingOnly = T)
data.dir <- "~/AthGene/data"
setwd(data.dir)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
athgene <- paste("I", 1:528, sep = "")

fl <- paste(data.dir, "custom.eigenvec", sep = "/")
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

fl <- paste(data.dir, "custom.eigenval", sep = "/")
eigenval <- read.table(file = fl, header = F)$V1
vals_df <- data.frame(value = eigenval, component = 1:NF)


density_computer <- function(df, col1, col2) {
  result <- fields::interp.surface(MASS::kde2d(df[[col1]], df[[col2]]), df[,c(col1, col2)])
}

dot_density <- density_computer(df, "C1", "C2")
df$density <- dot_density
#df$density <- 0.2

df$batch <- factor(rep(1:11, each = 48))

#ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
#df$alpha <- ifelse(df$super == "AthGene", 0.8, 0.2)
ggplot(data = df,
              mapping = aes(x = C1, y = C2, col = batch)) +
  geom_text(aes(label = batch)) +
  #geom_point(aes(alpha = 1/density), show.legend = F) +
  theme(text = element_text(size = 20)) +
  labs(x = paste("PC1 ", round(vals_df$value[1] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""),
       y = paste("PC2 ", round(vals_df$value[2] / sum(vals_df$value), digits = 2) * 100, "% var", sep = "")) +
  scale_alpha(range = c(.4, 1))
  
ggsave(plot = pca, "../plots/c1c2.png",
       width = 20, height = 10, dpi = 288) 
# myIndividual <- "I370"
# 
# 
# ggplot(data = df,
#        mapping = aes(x = C1, y = C2,
#                      col = super, alpha = 1/density)) +
#   geom_point() +  scale_colour_manual(values = cbPalette,
#                                       name = "Population") +
#   guides(alpha = F) +
#   labs(x = paste("PC1 ", round(vals_df$value[1] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""),
#        y = paste("PC2 ", round(vals_df$value[2] / sum(vals_df$value), digits = 2) * 100, "% var", sep = "")) +
#   geom_label(data = filter(df, sample == myIndividual),
#              mapping = aes(x = C1, y = C2, label = sample), show.legend = F)



ggsave("../plots/c1c2.svg", limitsize = F, dpi = 720)

df_plotted$density <- fields::interp.surface(
  MASS::kde2d(df_plotted$C1, df_plotted$C2), df_plotted[,c("C1", "C2")])


ggplot(data = filter(df_plotted, pop == "AthGene", C2 < 0.02 | batch != 8),
       mapping = aes(x = C1, y = C2, col = batch)) +
  geom_text(mapping = aes(label = batch), size = 7, show.legend = F) +
  geom_text_repel(data = filter(df_plotted, pop == "AthGene", C2 > 0.02 & batch == 8),
                  aes(label = batch), size = 7, show.legend = F) +
  scale_alpha(range = c(0.25, 1))

ggsave("../plots/batches.png",
       height = 5, width = 5)


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

ggplot(data = rbind(c(0, 0), vals_df),
       mapping = aes(x = component, y = cumsum(value)/sum(value))) +
  geom_line() + geom_point(data = vals_df,
                           size = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200)) +
  labs(x = "Component", y = "Cumulative variance")

ggsave("../plots/cumulative_variance.png", height = 5, width = 10)
