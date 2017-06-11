library("waffle")
library("ggplot2")
library("scales")
library("plyr")
library("dplyr")
library("reshape2")
data.dir <- arguments[[1]]


athgene <- paste("I", 1:528, sep = "")


panel <- read.table(file = "panel", header = T)
panel <- rbind(panel,
               data.frame(sample = athgene,
                 pop = "AthGene",
                 super_pop = "AthGene",
                 gender = NA))

fl <- paste(data.dir, "merged_clean.fam", sep = "/")
individuals <- read.table(file = fl,
                          col.names = c("sample", "family", "father", "mother", "gender", "pheno"))

individuals <- merge(individuals, select(panel, sample, super_pop),
                     by = "sample") %>%
  select(sample, super_pop)



fl <- paste(data.dir, "merged_clean.5.Q", sep = "/")
tbl <- read.table(fl)
tbl <- cbind(individuals, tbl)


mean_component <- ddply(tbl, .(super_pop), .fun = function(x) {
  select(x, V1:V5) %>% colMeans()
  }
)

component_names <- mean_component$super_pop[mean_component %>% select(-super_pop) %>%
                                        apply(MARGIN = 2, FUN = which.max) %>% unlist] %>%
  as.character


data <- tbl %>% filter(super_pop == "AthGene") %>% 
  select(V1:V5) %>%
  arrange(V1 > 0.5, V2 > 0.5, V3 > 0.5, V4 > 0.5, V5 > 0.5)
colnames(data) <- component_names

super_pop_name <- c("AFR", "AMR", "AthGene", "EAS", "EUR", "SAS", "OTHER")
data <- cbind(data, AthGene = 0)
data <- select(data, AFR, AMR, AthGene, EAS, EUR, SAS) %>%
  arrange(AFR > 0.5, AMR > 0.5, EAS > 0.5 ,EUR > 0.5, SAS > 0.5)


data <- melt(cbind(sample = athgene, data),
             id.vars = "sample", 
             variable_name = "super_pop")

ancestry <- ddply(data, .(sample), function(x) super_pop_name[which.max(x$value)])
ancestry$sample <- factor(as.character(ancestry$sample), levels = athgene)

write.table(arrange(ancestry, sample), file = "../admixture_ancestry.txt", quote = F, row.names = F, col.names = F)

#group_by(ancestry, V1) %>% dplyr::summarise(count = n())

data$sample <- factor(data$sample, levels = athgene)
fl <- paste(data.dir, "waffle_data", sep = "/")
write.table(x = data, file = fl, quote = F, row.names = F)
  
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = filter(data2, super_pop != "AthGene"),
        mapping = aes(x = sample, y = value, fill = super_pop, group = reorder(super_pop, value))) +
  geom_bar(stat = "identity",
           colour = "transparent", width = 1) +
  scale_fill_manual(name = "Ancestry",
                    values = cbPalette[levels(data$super_pop) != "AthGene"]) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


ggsave("../plots/admixture_plot.pdf",
       width = 1000, height = 300, limitsize = F)


