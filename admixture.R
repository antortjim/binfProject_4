#data.dir <- commandArgs(trailingOnly = T)[[1]]
data.dir <- "~/MEGA/data"
setwd(data.dir)

library("waffle")
library("ggplot2")
theme_set(theme_bw())
library("scales")
library("plyr")
library("dplyr")
library("reshape2")

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#athgene <- paste("I", 1:528, sep = "")


super_pop_name <- c("AFR", "AMR", "AthGene", "EAS", "EUR", "SAS", "OTHER")
super_pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")

sample.fl <- paste(data.dir, "representative_subset.fam", sep = "/")
samples <- read.table(sample.fl,
                      col.names = c("sample", "sample", "father", "mother", "gender", "phenotype"))

panel <- read.table(file = "panel", header = T)
athgene <- filter(samples, sample.1 == 1) %>% .$sample
panel <- rbind(panel,
               data.frame(sample = athgene,
                          pop = "AthGene",
                          super_pop = "AthGene",
                          gender = NA))


myPanel <- full_join(samples, select(panel, sample, pop, super_pop), by = "sample") %>%
  select(sample, pop, super_pop)



tbl <- cbind(sample = samples$sample, read.table("representative_subset.5.Q"))
tbl <- merge(tbl, myPanel, by = "sample")



mean_component <- ddply(tbl, .(super_pop), .fun = function(x) {
  select(x, V1:V5) %>% colMeans()
})

component_names <- mean_component$super_pop[mean_component %>% select(-super_pop) %>%
                                              apply(MARGIN = 2, FUN = which.max) %>% unlist] %>%
  as.character
colnames(tbl)[colnames(tbl) %in% paste("V", 1:5, sep = "")] <- component_names


plot_data <- melt(as.data.frame(mean_component), id.vars = "super_pop")

plot_data$variable <- factor(x = as.character(plot_data$variable),
                             levels = levels(plot_data$variable)[match(component_names, super_pops)])


ggplot(data = plot_data,
       mapping = aes(x = super_pop, y = value, fill = variable,
                     group = variable)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = F) +
  scale_fill_manual(values = cbPalette[sort(levels(tbl$super_pop)) != "AthGene"]) +
  labs(x = "Population group", y = "Component fraction")

ggsave("../plots/component_fraction.png")

component_names <- levels(panel$super_pop)[levels(panel$super_pop) != "AthGene"]

data <- tbl %>%
  #filter(super_pop == "AthGene") %>% 
  arrange(EAS > 0.5, AMR > 0.5, AFR > 0.5, EUR > 0.5, SAS > 0.5)

data_melt <- reshape2::melt(data,
             id.vars = c("sample", "super_pop", "pop"), 
             variable.name = "admixture") %>%
  arrange(sample, admixture)


ancestry <- ddply(.data = data_melt,
                  .variables =  .(sample),
                  .fun = function(x) levels(data_melt$admixture)[which.max(x$value)])
colnames(ancestry)[2] <- "admixture"

ancestry <- merge(ancestry, select(tbl, sample, super_pop, pop), by = "sample")
filter(ancestry, super_pop != "AthGene") %>% .[,2:3] %>%
  apply(MARGIN = 1, FUN = function(x) {x[1] == x[2]}) %>%
  table


#ancestry$sample <- factor(as.character(ancestry$sample), levels = athgene)

admixture.fl <- paste(data.dir, "admixture_ancestry.txt", sep = "/")
write.table(ancestry, file = admixture.fl,
            quote = F, row.names = F, col.names = F)

waffle.fl <- paste(data.dir, "waffle_data.txt", sep = "/")
write.table(x = data_melt, file = waffle.fl, quote = F, row.names = F)
data2 <- ddply(data_melt, .(sample), .fun = function(x) arrange(x, -value, sample))
  

data2 <- merge(data2, select(ancestry, sample, result = admixture), by = "sample")

data2 <- arrange(data2, super_pop)

A <- filter(data2, super_pop != "AthGene")
B <- filter(data2, super_pop == "AthGene")
B <- arrange(B, result)
data2 <- rbind(A, B)


data2$thous_genom <- factor(ifelse(data2$super_pop == "AthGene", "Company", "1000 Genomes"))


#data2 %>% filter(admixture == "AMR") %>% .$value %>% plot
data2$sample <- factor(as.character(data2$sample), levels = unique(as.character(data2$sample)))

temp <- ddply(arrange(data2, super_pop, sample, admixture), .(sample),
              function(x) levels(data2$admixture)[which.max(x$value)])

p <- ggplot(data = data2,
        mapping = aes(x = sample, y = value,
                      fill = admixture)) +

  geom_bar(stat = "identity",
           colour = "transparent", width = 1, show.legend = F) +
  facet_grid(. ~ thous_genom, scales = "free", space = "free") +
  scale_fill_manual(name = "Ancestry",
                    values = cbPalette[match(levels(data2$admixture),
                                             super_pop_name)]) +
  theme_void()


ggsave(plot = p, filename = "../plots/admixture_plot.pdf",
       width = 1000, height = 300, limitsize = F)
ggsave(plot = p, filename = "../plots/admixture_plot.png",
       width = 100, height = 30, limitsize = F, dpi = 72)

  #scale_y_continuous(expand = c(0, 0))
  
# scaler <- 1000
# 
# individuals <- arrange(individuals, super_pop)

# group_by(filter(ancestry, super_pop != "AthGene"), admixture) %>%
#   dplyr::summarise(count = n())
# 
# x <- group_by(filter(ancestry, super_pop != "AthGene"), super_pop) %>%
#   dplyr::summarise(count = n())
# 
# xmin <- c(1, cumsum(x$count)[-nrow(x)])
# xmax <- cumsum(x$count)
# xmin <- scaler * xmin / nrow(individuals)
# xmax <- scaler * xmax / nrow(individuals)
# 
# grob_data <- data.frame(
#   label = levels(ancestry$super_pop)[levels(ancestry$super_pop) != "AthGene"],
#   ymin = -10, ymax = -10,
#   xmin = xmin, xmax = xmax,
#   thous_genom = factor("1000 Genomes", levels = c("1000 Genomes", "Company"))
# )
# 
# q <- p
# for (i in 1:nrow(grob_data)) {
#   # q <- q + annotation_custom(grob=textGrob(label = grob_data$label[i],
#   #                                          hjust = 0,
#   #                                          gp = gpar(fontsize = 8)),
#   #                            xmin = grob_data$xmin[i], xmax = grob_data$xmax[i],
#   #                            ymin = -Inf, ymax = 0) + facet_grid(thous_genom)
#   q <- q + geom_text(label = grob_data$label[i],
#                      x = mean(c(grob_data$xmin[i],
#                                 grob_data$xmax[i])),
#                      y = -Inf, vjust = -1) +
#     facet_grid(. ~ grob_data$thous_genom[i])
# }
# 
# q

# 
# for (i in 1:nrow(grob_data)) {
# 
#   q <- q + geom_text(x = mean(c(grob_data$xmin[i],
#                                 grob_data$xmax[i])),
#                      label = grob_data$label[i],
#                      y = -Inf)
# }



# Waffle
myIndividuals <- c("I370", "I400", "I5", "I110", "I215")
myPops <- c("EAS", "AMR", "EUR", "AFR", "SAS")
for (i in 1:length(myIndividuals)) {
  myIndividual <- myIndividuals[i]
  
  myData <- filter(data_melt, sample == myIndividual, value != 0)
  #myData <- rbind(data, data.frame(sample = myIndividual, super_pop = "OTHER", value = 1 - sum(myData$value)))
  
  myData$value <- round(myData$value * 100, digits = 0)
  myData <- select(myData, name = admixture, vals = value)
  myData$name <- myData$name %>% as.character %>% as.factor
  #myData <- filter(myData, vals > 0)
  
  waffle(parts = myData, rows = 20, equal = F,
         colors = cbPalette[match(myData$name, super_pop_name)],
         legend_pos = "none") +
    theme(plot.title = element_text(size = 60, face = "bold")) +
    #ggtitle(paste(myIndividual, "(", myPops[i], ")", sep = ""))
    ggtitle(myIndividual)
  fl <- paste("../plots/", myPops[i], ".png", sep = "")
  ggsave(fl)
}

