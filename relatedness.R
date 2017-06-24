library("ggvis")
library("ggplot2")
theme_set(theme_bw())
library("dplyr")
data.dir <- "~/AthGene/data"
setwd(data.dir)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

individuals.fl <- paste(data.dir, "individuals_athgene.txt", sep = "/")
people <- read.table(file = individuals.fl,
                     col.names = c("sample", "pop", "super", "gender", "batch"),
                     fill = T, header = T)
people$batch <- factor(people$batch)
#

genome.fl <- paste(data.dir, "maf0.05.genome", sep = "/")

genome <- read.table(file = genome.fl, header = T)

genome <- merge(genome, select(people, FID1 = sample, super), by = "FID1")
filter(genome, PI_HAT > 0.1) %>% nrow


bias <- filter(genome, PI_HAT > 0.1) %>%
  group_by(super) %>% summarize(count = n())

xtable::xtable(x = bias)
ggplot(data = bias,
       mapping = aes(x = super, y = count/1e3, fill = super)) +
  geom_bar(stat = "identity") +
  ylab("kCount") +
  scale_fill_manual(values = cbPalette) +
  guides(fill = guide_legend(title = "Population"))


# z1 <- ggplot(data = filter(genome, PI_HAT > 0.1 & Z1 > 0.005),
#        mapping = aes(x = Z1, fill = super)) +
#   geom_histogram(position = "dodge", aes(y = (..count..) / 1e3)) +
#   scale_fill_manual(values = cbPalette) +
#   labs(y = "kCounts") +
#   scale_x_continuous(limits = c(0, 0.20)) +
#   guides(fill = guide_legend(title = "Population"))
# 
# z2 <- ggplot(data = filter(genome, PI_HAT > 0.1 & Z1 > 0.005),
#        mapping = aes(x = Z2, fill = super)) +
#   geom_histogram(position = "dodge", aes(y = (..count..) / 1e3)) +
#   scale_fill_manual(values = cbPalette) +
#   labs(y = "kCounts") +
#   scale_x_continuous(limits = c(0, 0.20)) +
#   guides(fill = guide_legend(title = "Population"))


pi_hat <- ggplot(data = filter(genome, PI_HAT > 0),
             mapping = aes(x = PI_HAT, fill = super)) +
  geom_histogram(position = "dodge", aes(y = (..count..) / 1e3)) +
  scale_fill_manual(values = cbPalette) +
  labs(y = "kCounts") + 
  scale_x_continuous(limits = c(0, 0.2)) +
  guides(fill = guide_legend(title = "Population", nrow = 1)) +
  xlab(expression(hat(pi))) +
  theme(legend.position = "top", text = element_text(size = 20))

ggsave(filename = "../plots/pi_hat.png", plot = pi_hat,
       height = 5, width = 12)

# ggplot(data = genome,
#        mapping = aes(x = super, y = PI_HAT, fill = super)) +
#   geom_boxplot() +
#   scale_fill_manual(values = cbPalette) +
#   labs(y = "kCounts") + 
#   guides(fill = guide_legend(title = "Population"))


genome <- filter(genome, PI_HAT > 0.1)

genome$density <- fields::interp.surface(
  MASS::kde2d(genome$Z1, genome$Z2), genome[,c("Z1", "Z2")])

ibd <- ggplot(data = genome,
       mapping = aes(x = Z1, y = Z2, col = super, alpha = 1/density)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = cbPalette) +
  guides(col = guide_legend(title="Population", nrow = 1)) +
  scale_alpha(range = c(.25, .6)) +
  guides(alpha = FALSE) +
  theme(legend.position = "top", text = element_text(size = 20))
  
ibd
ggsave(filename = "../plots/IBD.png", plot = ibd)

#
genome.fl <- paste(data.dir, "athgene_relatives.genome", sep = "/")

genome <- read.table(file = genome.fl, header = T, stringsAsFactors = F)

genome <- merge(genome, select(people, FID1 = sample, super), by = "FID1")
filter(genome, PI_HAT > 0.1) %>% nrow
genome %>% nrow

#
#genome %>% filter(Z2 == 1) %>% View
visual_data <- filter(genome, PI_HAT > 0)





relatives <- ggplot(data = genome,
       mapping = aes(x = Z1, y = Z2, alpha = PI_HAT)) +
  geom_point(show.legend = F) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = cbPalette) +
  annotate("label", x = c(0.1, .95, 0.5, 0.45, 0.25, 0.05),
           y = c(0.95, 0.10, .3, .05, 0.05, 0.05),
           label = c("MZ / Dup", "PO", "FS", "HS", "C1", "U"),
           col = RColorBrewer::brewer.pal(n = 6, name = "Set1")) +
  theme(text = element_text(size = 20))

ggsave(filename = "../plots/relatedness.png", plot = relatives,
       height = 5, width = 10)



#
# genome.fl <- paste(data.dir, "athgene_relatives.genome", sep = "/")
# 
# genome <- read.table(file = genome.fl, header = T, stringsAsFactors = F)
# 
# genome <- merge(genome, select(people, FID1 = sample, super), by = "FID1")
# 
# 
# visual_data <- filter(genome, PI_HAT > 0.1)
# visual_data$pair <- 
#   paste(gsub(pattern = "(I)([0-9]*)", "\\2", x = visual_data$FID1),
#         "+",
#         gsub(pattern = "(I)([0-9]*)", x = visual_data$FID2, "\\2"),
#         sep = "")
# 
# 
# ggvis(data = visual_data, 
#       ~Z1, ~Z2,
#       opacity:= ~PI_HAT,
#       key:= ~pair
# ) %>%
#   layer_points() %>%
#   add_tooltip(html = function(data) paste("Customers: ", data$pair),
#               on = "hover")
