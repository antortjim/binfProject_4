library("shiny")
library("ggplot2")
library("reshape")
library("dplyr")
data.dir <- "."

return_data <- function(input) {
  fl <- paste(data.dir, "tidy_data_scores.txt", sep = "/")
  tidy_data <- read.table(file = fl,
                          header = T, sep = "\t", na.strings = "NA", fill = T, quote = "\"" ,
                          stringsAsFactors = F)
  
  fl_subset <- paste(data.dir, "the_subset.rds", sep = "/")
  the_subset <- readRDS(file = fl_subset)
  
  
  # Focus on VO2max
  vo2max_snps <- c("rs11549465", "rs8192678", "rs1870377", "rs11091046")
  vo2max <- tidy_data[tidy_data$SNP %in% vo2max_snps,] 
  vo2max <- select(vo2max, gene, snp = SNP, exm.number)
  vo2max$gain <- c("C", "A", "A", "C")
  vo2max$loss <- c("T", "G", "T", "A")
  vo2max$homo_gain <- c(2, 2.59 - 1.33, 57.9 - 56.4, 62.3 - 57.4)
  
  data <- the_subset[match(vo2max$exm.number, rownames(the_subset)), ] %>% as.matrix
  
  gain_alleles <- 1:4 %>% lapply(function(x) gsub(pattern = vo2max$loss[x],
                                                  x = data[x,], replacement = "\\1")) %>%
    lapply(function(x) nchar(x))
  
  gain_alleles_count <- do.call(what = rbind, args = gain_alleles) %>% t %>% rowSums()
  
  data <- as.data.frame(cbind(sample = 1:528, gain = gain_alleles_count))
  
  data$fill <- factor(ifelse(data$gain == filter(data, sample == input$mySample)$gain, 1, 0))
  return(data)
}

# data2 <- 1:4 %>%
#   lapply(function(x) gain_alleles[[x]] * vo2max$homo_gain[x]) %>%
#   do.call(what = rbind) %>% colSums()


ui <- fluidPage(
  
  fluidRow(
    column(3),
    column(6,
  plotOutput('plot1')),
  column(3),
  
  fluidRow(
    column(4),
    column(4,
    textInput("mySample", "Select a customer", 1)),
    column(4)
    )
  ) 
)


server <- function(input, output) {
  
    
    mySample <- reactive({
      mySample <- input$mySample
    })
    
    
    data <- reactive({
      return_data(input)
    })
    
    observe({
    output$plot1 <- renderPlot({
      ggplot(data = data(),
           mapping = aes(x = gain, fill = fill)) +
      geom_bar() +
      theme_bw() +
      expand_limits(x = c(0, 8)) +
      scale_y_continuous(limits = c(0, max(table(data()$gain)) + 10),
                         expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 8 , 1)) +
      scale_fill_manual(values = c("#606060", "#00FF00")) +
      guides(fill = F) +
      ggtitle(label = "VO2max", subtitle = "4 SNPs with a max gain of 20 ml min-1 kg-1") +
      labs(x = "Gain alleles", y = "# Customers") +
        theme(text = element_text(size = 20))
    })
    })
}

shinyApp(ui = ui, server = server)




# ranks <- percent_rank(x = arrange(data, gain) %>% .[,"gain"])
# names(ranks) <- dense_rank(x = arrange(data, gain) %>% .[,"gain"])              
# 
# better_than_you <- sum(names(ranks) > filter(data, sample == mySample) %>%
#                          .$gain) / length(ranks) * 100
# xintercept <- 8 * (100 - better_than_you) / 100
# 
# 
# polygon_data <- data.frame(id = factor(c(rep(0, 4), rep(1, 4))),
#                                   x = c(xintercept, 0, 0, xintercept,
#                                         8.5, xintercept, xintercept, 8.5),
#                                   y = c(0, 0, max(table(data$gain)) + 10, max(table(data$gain)) + 10,
#                                         0, 0, max(table(data$gain)) + 10, max(table(data$gain)) + 10))
# 
# p <- p + geom_vline(xintercept = xintercept) 
# p


#p <- geom_polygon(data = polygon_data, mapping = aes(group = id, x = x,
#                                                  y = y, fill = id), alpha = 0.1)
  

##

# Vo2 max data
# rs8192678 (ml min-1 kg lbm-1  ) lbm = lean body mass
# PGC1 alpha (PGC1A, ppar alpha)
# G/A. G major MINUS STRAND
# G/G (Gly/Gly) -> 1.33 +- 0.79 
# X/A (X/Ser) -> 2.59 +- 0.64 

# rs1870377 ml/min/kg
# VEGFR2 vascular endothelium growth factor receptor
# T/A. T major  PLUS STRAND. Gln -> His
# T/T (His/His) -> 56.4 +- 1.2
# T/A and AA -> 57.9 +- 1.5

# rs11091046
# AGTR2
# C allele -> 62.3 ml min-1 kg-1
# A allele -> 57.4 ml min-1 kg-1

# rs11549465
# CC: 2
# CT: 0
# TT: 0


