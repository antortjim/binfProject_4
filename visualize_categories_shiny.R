#data.dir <- commandArgs(trailingOnly = T)[[1]]
data.dir <- "~/AthGene/data"
setwd(data.dir)

library("shiny")

scores.fl <- paste(data.dir, "scores_99.txt", sep = "/")
scores <- read.table(file = scores.fl, header = T, stringsAsFactors = F)
tidy_scores.fl <- paste(data.dir, "tidy_data_scores.txt", sep = "/")
tidy_data <- read.table(file = tidy_scores.fl,
                        header = T, sep = "\t", na.strings = "NA", fill = T, quote = "\"" ,
                        stringsAsFactors = F)
snps_categories.fl <- paste(data.dir, "snps_categories.rds", sep = "/")
snps_categories <- readRDS(file = snps_categories.fl)

the_subset.fl <- paste(data.dir, "the_subset.rds", sep = "/")
the_subset <- readRDS(file = the_subset.fl)
categories <- names(snps_categories)
category_selector <- 1:length(categories)
names(category_selector) <- categories

ui <- fluidPage(
  titlePanel("Scores distribution"),
  textOutput('number_snps'),
  mainPanel(
    plotOutput('plot1'),
    selectInput("category", "Select a category", categories,
                "VO2max")

  )  # end mainPanel
)


# Initalize parameters of the simulation and data structures

server <- function(input, output, session) {
  
  library("ggplot2")
  theme_set(theme_bw())
  library("reshape")
  library("dplyr")
  

  filter_data <- function(input) {
    
    max_scores <- apply(cbind(tidy_data[["maj.maj"]] %>%
                                as.numeric,
                              tidy_data[["maj.min"]] %>%
                                as.numeric,
                              tidy_data[["min.min"]] %>%
                                as.numeric),
                        1, max)
    
    names(max_scores) <- tidy_data[["exm.number"]]
    
    
    myCategory <- input$category
    #print(myCategory)
    mySNPs <- snps_categories[[myCategory]]
    number_snps <- length(mySNPs)
    #max_scores[names(max_scores) %in% mySNPs]
    #sum(colnames(scores) %in% make.names(mySNPs))
    max_score <- sum(max_scores[mySNPs])
    
    category_scores <- scores[,colnames(scores) %in% make.names(mySNPs)]
    #print(category_scores)
    
    
    if(length(mySNPs) == 1) {
      myData <- data.frame(ind = 1:nrow(scores),
                           score = category_scores)
    } else {
      myData <- data.frame(ind = 1:nrow(scores),
                           score = rowSums(category_scores))
    }
    
    return(list(myData, max_score, number_snps))
  }
  
  result <- reactive({
    filter_data(input)
  })
  
  myData <- reactive({
    result()[[1]]
    })
  
  max_score <- reactive({
    result()[[2]]
    })
  
  number_snps <- reactive({
    result()[[3]]
  })
  

  observe({
    output$number_snps <- renderText({
      paste("Number of SNPs: ", number_snps(), sep = "")
    })
    
    output$plot1 <- renderPlot({
      # if(myData()$score %>% is.null) {
      #   p <- ggplot()
      # } else {
        #print(myData())
        p <- ggplot(data = myData(),
                      mapping = aes(x = score)) +
          geom_histogram(position = "dodge", binwidth = .5) +
          expand_limits(x = max_score()) + 
          ggtitle(label = "",
                  subtitle = paste("Max score: ", max_score(), sep = ""))
      #
      p
    })
    
  })
  

  
  
  # 
  # ggplot(data = as.data.frame(table(myData$score)),
  #        mapping = aes(x = Var1, y = Freq, group = factor(1))) +
  #   geom_line() +
  #   labs(x = "Score", y = "Count")

}

shinyApp(ui = ui, server = server)
