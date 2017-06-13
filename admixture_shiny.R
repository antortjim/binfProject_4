library("ggplot2")
library("dplyr")
library("waffle")
library("shiny")

#data.dir <- commandArgs(trailingOnly = T)[[1]]
data.dir <- "~/MEGA/AthGene/data"
setwd(data.dir)



ui <- fluidPage(
  titlePanel("Admixture result"),
  mainPanel(
    plotOutput('plot1'),
    textInput("myIndividual", "Select an individual", "I1")
  )  # end mainPanel
)

server <- function(input, output) {
  
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    filter_data <- function(input) {
    waffle.fl <- paste(data.dir, "waffle_data.txt", sep = "/")
    
    data <- read.table(file = waffle.fl, header = T)
    print(data)
    myData <- filter(data, sample == input$myIndividual, value != 0)
    
    myData$value <- round(myData$value * 100, digits = 0)
    myData <- select(myData, name = super_pop, vals = value)
    myData$name <- myData$name %>% as.character %>% as.factor
    return(myData)
  }
  myData <- reactive({
    filter_data(input)
  })
    
    observe({
      output$plot1 <- renderPlot({
        super_pop_name <- c("AFR", "AMR", "AthGene", "EAS", "EUR", "SAS", "OTHER")
        waffle(parts = myData(), rows = 10,
               colors = cbPalette[match(myData()$name, super_pop_name)]) +
          ggtitle(input$myIndividual)
      })
    })
}

shinyApp(ui = ui, server = server)
