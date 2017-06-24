---
  title: "Relatedness"
author: "Antonio Ortega Jiménez"
date: "June 17, 2017"
output: html_document
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library("ggvis")
library("ggplot2")
library("dplyr")
library("shiny")
data.dir <- "~/AthGene/data"
setwd(data.dir)

individuals.fl <- paste(data.dir, "individuals_athgene.txt", sep = "/")
people <- read.table(file = individuals.fl,
                     col.names = c("sample", "pop", "super", "gender", "batch"),
                     fill = T, header = T)
people$batch <- factor(people$batch)
#
genome.fl <- paste(data.dir, "athgene_relatives.genome", sep = "/")

genome <- read.table(file = genome.fl, header = T, stringsAsFactors = F)

genome <- merge(genome, select(people, FID1 = sample, super), by = "FID1")


visual_data <- filter(genome, PI_HAT > 0.1)

visual_data$pair <- 
  paste(gsub(pattern = "(I)([0-9]*)", "\\2", x = visual_data$FID1),
        "+",
        gsub(pattern = "(I)([0-9]*)", x = visual_data$FID2, "\\2"),
        sep = "")


ui <- fluidPage(
  mainPanel(
    plotOutput('plot1'),
                "Relatedness")
  )

# Initalize parameters of the simulation and data structures

server <- function(input, output, df = visual_data) {
  
  print(df)
  
  output$plot1 <- renderPlot({
    ggvis(data = df, 
        ~Z1, ~Z2,
        opacity:= ~PI_HAT,
        key:= ~pair
  ) %>%
    layer_points() %>%
    add_tooltip(html = function(data) paste("Customers: ", data$pair),
                on = "hover")
  })
  
}

shinyApp(ui = ui, server = server)
