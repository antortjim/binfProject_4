library("ggplot2")
library("dplyr")
library("waffle")
library("shiny")
theme_set(theme_bw())
#data.dir <- commandArgs(trailingOnly = T)[[1]]
data.dir <- "~/MEGA/AthGene/data"
setwd(data.dir)

density_computer <- function(df, col1, col2) {
  result <- fields::interp.surface(MASS::kde2d(df[[col1]], df[[col2]]), df[,c(col1, col2)])
}

prepare_pca <- function() {
  #data.dir <- commandArgs(trailingOnly = T)
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
  
  
  
  dot_density <- density_computer(df, "C1", "C2")
  df$density <- dot_density
  
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
  return(list(df, vals_df))
}

plot_pca <- function(input, df, vals_df) {
  #ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
  #df$alpha <- ifelse(df$super == "AthGene", 0.8, 0.2)
  plot2 <- ggplot(data = df,
                  mapping = aes(x = C1, y = C2,
                                col = super, alpha = 1/density)) +
    geom_point() +  scale_colour_manual(values = cbPalette,
                                        name = "Population") +
    guides(alpha = F) +
    labs(x = paste("PC1 ", round(vals_df$value[1] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""),
         y = paste("PC2 ", round(vals_df$value[2] / sum(vals_df$value), digits = 2) * 100, "% var", sep = "")) +
    geom_label(data = filter(df, sample == input$myIndividual),
               mapping = aes(x = C1, y = C2, label = sample), show.legend = F) +
    scale_alpha(range = c(.4, 0.75))
  return(plot2)
}

filter_data <- function(input) {
  waffle.fl <- paste(data.dir, "waffle_data.txt", sep = "/")
  
  data <- read.table(file = waffle.fl, header = T)
  print(data)
  myData <- filter(data, sample == input$myIndividual, value != 0)
  
  myData$value <- round(myData$value * 100, digits = 0)
  myData <- select(myData, name = admixture, vals = value)
  myData$name <- myData$name %>% as.character %>% as.factor
  return(myData)
}




ui <- fluidPage(
  titlePanel("Principal Component Analysis vs. Admixtrue"),
  
  fluidRow(
    splitLayout(cellWidths = c("30%", "70%"),
                plotOutput("plot1"),
                plotOutput("plot2"))
  ),
  
    wellPanel(
    textInput(inputId = "myIndividual",
              label = "Select an individual",
              value = "I1"),
    helpText("Enter individual id from I1 to I528.\n I370 is an East Asian example.")
    )
)

server <- function(input, output) {
  
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442","#0072B2", "#D55E00", "#CC79A7")
    result <- prepare_pca()
    df <- result[[1]]
    vals_df <- result[[2]]
    
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
    
    observe({
      output$plot2 <- renderPlot({
        plot_pca(input, df, vals_df)
      })
    })
}

shinyApp(ui = ui, server = server)