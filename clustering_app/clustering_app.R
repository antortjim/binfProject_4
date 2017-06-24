library("ggplot2")
library("dplyr")
library("waffle")
library("shiny")
theme_set(theme_bw())
#data.dir <- commandArgs(trailingOnly = T)[[1]]
data.dir <- "."
#setwd(data.dir)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

density_computer <- function(df, col1, col2) {
  result <- fields::interp.surface(MASS::kde2d(df[[col1]], df[[col2]]), df[,c(col1, col2)])
}

prepare_pca <- function() {
  #data.dir <- commandArgs(trailingOnly = T)
  athgene <- paste("I", 1:528, sep = "")
  
  fl <- paste(data.dir, "main_two.eigenvec", sep = "/")
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
  
  
  
  df_plotted <- filter(df, pop %in%  c("ITU", "ACB", "CEU", "PEL", "AthGene", "CHS"))
  
  df_plotted$density <- fields::interp.surface(
    MASS::kde2d(df_plotted$C1, df_plotted$C2), df_plotted[,c("C1", "C2")])
  
  df_plotted[df_plotted$super == "AthGene", "density"] <- 1
  return(list(df_plotted, vals_df))
  
  
  #ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
  #df$alpha <- ifelse(df$super == "AthGene", 0.8, 0.2)
  
}

plot_pca <- function(input, df_plotted, vals_df) {
  #ggplot(data = filter(df, pop %in% c("CEU", "PEL", "CHB", "AthGene", "ITU", "ACB")),
  #df$alpha <- ifelse(df$super == "AthGene", 0.8, 0.2)
  plot2 <- ggplot(data = df_plotted,
                  mapping = aes(x = C1, y = C2,
                                col = super, alpha = 1/density)) +
    geom_point() +
    theme(text = element_text(size = 20),
          legend.position = "top") +
    scale_colour_manual(values = cbPalette) +
    guides(alpha = F, col = guide_legend(nrow = 1, title = "Population", override.aes = list(size = 10))) +
    labs(x = paste("PC1 ", round(vals_df$value[1] / sum(vals_df$value), digits = 2) * 100, "% var", sep = ""),
         y = paste("PC2 ", round(vals_df$value[2] / sum(vals_df$value), digits = 2) * 100, "% var", sep = "")) +
    scale_alpha(range = c(.4, 0.75)) +
    scale_x_continuous(breaks = seq(-0.02, 0.04, 0.02)) +
    scale_y_continuous(breaks = seq(-0.02, 0.04, 0.02)) +
   geom_label(data = filter(df_plotted, sample == input$myIndividual),
             mapping = aes(x = C1, y = C2, label = sample), show.legend = F)
  return(plot2)
}

filter_data <- function(input) {
  waffle.fl <- paste(data.dir, "waffle_data.txt", sep = "/")
  
  data <- read.table(file = waffle.fl, header = T)
  myData <- filter(data, sample == input$myIndividual, value != 0)
  
  myData$value <- round(myData$value * 100, digits = 0)
  myData <- select(myData, name = admixture, vals = value)
  myData$name <- myData$name %>% as.character %>% as.factor
  return(myData)
}


ui <- fluidPage(
  titlePanel("Admixture vs. Principal Component Analysis"),
  fluidRow(
    column(4,
           textInput(inputId = "myIndividual",
                     label = "Select an individual",
                     value = "I1"),
           helpText("Enter individual id from I1 to I528.\n I370 is an East Asian example."), # end of sideBar
           plotOutput("plot1")),
    column(8, plotOutput("plot2", height = 740, width = 740))
    )
)

    
server <- function(input, output) {
  
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442","#0072B2", "#D55E00", "#CC79A7")
    result <- prepare_pca()
    df_plotted <- result[[1]]
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
        plot_pca(input, df_plotted, vals_df)
      })
    })
}

shinyApp(ui = ui, server = server)