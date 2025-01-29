#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(viridis)
library(magrittr)
library(dplyr)
library(tidyverse)
library(tibble)
library(knitr)
library(uwot)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Novaseq UMAP"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput('min.count', "Minimum Count:", min=1, max=20, value=5, step=1),
      sliderInput('min.total.count', "Minimum Total:", min=1, max=20, value=10, step=1),
      sliderInput('large.n', "Large N:", min=1, max=20, value=5, step=1),
      sliderInput('min.prop', "Minimum Proportion:", min=0, max=1, value=.125, step=.0625),
      selectInput("colour.by", "Colour by:", c("Patient", "Timepoint", "Plate", "Empty Wells", "Gene Expression (Log CPM)"), selected='Plate'),
      selectizeInput('gene.id', "Gene:", NULL, selected="THY1"),
      actionButton("generate", "generate")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
          plotOutput('umapPlot', width='900px', height='700px')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  empty.wells <- read.csv('../../data/novaseq_metadata/empty_wells.csv')
  index.data <- read.csv('../../data/novaseq_metadata/index_data.csv')
  tagementation.layout <- read.csv('../../data/novaseq_metadata/tagmentation_layout.csv')
  
  # novaseq.counts <- read.table('../../results/novaseq_featureCounts_bothMatch_countReadPairs_geneName.txt', 
  #                              header=T, skip=1)[,-c(2:6)] |> 
  #   tibble::column_to_rownames('Geneid') |> as.matrix()
  # novaseq.counts <- novaseq.counts[matrixStats::rowSums2(novaseq.counts) > 1, ]
  
  # sample.df <- data.frame(id=stringr::str_extract(colnames(novaseq.counts), 'Plate_[0-9]+_[A-Z][0-9]+')) |> left_join(tagementation.layout |> dplyr::mutate(id=paste0(tag_plate, '_', location)), by='id')
  sample.df <- read.csv('sample-df.csv')
  
  design <- model.matrix(~tag_plate+timepoint, data=sample.df)
  #novaseq.dgelist <- edgeR::DGEList(novaseq.counts, samples=sample.df, remove.zeros=T) |> edgeR::calcNormFactors() |> edgeR::estimateCommonDisp(design)
  novaseq.dgelist <- readRDS('novaseq-dgelist.rds')
  novaseq.cpm <- edgeR::cpm(novaseq.dgelist, log=T) |> t()
  
  print('Startup Complete')
  
  updateSelectizeInput(session, 'gene.id', choices = rownames(novaseq.dgelist$counts) |> sort(), server = TRUE)
  
  umap.fit <- reactive({
    genes.to.keep <- edgeR::filterByExpr(novaseq.dgelist, 
                                         design=design, 
                                         min.count=input$min.count,
                                         min.total.count=input$min.total.count, 
                                         large.n=input$large.n, 
                                         min.prop=input$min.prop)
    print(paste0("Number of genes after filtering: ", sum(genes.to.keep)))
    print("Fitting UMAP...")
    .tmp <- uwot::umap2(edgeR::cpm(novaseq.dgelist, log=T)[genes.to.keep,] |> t(), seed=42) |> data.frame()
    print("Complete")
    .tmp
    }) %>% bindCache(input$min.count, input$min.total.count, input$large.n, input$min.prop)
    
  
  output$umapPlot <- renderPlot({
    point.size <- 2
    if(input$colour.by == 'Patient'){
      p <- ggplot(cbind(umap.fit(), sample.df), aes(x=X1, y=X2, colour=PID))+geom_point(size=point.size)
    }else if(input$colour.by == 'Timepoint'){
      p <- ggplot(cbind(umap.fit(), sample.df), aes(x=X1, y=X2, colour=timepoint))+geom_point(size=point.size)
    }else if(input$colour.by == 'Plate'){
      p <- ggplot(cbind(umap.fit(), sample.df), aes(x=X1, y=X2, colour=tag_plate))+geom_point(size=point.size)
    }else if(input$colour.by == 'Empty Wells'){
      p <- ggplot(cbind(umap.fit(), sample.df), aes(x=X1, y=X2, colour=empty_well))+geom_point(size=point.size)
    }else{
      gene <- input$gene.id 
      .sample.df <- sample.df %>% dplyr::mutate(!!gene := novaseq.cpm[,gene])
      sym.gene <- sym(gene)
      p <- ggplot(cbind(umap.fit(), .sample.df), aes(x=X1, y=X2, colour=!!sym.gene))+geom_point(size=point.size)
    }
    p 
  }, width='900px', height='700px') %>% bindCache(input$min.count, input$min.total.count, input$large.n, input$min.prop, input$colour.by, input$gene.id) %>% bindEvent(input$generate)
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
