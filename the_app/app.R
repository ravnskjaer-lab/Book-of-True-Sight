.libPaths( c( .libPaths(), "/work/R/") )

.libPaths()

library(BiocManager)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(matrixStats)
library(MatrixGenerics)
library(Biobase)
library(SummarizedExperiment)
library(DESeq2)
library(ggpubr)
library(ggsci)
library(shinythemes)
library(dplyr)
library(tidyverse)
library(shiny)



#Path to datasets and list of datasets named files
files <- dir(path = "/work/bookoftruesight/data/",pattern="*\\.rds$")    
files <- substr(files, 1,nchar(files)-4)


### USER INTERFACE ###

ui <- fluidPage(
  column(8,offset = 4, titlePanel("Book of True Sight")),
  mainPanel(
  h3("LG/KR bulk RNA-seq datasets"),
  h6("Made by FT (2022)®")),
  
  
  #Gene input box
  fluidRow(column(6,
                  selectInput(
                    "Dataset", "Choose a dataset", files,
                    selected = "WD_52weeks",
                    multiple = FALSE),
                  uiOutput("tab")
  ),
  column(6,
         textOutput("ownership"),
         textOutput("design"),
         actionButton("compute","Compute differential expression analysis")
  )
  ),
  
  #Gene input box
  fluidRow(column(6, "Variables to plot",
                  selectizeInput(
                    "gene", "Gene of choice", choices = NULL,
                    multiple = FALSE,
                    selected = NULL)
                  ,
                  
                  #Clinical variable input box 
                  selectInput(
                    "variables", "Condition", choices = NULL,
                    multiple = FALSE,
                    selected = NULL)
  ),
  
  #Contrasts to use in DESeq results table
    column(6, "Contrasts for differential expression analysis",
         uiOutput("contrast1"),
         uiOutput("contrast2"),
         uiOutput("covariate")
  )),
  
  #Plot input box
  fluidRow(column(12, tableOutput("data")),
           column(5, "Boxplot", plotOutput("boxplot")),
           column(5, "PCA", plotOutput("PCA")),
           column(12, "Confounding factors (N.B. PCs not comparable to PCA plot)", plotOutput("confounding.factors"))
             ),
  theme = shinytheme("journal"),
)

### SERVER LOGIC ####

server <- function(input, output, session) {
  
  #Dataset to use based on input dataset 
  dataset <- reactive({
    filename <- paste0("/work/bookoftruesight/data/",
                       input$Dataset,
                       ".rds")
    expression.data <- readRDS(filename)
    
    return(expression.data)
  })
  
  
  observeEvent(dataset(),{
    freezeReactiveValue(input, "gene")
    updateSelectizeInput(inputId = "gene", choices = rownames(assay(dataset())),server = T)
  })
  
  observeEvent(dataset(),{
    freezeReactiveValue(input, "variables")
    updateSelectInput(inputId = "variables", choices = colnames(colData(dataset()))[c(1:paste(length(colnames(colData(dataset())))))])
  })
  

  output$contrast1<-renderUI({
    validate(
    need(input$variables, "Have patience young one")
  )
      selectInput(inputId = "contrast1",
                  "Contrast 1",
                  choices = levels(factor(colData(dataset())[,input$variables])),
                  selected = NULL,multiple = F)})
  
  output$contrast2<-renderUI({
    validate(
      need(input$variables, "Have patience young one")
    )
    selectInput(inputId = "contrast2",
                "Contrast 2",
                choices = levels(factor(colData(dataset())[,input$variables])),
                selected = NULL,multiple = F)})  

  output$covariate<-renderUI({
    validate(
      need(input$variables, "Have patience young one")
    )
    selectInput(inputId = "covariate",
                "Covariate",
                choices = c("NONE", colnames(colData(dataset())[!colnames(colData(dataset())) %in% input$variables])),
                selected = NULL,multiple = F)
    }
    )  
  
    #Creating the dataframe to plot from based on input dataset 
  df <- reactive({
    df.1 <- data.frame(t(assay(vst(dataset(), blind = F))))
    name <- paste0(input$variables)
    df.1[[name]] <- colData(dataset())[[input$variables]]
    
    return(df.1)
  })
  
  #Creating a dataframe of PCs and metadata to test for confounding factors
  confounding.factors <- reactive({
    confounding.factors.1 <- as.data.frame(
      prcomp(assay(vst(dataset(), blind = F)))[2]$rotation)[,1:6]
    
    variables <- colnames(colData(dataset())[1:length(colnames(colData(dataset())))])
    
    pca.prcomp = as.data.frame(prcomp(assay(vst(dataset(), blind = F)))[2]$rotation)[,1:6]
    covariate_list <- list()
    
    for (i in 1:length(variables)) {
      variable.name <- paste0(variables[[i]])
      pca.prcomp[[variable.name]] <- colData(dataset())[[variable.name]]
      
      covariate_list[[variable.name]] <- summary(aov(cbind(PC1,PC2,PC3,PC4,PC5,PC6) ~ pca.prcomp[[variable.name]], data = pca.prcomp))
      covariate_list[[variable.name]] <- unlist(lapply(covariate_list[[variable.name]],function(x) x[1,5]))
      covariates.df <- data.frame(do.call(cbind, covariate_list))
      covariates.df$PC <- substr(rownames(covariates.df), nchar(rownames(covariates.df))-3,nchar(rownames(covariates.df)))
      covariates.df <- data.frame(pivot_longer(covariates.df, cols = 1:ncol(covariates.df)-1, names_to = "effect", values_to = "pval"))
    }
    return(covariates.df)
  })
  
  #Table from differential expression analysis 
  output$design <- renderText({paste("Model: ", paste("~",ifelse(input$covariate == "NONE","",paste(input$covariate,"+")),input$variables), ": Contrasts = ", input$contrast2, "vs.", input$contrast1)})
  output$ownership <- renderText({paste("Ownership of data: ",unique(mcols(dataset())$owner))})
  
  url <- a("LG/KR bulk RNA-seq data", href="https://syddanskuni.sharepoint.com.mcas.ms/:x:/r/Sites/Hepatic_fanatics/_layouts/15/Doc.aspx?sourcedoc=%7B9C879D54-4F98-4CC4-A090-50A64DB5B9CD%7D&file=LG.KR_datasets.xlsx&action=default&mobileredirect=true",
           target = "_blank")
  output$tab <- renderUI({
    tagList("Dataset information:", url)
  })
  #Box plot  
  output$boxplot <- renderPlot(width = 300, height = 400,
                               {
                                 validate(
                                   need(input$gene, 'Select at least one gene!'),
                                   need(input$variables != '', 'Please choose a variable')
                                 )
                                 ggboxplot(df(), x = input$variables, y = input$gene,
                                           fill = input$variables,
                                           add = "jitter")+
                                   theme(legend.position = "none")+xlab("")+ylab(paste(input$gene, " expression (normalised counts)"))+
                                   stat_compare_means()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                                   scale_fill_manual(values = pal_jama("default")(length(unique(colData(dataset())[[input$variables]]))))
                               })
  
  #PCA plot
  output$PCA <- renderPlot({
    validate(
      need(input$variables != '', 'Please choose a variable')
    )
    plotPCA(vst(dataset(), blind = F), intgroup = input$variables)+theme_minimal()+
      scale_color_manual(values = pal_jama("default")(length(unique(colData(dataset())[[input$variables]]))))
  })  
  
  
  #Plot of confounding factors ANOVA
  output$confounding.factors <- renderPlot({
    
    dyrup <- colorRampPalette(c('#25557B', '#F4CE80', '#C4031D'))(100) #Color scale from dyRup (Terkelsen et al.(2021))
    
    ggplot(confounding.factors(), aes(x=PC, y=effect, size = -log10(pval), fill = -log10(pval)))+ #Using log scale 
      geom_point(alpha=1,shape=21, color="black")+ 
      scale_size(range = c(4, 20), limits = c(1.30103,max(-log10(confounding.factors()$pval))))+ #Scaling the limit from 1.3 excludes p values below 0.05 
      scale_fill_gradientn(colours = dyrup, guide = guide_colorbar(frame.colour = "black", ticks = T, frame.linewidth = 1, ticks.linewidth = 1, ticks.colour = "black"))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      guides(size = "none")
    
  })
  
  table <- eventReactive(input$compute,{
    
    validate(need(input$variables != "", "Please choose a variable"))
    dds <- DESeqDataSetFromMatrix(countData = assay(dataset()),
                             colData = colData(dataset()),
                             design = formula(paste("~", if(input$covariate == "NONE"){}
                                                    else{paste(input$covariate)},"+", input$variables)))
    de <- DESeq(dds)
    results.table <- data.frame(results(de, contrast = c(input$variables, input$contrast2, input$contrast1)))
    results.table[,"pvalue"] <- sprintf('%.0e',results.table[,"pvalue"])
    results.table[,"padj"] <- sprintf('%.0e',results.table[,"padj"])
    
    return(results.table)
    }
  )
  
  
  output$data <- renderTable({
    validate(
      need(input$contrast1 != input$contrast2, "Please choose different contrasts"))
    table()[input$gene,]
    }, rownames = T)
  
  
}

shinyApp(ui, server)  


