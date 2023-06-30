#######################################################
## This file contains snippets of the proteomics app.##
## The server function code contains all of the      ##
## statistical analyses/processing.                  ##
## Input$x refers to user input from the UI function.##
## Parts of the code may refer to some functions     ##
## that are not provided here.                       ##
#######################################################
## The app has more functionality than included      ## 
## here, but parts of the modules (data processing,  ##
## Venn diagram, PCA, volcano plot, and gene         ##
## ontology) are present.                            ##
#######################################################
library(shiny)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(shiny)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(rgl)
library(dplyr)
library(cluster)
library(rlist)
library(tibble)
library(ggvenn)
library(upsetjs)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(BiocManager)
library(msigdbr)
library(rlang)
library(tidyr)
library(enrichplot)
library(europepmc)
library(SPIA)
library(pathview)
library(grid)
library(png)
library(DT)
library(preprocessCore)
library(rentrez)
library(M3C)
library(GOSemSim)
library(circlize)

################## UI #####################
ui <- fluidPage(
  
  titlePanel("Proteomics"),
  
  sidebarLayout(
    ############## Sidebar ##############
    sidebarPanel(
      ### Select files
      h4("Upload Files"),
      fileInput("fileone", "Expression File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      fileInput("filetwo", "Metadata File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      tags$hr(),
      
      ## Extreme values
      h4("Extreme Values"),
      radioButtons("na_processing", "Proteins with any NA values: ",
                   choices = c("Omit", "Replace with 0", "Replace with 0.0001"),
                   selected = "Omit"),
      radioButtons("zero_processing", "Proteins with any 0 values:",
                   choices = c("Keep", "Omit", "Replace with 0.0001"),
                   selected = "Keep"),
      sliderInput("filter_small", "Exclude __% of lowest expressed proteins: ",
                  0,100,0,1),
      tags$hr(),
      
      ## Protein ID
      h4("Protein Identification"),
      radioButtons("organism", "Organism: ",
                   choices = c("Human (Homo sapiens)" = "Homo sapiens", "Mouse (Mus musculus)" = "Mus musculus"),
                   selected = "Homo sapiens"),
      radioButtons("keytype", "Accession used: ",
                   choices = c("UNIPROT", "ENTREZID", "Gene Symbol" = "SYMBOL"),
                   selected = "UNIPROT"),
      radioButtons("keytype_split", "Multiple accessions combined using:",
                   choices = c(";", "_", ",", "None"),
                   selected = ";"),
      radioButtons("duplicates", "If a protein is present in multiple rows: ",
                   choices = c("Keep first protein only" = "Omit", "Keep all proteins" = "Keep"),
                   selected = "Omit"),
      tags$hr(),
      
      ## Normalization
      h4("Normalization"),
      radioButtons("centermeans", "Normalize Samples:",
                   choices = c("None" = "None", Quantile = "Quantile"),
                   selected = "None"),
      tags$hr(),
      
      ## Stats params
      h4("Statistical Parameters"),
      sliderInput("sig_cutoff", "Significance",
                  0.001, 1, 0.05, 0.001),
      sliderInput("fc_cutoff", "Fold Change Cutoff",
                  1, 10, 2, 0.05),
      radioButtons("adjustment", "Adjustment",
                   choices = c(None = "none",
                               BenjaminiHochberg = "BH",
                               Bonferroni = "bonferroni",
                               Holm = "holm"),
                   selected = "none"),
      
      actionButton("submit", "Submit"),
      actionButton("downloadoptions", "Download Report")
    ), ## close sidebarPanel
    
    ################ Main Panel ####################
    mainPanel(
      tabsetPanel(type="tab",
                  ############# Instructions #############
                  tabPanel("Instructions", 
                           tags$br(),
                           h3("Before Uploading"),
                           h5("Make sure that the sample names are spelled and ordered identically in the expression and metadata file."),
                           h5("All names must start with characters only and contain no spaces."),
                           h5("Each group must have at least two samples or replicates."),
                           h5("The first column must be a list of proteins that are identified by Gene Symbol, UNIPROT, or ENTREZID."),
                           h5("If there are multiple protein identifiers for a single protein, please note that only the first one will be used!"),
                           tags$hr(),
                           h3("Formatting"),
                           h4("Expression file"),
                           tags$img(src='expression_format.png', height=150, width=600),
                           tags$br(),
                           h4("Metadata file (no headers)"),
                           tags$img(src='metadata_format.png', height=150, width=200),
                           tags$hr(),
                           h3("Pipeline"),
                           tags$img(src='proteomics_pipeline.png', height=325, width=600)
                  ),
                  
                  ########## DATA #########
                  tabPanel("Data",
                           conditionalPanel(condition = "input.submit > 0",
                                            h3("Summary"),
                                            htmlOutput("length"),
                                            tags$hr(),
                                            
                                            h3("Data"),
                                            dataTableOutput("contents"),
                                            downloadButton("contents_download", "Download Table"),
                                            tags$hr(),
                                            
                                            h3("Statistics"),
                                            dataTableOutput("contents_stats"),
                                            downloadButton("contents_stats_download", "Download Table"),
                                            tags$hr(),
                                            
                                            h3("Grouping"),
                                            dataTableOutput("meta_check")
                           )
                  ),
                  
                  ############ EDA ############
                  tabPanel("Exploratory Data Analysis",
                           tabsetPanel(type="pills",
                                       
                                       tabPanel("Univariate",
                                                tags$br(),
                                                h5("The Univariate tab shows statistics for a specific protein across all samples."),
                                                textInput("protein_viewer_select", "Search Protein: "),
                                                actionButton("protein_viewer_run", "Run"),
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.protein_viewer_run > 0",
                                                                 h3("Density Plot"),
                                                                 plotOutput("density_plot"),
                                                                 tags$hr(),
                                                                 
                                                                 h3("Quantile-Quantile Plot"),
                                                                 plotOutput("qq_plot"),
                                                                 tags$hr(),
                                                ),
                                       ), ## close Univariate
                                       
                                       tabPanel("Venn",
                                                tags$br(),
                                                h5("The Venn diagram shows how many proteins are present in a group or intersection of groups."),
                                                checkboxGroupInput("venn_choices", "Select up to 4 groups: "),
                                                sliderInput("filter_presence", "Protein is absent in group if detected in less than __% of samples: ",
                                                            1,100,100,1),
                                                radioButtons("venn_enable_upload", "Upload a separate protein list to find overlaps?",
                                                             choices = c("Disable", "Enable"), selected = "Disable"),
                                                
                                                conditionalPanel(condition = "input.venn_enable_upload == 'Enable'",
                                                                 radioButtons("venn_file_keytype", "Accession: ",
                                                                              choices = c("Gene Symbol" = "SYMBOL", "UNIPROT" = "UNIPROT", "ENTREZID" = "ENTREZID"), selected = "SYMBOL"),
                                                                 fileInput("venn_uploadfile", "Upload File: ",
                                                                           multiple = TRUE,
                                                                           accept = c("text/csv",
                                                                                      "text/comma-separated-values,text/plain",
                                                                                      ".csv"))
                                                ),
                                                actionButton("venn_run", "Run"),
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.venn_run > 0",
                                                                 h3("Venn Diagram"),
                                                                 plotOutput("venndiagram")
                                                ),
                                       ), ## close Venn 
                                       
                                       tabPanel("PCA",
                                                tags$br(),
                                                h5("Choose the number of proteins to include in the PCA. The proteins selected have the highest coefficient of variation."),
                                                numericInput("pca_topcv", "Number of Proteins:", 100, min = 2),
                                                radioButtons("pca_samples", "Select samples:", choices = c("All", "Select"), selected = "All"),
                                                conditionalPanel(condition = "input.pca_samples == 'Select'",
                                                                 h5("Select which samples should be included in the PCA."),
                                                                 checkboxGroupInput("pca_sample_choices", "Select Samples:")
                                                ),
                                                actionButton("pca_run", "Run"),
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.pca_run > 0",
                                                                 h3("2-D Visualization"),
                                                                 plotOutput("pca"),
                                                                 tags$hr()
                                                                 )
                                                ),
                                       ), ## close PCA
                           ), ## close tabsetPanel
                  ), ## close tabPanel
                  
                  
                  ########### DIFFERENTIAL EXPR ############
                  tabPanel("Differential Expression",
                           tabsetPanel(type="pills",
                                       
                                       tabPanel("Volcano",
                                                tags$br(),
                                                h5("Compare two groups. The log fold change reflects the ratio of expression between Group 2 and Group 1."),
                                                fixedRow(
                                                  column(6, checkboxGroupInput("groups1", "Group 1 (Control): ")),
                                                  column(6, checkboxGroupInput("groups2", "Group 2 (Experimental): "))
                                                ),
                                                tags$hr(),
                                                
                                                h5("Select parameters to perform the difference of means test for each protein."),
                                                radioButtons("volcano_test_type", "Significance Test:", 
                                                             choices = c("Parametric (T-test)" = "Parametric", "Nonparametric (Wilcoxon)" = "Nonparametric"), 
                                                             selected = "Parametric"),
                                                tags$hr(),
                                                
                                                h5("Select proteins to view in the plot and table. Select 'All' to see all proteins, 'Significant' to see only significant proteins, and 'Protein List' to see a specific list of proteins."),
                                                radioButtons("volcano_select_proteins", "Show proteins:",
                                                             choices = c("Significant" = "Significant", "All" = "All", "Protein List" = "Upload"), selected = "Significant"),
                                                conditionalPanel(condition = "input.volcano_select_proteins == 'Upload'",
                                                                 h5("Show proteins that are found in the uploaded list. Upload a csv file with only a list of proteins, no headers."),
                                                                 radioButtons("volcano_file_keytype", "Accession: ",
                                                                              choices = c("Gene Symbol" = "SYMBOL", "UNIPROT" = "UNIPROT", "ENTREZID" = "ENTREZID"), selected = "SYMBOL"),
                                                                 fileInput("volcano_file", "Select List to Upload: ",
                                                                           multiple = TRUE,
                                                                           accept = c("text/csv",
                                                                                      "text/comma-separated-values,text/plain",
                                                                                      ".csv"))),
                                                actionButton("volcano_run", "Run"),
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.volcano_run > 0",
                                                                 h3("Volcano Plot"),
                                                                 plotOutput("volcano"),
                                                                 tags$hr(),
                                                ),
                                       ),
                           ), ## close tabSetPanel
                  ), ## close tabPanel
                  
                  ########### FUNCTIONAL ##########
                  tabPanel("Functional Analysis",
                           tabsetPanel(type="pills",
                                
                                       tabPanel("Gene Ontology",
                                                h5("Over-representation Test: please select a list of proteins."),
                                                radioButtons("go_proteinlist", "Select protein list: ",
                                                             choices = c("Heatmap", "Volcano", "Upload File"),
                                                             selected = "Heatmap"),  
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.go_proteinlist == 'Heatmap'",
                                                                 h5("Use proteins from the heatmap tab that have been organized by cluster."),
                                                                 checkboxGroupInput("go_heatmap_clusters_choices", "Select Clusters:")
                                                ),
                                                conditionalPanel(condition = "input.go_proteinlist == 'Volcano'",
                                                                 radioButtons("go_volcano_choices", "Filter by Fold Change:", choices = c("Positive", "Negative")),
                                                                 uiOutput("go_volcano_groups")
                                                ),
                                                conditionalPanel(condition = "input.go_proteinlist == 'Upload File'",
                                                                 h5("Upload a csv file with only a list of proteins, no headers."),
                                                                 radioButtons("go_file_keytype", "Accession: ",
                                                                              choices = c("Gene Symbol" = "SYMBOL", "UNIPROT" = "UNIPROT", "ENTREZID" = "ENTREZID"), 
                                                                              selected = "SYMBOL"),
                                                                 fileInput("go_file", "Select List to Upload: ",
                                                                           multiple = TRUE,
                                                                           accept = c("text/csv", "text/comma-separated-values,text/plain",".csv"))
                                                ),
                                                tags$hr(),
                                                
                                                radioButtons("GO_ont", "Select Ontology: ", choices = c("Biological Process" = "BP", "Cellular Compartment" = "CC", "Molecular Function" = "MF", "All" = "ALL"),
                                                             selected = "BP"),
                                                radioButtons("go_universe", "Reference List: ", choices = c("All genes" = "All", "This dataset" = "DF"),
                                                             selected = "All"),
                                                actionButton("go_calculate", "Run"),
                                                tags$hr(),
                                                
                                                conditionalPanel(condition = "input.go_calculate > 0",
                                                                 h3("Enrichment Results"),
                                                                 plotOutput("go_dotplot")
                                                                 ), # close heatmap conditional
                                       ), # close GO
                           ), # close Functional Analysis tabsetPanel
                  ), # close tabPanel
    #################
    ) # close mainPanel
  ) # close sidebarLayout
) # close fluidPage
    

################## SERVER ##################
server <- function(input, output, session) {
  while (!is.null(dev.list()))  dev.off()
  ##################### INPUT########################
  
  # Select organism
  orgdb <- reactiveValues(org=org.Hs.eg.db)
  orgdb_string <- reactiveValues(org='org.Hs.eg.db')
  
  # Update checkboxinputs
  observe({
    if(input$organism == "Homo sapiens"){
      orgdb$org <- org.Hs.eg.db
      orgdb_string$org <- 'org.Hs.eg.db'
    }
    else if(input$organism == "Mus musculus"){
      orgdb$org <- org.Mm.eg.db
      orgdb_string$org <- 'org.Mm.eg.db'
    }
    
    choice <- names()
    select1 <- choice[1]
    select2 <- choice[2]
    
    updateCheckboxGroupInput(session = session,
                             inputId = "groups1",
                             choices = choice,
                             selected = select1)
    updateCheckboxGroupInput(session = session,
                             inputId = "groups2",
                             choices = choice,
                             selected = select2)
    updateCheckboxGroupInput(session = session,
                             inputId = "venn_choices",
                             choices = choice,
                             selected = c(select1, select2))
    updateCheckboxGroupInput(session = session,
                             inputId = "groups_heatmap",
                             choices = choice,
                             selected = choice)
    updateCheckboxGroupInput(session = session,
                             inputId = "go_heatmap_groups",
                             choices = choice,
                             selected = choice)
    updateCheckboxGroupInput(session = session,
                             inputId = "gsea_heatmap_groups",
                             choices = choice,
                             selected = choice)
  })
  
  # main data frame, after all data processing
  df <- reactive({
    withProgress(message="Uploaded dataset!", value = 0.3, {
      req(input$fileone)
      df <- read.csv(input$fileone$datapath,
                     header = TRUE,
                     sep = ",",
                     quote = '"',
                     encoding = "UTF-8")
      
      if(is.null(input$fileone) && input$demo > 0) {
        df <- read.csv(demoDataExpression,
                       header = TRUE,
                       sep = ",",
                       quote = '"',
                       encoding = "UTF-8")
      }
      
      incProgress(0.3, "Cleaning data...")
      if (input$na_processing == "Omit") {df <- na.omit(df)}
      else if (input$na_processing == "Replace with 0") {df[is.na(df)] <- 0}
      else if (input$na_processing == "Replace with 0.0001") {df[is.na(df)] <- 0.0001}
      
      # remove a protein if all sample values are ~0
      df <- df[rowSums(df[2:ncol(df)])>0.0001*(ncol(df)-1),]
      
      if(input$zero_processing == "Omit") {df <- df[apply(df!=0, 1, all),]}
      else if(input$na_processing == "Replace with 0.0001") {df[df == 0] <- 0.0001}
      
      if(input$filter_small > 0){
        means <- apply(df[,2:ncol(df)], 1, mean)
        cutoff <- quantile(means, c(input$filter_small/100))
        df2 <- df[,2:ncol(df)]
        df <- df[rowMeans(df2)>cutoff,]
        print(head(df))
      }
      
      incProgress(0.2, "Processing protein names...")
      Protein_names <- df[,1]
      Protein <- c()
      for (i in 1:length(Protein_names)){
        Protein <- append(Protein, strsplit(Protein_names[i], split=input$keytype_split)[[1]][1])
      }
      df <- df[,2:ncol(df)]
      df <- df[,order(metadf()$V2)]
      df <- cbind.data.frame(Protein,df)
      
      if(input$duplicates == "Omit"){
        df <- df %>% distinct(Protein, .keep_all=TRUE)
      }
      colnames(df) <- gsub('\\.','_', colnames(df))
      colnames(df) <- gsub('-','_', colnames(df))
      df[,1] <- gsub(' ','_', df[,1])
      df[,1] <- gsub(',','_', df[,1])
      
      if(input$centermeans == "Quantile"){
        incProgress(0.1, "Normalizing...")
        dfonly <- df[,2:ncol(df)]
        cols <- colnames(dfonly)
        centered_df <- normalize.quantiles(as.matrix(dfonly))
        centered_df <- as.data.frame(centered_df)
        
        df <- cbind.data.frame(df[,1], centered_df)
        colnames(df) <- c("Protein", cols)
      }
    })
    return(df)
  })
  
  # Convert original to a "dictionary" of gene names - by column: UNIPROT, SYMBOL, ENTREZID
  dfprots <- eventReactive(input$submit, {
    withProgress(message="Converting protein IDs...", value=0.1, {
      if(input$keytype == "UNIPROT"){
        current <- df()[,1]
        # convert to other 2 types
        conv <- bitr(current, fromType=input$keytype, toType=c("SYMBOL", "ENTREZID"), OrgDb=orgdb$org, drop=FALSE)
        conv[is.na(conv)] <- '#N/A'
        
        if(input$duplicates == "Omit"){
          symbol <- conv[,c(1,2)] %>% distinct(UNIPROT, .keep_all=TRUE)
          symbol <- symbol$SYMBOL
          entrez <- conv[,c(1,3)] %>% distinct(UNIPROT, .keep_all=TRUE)
          entrez <- entrez$ENTREZ
        }
        else if(input$duplicates == "Keep"){
          symbol <- c()
          entrez <- c()
          
          ### slow version - but avoids cutting unique proteins
          for (i in 1:length(current)){
            incProgress(0.8/length(current), "Checking accessions...")
            match <- conv[current[i] == conv$UNIPROT,]
            if(nrow(match) > 1){
              match <- match[1,]
            }
            symbol <- append(symbol, match[,2])
            entrez <- append(entrez, match[,3])
          }
        }
        conv <- data.frame(UNIPROT = current, SYMBOL = symbol, ENTREZID = entrez)
      }
      else if(input$keytype == "ENTREZID"){
        current <- df()[,1]
        conv <- bitr(current, fromType=input$keytype, toType=c("SYMBOL", "UNIPROT"), OrgDb=orgdb$org, drop=FALSE)
        conv[is.na(conv)] <- '#N/A'
        
        if(input$duplicates == "Omit"){
          symbol <- conv[,c(1,2)] %>% distinct(ENTREZID, .keep_all=TRUE)
          symbol <- symbol$SYMBOL
          uniprot <- conv[,c(1,3)] %>% distinct(ENTREZID, .keep_all=TRUE)
          uniprot <- uniprot$UNIPROT
        }
        else if (input$duplicates == "Keep"){
          symbol <- c()
          uniprot <- c()
          
          for (i in 1:length(current)){
            incProgress(0.8/length(current), "Checking accessions...")
            match <- conv[current[i] == conv$ENTREZID,]
            if(nrow(match) > 1){
              match <- match[1,]
            }
            symbol <- append(symbol, match[,2])
            uniprot <- append(uniprot, match[,3])
          }
        }
        conv <- data.frame(UNIPROT = uniprot, SYMBOL = symbol, ENTREZID = current)
      }
      else if(input$keytype == "SYMBOL"){
        current <- df()[,1]
        conv <- bitr(df()[,1], fromType=input$keytype, toType=c("UNIPROT", "ENTREZID"), OrgDb=orgdb$org, drop=FALSE)
        conv[is.na(conv)] <- '#N/A'
        
        if(input$duplicates == "Omit"){
          uniprot <- conv[,c(1,2)] %>% distinct(SYMBOL, .keep_all=TRUE)
          uniprot <- uniprot$UNIPROT
          entrez <- conv[,c(1,3)] %>% distinct(SYMBOL, .keep_all=TRUE)
          entrez <- entrez$ENTREZ
        }
        else if(input$duplicates == "Keep"){
          uniprot <- c()
          entrez <- c()
          
          for (i in 1:length(current)){
            incProgress(0.8/length(current), "Checking accessions...")
            match <- conv[current[i] == conv$SYMBOL,]
            if(nrow(match) > 1){
              match <- match[1,]
            }
            uniprot <- append(uniprot, match[,2])
            entrez <- append(entrez, match[,3])
          }
        }
        conv <- data.frame(UNIPROT = uniprot, SYMBOL = current, ENTREZID = entrez)
      }
    })
    
    conv[is.na(conv)] <- '#N/A'
    return(conv)
  })
  
  # matrix only - no names of proteins (for quick calculation)
  dfres <- eventReactive(input$submit, {
    obs <- observations()
    dfr <- df()[2:(obs+1)]
    return(dfr)
  })
  
  # transposed dataframe - samples are in rows. df_trans()$group gives the grouping of each sample.
  df_trans <- eventReactive(input$submit, {
    df_t <- base::as.data.frame(base::t(dfres()))
    df_t$group <- rep(names(), times())
    return(df_t)
  })
  
  # pooled standard deviation for each protein
  pooled_sd <- reactive({
    sds <- data.frame(Protein = df()[,1])
    for(i in 1:length(names())){
      mini_df <- df_trans()[df_trans()$group == names()[i], 1:nrow(dfres())]
      sd <- apply(mini_df, 2, sd)
      sd <- sd^2*(times()[i]-1)
      sds <- cbind.data.frame(sds, sd)
    }
    
    weighted <- rowSums(sds[2:ncol(sds)])
    # define pooled var as (n1-1)s1^2 + (n2-1)s2^2 / (n1+n2-2)
    weighted <- weighted/(observations()-length(times()))
    # define pooled sd as sqrt(pooled_var)
    weighted <- sqrt(weighted)
    pooled_sd <- data.frame(Protein=df()[,1], Pooled=weighted)
    
    return(pooled_sd)
  })
  
  ################# DATA AND QUALITY CHECK #################
  # show contents for check
  summary_stats <- reactive({
    
    sw <- function(x){
      if(length(unique(x)) == 1){
        # can't allow all samples to have the same value
        x[1] = x[1] + 0.001
      }
      shapiro <- shapiro.test(x)
      return(shapiro$p.value)
    }
    
    quantiles <- apply(dfres(), 1, quantile)
    
    means <- apply(dfres(), 1, mean)
    stdev <- apply(dfres(), 1, sd)
    pcv <- stdev/means*100
    min <- quantiles[1,]
    q25 <- quantiles[2,]
    median <- quantiles[3,]
    q75 <- quantiles[4,]
    max <- quantiles[5,]
    shapiro <- apply(dfres(), 1, sw)
    normality <- shapiro > 0.05
    
    summary <- cbind.data.frame(df()[,1], means, stdev, pcv, min, q25, median, q75, max, shapiro, normality)
    base::colnames(summary) <- c("Protein", "Average", "SD", "%CV", "Minimum", "Q1", "Median", "Q3", "Maximum", "Shapiro-Wilk P-value", "Normality")
    
    return(summary)
  })
  
  # ratio of proteins that are normally distributed
  ratio_normal <- eventReactive(input$submit, {
    ratio <- 1-(sum(summary_stats()$Normality)/nrow(summary_stats()))
    return(ratio)
  })
  
  ##################### UNIVARIATE ########################
  # density plot of proteins across groups and altogether (merged)
  density_plot <- reactive({
    # samples categorized all together
    prot <- rep("Merge", length(selected_values()))
    combined_df <- cbind.data.frame(Group = prot, Values = t(selected_values()))
    colnames(combined_df) <- c("Group", "Values")
    
    # samples categorized by grouping
    prot2 <- rep(names(), times())
    combined_df2 <- cbind.data.frame(prot2, t(selected_values()))
    colnames(combined_df2) <- colnames(combined_df)
    combined_df <- rbind(combined_df, combined_df2)
    
    title <- paste("Density plot of ", input$protein_viewer_select, sep = "")
    p <- ggplot(combined_df, aes(x=Values, col=Group)) + geom_density() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                                                 axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'white'))
    return(p)
  })
  
  output$density_plot <- renderPlot({
    return(density_plot())
  })
  
  # quantile-quantile plot of proteins across groups and altogether (merged)
  qq_plot <- reactive({
    # samples categorized all together
    prot <- rep("Merge", length(selected_values()))
    combined_df <- cbind.data.frame(Group = prot, Values = t(selected_values()))
    colnames(combined_df) <- c("Group", "Values")
    
    # samples categorized by grouping
    prot2 <- rep(names(), times())
    combined_df2 <- cbind.data.frame(prot2, t(selected_values()))
    colnames(combined_df2) <- colnames(combined_df)
    
    combined_df <- rbind(combined_df, combined_df2)
    
    title <- paste("Quantile-Quantile plot of ", input$protein_viewer_select, sep = "")
    p <- ggplot(combined_df, aes(sample=Values, colour=Group)) + stat_qq() + stat_qq_line() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                                                                     axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'white'))
    return(p)
  })
  
  output$qq_plot <- renderPlot({
    return(qq_plot())
  })
  
  ################# VENN DIAGRAM #########################
  # table output: first column = protein names, rest of the columns = Present(T)/Absent(F) for each group
  venn_df <- reactive({
    req(input$fileone)
    df <- read.csv(input$fileone$datapath,
                   header = TRUE,
                   sep = ",",
                   quote = '"',
                   encoding = "UTF-8")
    
    colnames(df) <- gsub('\\.','_', colnames(df))
    colnames(df) <- gsub('-','_', colnames(df))
    
    # Replace 0s with NA
    df[df <= 0.001] <- NA
    return(df)
  })
  
  # if another file is to be uploaded
  venn_uploadfile <- reactive({
    req(input$venn_uploadfile)
    df <- read.csv(input$venn_uploadfile$datapath,
                   header = FALSE,
                   sep = ",",
                   quote = '"',
                   encoding = "UTF-8")
    return(df[,1])
  })
  
  # table to use for venn diagram
  venn_diagram_table <- reactive({
    withProgress(message="Generating Venn Diagram", value=0.2, {
      df <- venn_df()
      
      # subset dataframe into smaller dataframes
      mini_dfs <- list()
      for (i in 1:length(names())){
        s <- base::subset(metadf(), metadf()$V2 == names()[i])
        mini_dfs <- list.append(mini_dfs, df[,s$V1])
      }
      
      incProgress(0.4, "Examining NAs and 0s...")
      # each smaller dataframe is examined for NAs (T/F) in any of the replicates
      na_s <- list()
      
      # define function where presence must be at least x%
      rowNA <- function(x){
        nas <- !is.na(x)
        ratio <- sum(nas)/length(nas)
        if(ratio < input$filter_presence/100)
          return(TRUE)
        else
          return(FALSE)
      }
      
      for (i in 1:length(mini_dfs)){
        n <- !(apply(mini_dfs[[i]], 1, rowNA))
        na_s <- list.append(na_s, n)
      }
      names(na_s) <- names()
      
      na_s <- as.data.frame(do.call(cbind, na_s))
      d <- cbind.data.frame(df[,1], na_s)
      base::colnames(d) <- c("Protein", names())
      
      # is separate file uploaded?
      if(input$venn_enable_upload == "Disable"){
        File <- rep(FALSE, nrow(d))
      }
      # if yes, call "true" for all proteins in the separate file
      else if(input$venn_enable_upload == "Enable" & length(venn_uploadfile()) > 0){
        if(input$venn_file_keytype != input$keytype){
          uploadfile <- bitr(venn_uploadfile(), fromType=input$venn_file_keytype, toType=input$keytype, OrgDb = orgdb$org)
          uploadfile <- uploadfile[,2]
        }
        else{
          uploadfile <- venn_uploadfile()
        }
        File <- d[,1] %in% uploadfile
      }
      
      d <- cbind.data.frame(d, File)
      return(d)
    })
  })

  venn_diagram <- eventReactive(input$venn_run, {
    if(length(input$venn_choices) + as.numeric(input$venn_enable_upload == "Enable") > 4 | length(input$venn_choices) + as.numeric(input$venn_enable_upload == "Enable") < 2){
      showModal(venn_failedModal())
      return()
    }
    if(length(input$venn_choices) == 2) {
      p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2]), show_percentage=FALSE, stroke_size=0.5, set_name_size=3)
    }
    else if(length(input$venn_choices) == 3) {
      p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2], input$venn_choices[3]), show_percentage=FALSE, stroke_size=0.5, set_name_size=3)
    }
    else if(length(input$venn_choices) == 4) {
      p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2], input$venn_choices[3], input$venn_choices[4]), show_percentage=FALSE, stroke_size=0.5, set_name_size=3)
    }
    return(p)
  })
  
  output$venndiagram <- renderPlot({
    return(venn_diagram())
  })
  
  #################### PRINCIPAL COMPONENTS ANALYSIS ####################
  observeEvent(input$submit, {
    choice <- colnames(df())[2:ncol(df())]
    updateCheckboxGroupInput(session = session,
                             inputId = "pca_sample_choices",
                             choices = choice)
  })
  
  # orders dataframe by coefficient of variation
  cv_sorted <- reactive({
    # subset samples by grouping
    unq <- base::unique(metadf()$V2)
    
    subs <- base::subset(df_trans(), df_trans()$group == unq[1])
    subs <- apply(subs[1:(ncol(subs)-1)], 2, mean)
    
    for (i in 2:length(unq)) {
      s <- base::subset(df_trans(), df_trans()$group == unq[i])
      m <- apply(s[1:(ncol(s)-1)], 2, mean)
      subs <- base::rbind(subs, m)
    }
    
    subs <- t(subs)
    
    # calculate coefficient of variation
    means <- apply(subs, 1, mean)
    stdev <- apply(subs, 1, sd)
    pcv <- stdev/means*100
    
    dfr <- dfres()
    dfr$pcv <- pcv
    
    sorted <- dfr[base::order(-dfr$pcv),]
    return(sorted)
  })
  
  pca <- reactive({
    # only allow between 2 and maximum proteins in PCA
    if(input$pca_topcv > nrow(df()) | input$pca_topcv < 2){
      showModal(pca_failedModal())
      return()
    }
    pca <- prcomp(pca_matrix_t(), center=TRUE, scale.=TRUE)
    
    # takes the rotation value
    df_pca <- base::as.data.frame(pca$x)
    
    # assigns group to each replicate
    if(input$pca_samples == 'Select' & length(input$pca_sample_choices) > 0){
      md <- metadf()[metadf()$V1 %in% input$pca_sample_choices,]
      sum <- data.frame(times=summary(md$V2))
      combined <- data.frame(names=rownames(sum), times=sum$times)
      combined <- combined[combined$times > 0,]
    }
    else if(input$pca_samples == 'All'){
      combined <- data.frame(names=names(), times=times())
    }
    
    df_pca$group <- rep(combined$names, combined$times)
    df_pca$group <- factor(df_pca$group, levels=combined$names)
    
    # assigns a color to each group
    length_times <- length(times())
    color <- c()
    for(i in 1:length(df_pca$group)){
      index <- which(names() == df_pca$group[i])
      color <- append(color, rainbow(length_times)[index])
    }
    df_pca$color <- color
    
    if(input$pc_x == 1 & input$pc_y == 2) {p <- ggplot(df_pca, aes(x=PC1, y=PC2, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times()))) 
    if(input$pc_x == 1 & input$pc_y == 3) {p <- ggplot(df_pca, aes(x=PC1, y=PC3, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 1 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC1, y=PC4, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 1 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC1, y=PC5, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 2 & input$pc_y == 3) {p <- ggplot(df_pca, aes(x=PC2, y=PC3, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 2 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC2, y=PC4, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 2 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC2, y=PC5, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 3 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC3, y=PC4, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 3 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC3, y=PC5, col=group)) + geom_point(shape=19, size=5, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    if(input$pc_x == 4 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC4, y=PC5, col=group)) + geom_point(shape=19, size=53, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
    
    p <- p + theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'white'))
    
    return(p)
  })

  output$pca <- renderPlot({
      return(pca())
  }
  
  ########################### VOLCANO PLOT ############################
  # dataframe of all replicates in group 1
  v1 <- reactive({
    withProgress(message="Subsetting dataframe...", value=0.3, {
      df_1 <- base::subset(df_trans(), df_trans()$group %in% input$groups1)
      v1 <- t(df_1[1:(base::nrow(dfres()))])
    })
    return(v1)
  })
  
  # dataframe of all replicates in group 2
  v2 <- reactive({
    df_2 <- base::subset(df_trans(), df_trans()$group %in% input$groups2)
    v2 <- t(df_2[1:(base::nrow(dfres()))])
    return(v2)
  })
  
  # calculates log2fc between v1 and v2
  Log2FC <- reactive({
    withProgress(message="Calculating log fold changes...", value=0.2, {
      n1 <- apply(v1(), 1, mean)
      n2 <- apply(v2(), 1, mean)
      
      fc <- n2/n1
      Log2FC <- log2(fc)
    })
    return(Log2FC)
  })
  
  # bind pvalues of the difference of means test (parametric/nonparametric) to each protein
  df_pval <- eventReactive(input$volcano_run, {
    l_1 <- ncol(v1())
    l_2 <- ncol(v2())
    
    # sample values and pooled SD for each protein
    vfull <- cbind.data.frame(v1(), v2(), (pooled_sd()$Pooled))
    
    t_testing2 <- function(x, l1, l2) {
      x_1 <- x[1:l1]
      x_2 <- x[(l1+1):(l1+l2)]
      
      # won't allow when all sample values are the same
      if(length(unique(x_1)) == 1 & length(unique(x_2)) == 1){
        return(1)
      }
      if(length(unique(x_1)) == 1 | length(unique(x_2)) == 1){
        x_1[1] <- x_1[1]+0.002
        x_2[1] <- x_2[1]+0.001
      }
      
      # define t-stat as mean difference / pooled SD*(1/l1 + 1/l2)
      pooled <- x[ncol(vfull)]
      meandiff <- abs(mean(x_2) - mean(x_1))
      t_stat <- meandiff/(pooled*((1/l1)+(1/l2)))
      
      # two-sided test, take double of p-value observed
      result <- 2*pt(q=t_stat, df=(observations()-length(times())), lower.tail=FALSE)
      return(result)
    }
    
    wilcoxon_testing <- function(x, l1, l2){
      x_1 = x[1:l1]
      x_2 = x[(l1+1):(l1+l2)]
      
      w <- wilcox.test(x_1, x_2, exact=TRUE)
      return(w$p.value)
    }
    
    withProgress(message="Generating Volcano Plot...", value=0.1, {
      incProgress(0.5, detail = "Calculating p-values...")
      
      if(input$volcano_test_type == "Parametric"){
        Pval <- base::apply(vfull, 1, t_testing2, l1=l_1, l2=l_2)
      }
      else if(input$volcano_test_type != "Parametric"){
        Pval <- base::apply(vfull, 1, wilcoxon_testing, l1=l_1, l2=l_2)
      }
      
      if(input$adjustment != 'none') {
        Pval <- p.adjust(Pval, input$adjustment) 
      }
      
      LogPval <- -1*log10(Pval) 
      
      tt_log <- cbind(Pval, LogPval)
    })
    
    return(tt_log)
  })

  # combines p-value and log2fc into one dataframe
  df_tt <- eventReactive(input$volcano_run, {
    df_tt <- base::cbind(dfres(), Log2FC(), df_pval())
    df_tt_test <- df_tt[order(df_tt$Pval),]
    return(df_tt)
  })

  # combines p-value, log2fc, and significance (i.e., labeled or not?) into one dataframe
  df_sig <- reactive({
    if(input$volcano_select_proteins == "Significant"){
      fc <- df_tt()[,(base::ncol(df_tt())-2)]
      logpval <- df_pval()[,2]
      bind <- cbind(logpval, fc)
      sig <- bind[,1] >= -1*log10(input$sig_cutoff) & abs(bind[,2]) > log2(input$fc_cutoff)
    }
    else if(input$volcano_select_proteins == "All"){
      sig <- rep(TRUE, nrow(df_tt()))
    }
    else if(input$volcano_select_proteins == "Upload"){
      if(input$volcano_file_keytype == "UNIPROT"){
        protein <- dfprots()[,1]
      }
      else if(input$volcano_file_keytype == "SYMBOL"){
        protein <- dfprots()[,2]
      }
      else if(input$volcano_file_keytype == "ENTREZID"){
        protein <- dfprots()[,3]
      }
      sig <- protein %in% volcano_listupload()
    }
    
    df_sig <- cbind(df_tt(), sig)
    
    return(df_sig)
  })

  # plot volcano
  volcano <- eventReactive(input$volcano_run, {
    df_volcano <- df_sig()
    
    g1 <- base::paste(input$groups1, collapse='/')
    g2 <- base::paste(input$groups2, collapse='/')
    title <- base::paste(g1, "vs", g2, sep=' ')
    p <- ggplot(df_volcano) +
      geom_point(mapping=aes(x=Log2FC(),y=LogPval,col=sig), show.legend=FALSE)+labs(x = "Log2(Fold Change)", y = "-Log(P-adj)")+ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'white'))
    return(p)
  })

  output$volcano <- renderPlot({
    return(volcano())
  })
  
  ########################## GENE ONTOLOGY ############################
  # use list = volcano
  go_listvolcano <- reactive({
    # get significant proteins only
    dfs <- cbind.data.frame(dfprots()[,3], df_sig())
    
    ls <- dfs[dfs$Pval <= input$sig_cutoff, ]
    # subset by logfc
    if(input$go_volcano_choices == "Positive"){
      ls <- ls[ls[,(ncol(dfs)-3)] >= log2(input$fc_cutoff),]
    }
    else if(input$go_volcano_choices == "Negative"){
      ls <- ls[ls[,(ncol(dfs)-3)] <= (-1*log2(input$fc_cutoff)),]
    }
    # get entrezid
    ent <- ls[,1]
    return(ent)
  })
  
  # input: list of proteins, output: enrichment object - contains pathway name, pval, gene ratio, and associated genes
  go_enrichment <- eventReactive(input$go_modalrun, {
    enrichment<- function(x, y){
      plot=enrichGO(x, orgdb$org,
                    keyType = "ENTREZID",
                    ont=y,
                    pvalueCutoff = input$sig_cutoff,
                    pAdjustMethod = input$adjustment,
                    qvalueCutoff = 0.05,
                    minGSSize=10,
                    maxGSSize = 2000,
                    readable = FALSE,
                    pool = FALSE)
      p <- plot
      print("enrich done")
      return(p)
    }
    
    enrichment_dataset <- function(x, y){
      plot=enrichGO(x, orgdb$org,
                    keyType = "ENTREZID",
                    ont=y,
                    universe=ds_full,
                    pvalueCutoff = input$sig_cutoff,
                    pAdjustMethod = input$adjustment,
                    qvalueCutoff = 0.05,
                    minGSSize=10,
                    maxGSSize = 2000,
                    readable = FALSE,
                    pool = FALSE)
      p <- plot
      print("enrich with dataset done")
      return(p)
    }
    
    ds_full <- bitr(venn_df()[,1], fromType=input$keytype, toType='ENTREZID', OrgDb=orgdb$org, drop=FALSE)
    ds_full <- ds_full$ENTREZ
    
    if (input$go_proteinlist == "Volcano"){
      print("volcano")
      if(is.null(go_listvolcano())){
        showModal(selection_fail_Modal())
        return()
      }
      list1 <- go_listvolcano()
    }
    
    if(input$go_universe == "All") {
      withProgress(message = 'Calculating enrichment', value = 0.1, {
        enrichment <- enrichment(list1, input$GO_ont)
        incProgress(0.9, detail = "Enrichment done")
      })
    }
    else if(input$go_universe == "DF"){
      withProgress(message = 'Calculating enrichment', value = 0.1, {
        enrichment <- enrichment_dataset(list1, input$GO_ont)
        incProgress(0.9, detail = "Enrichment done")
      })
    }
    return(enrichment)
  })
  
  # shows dotplot: x-axis/size=gene ratio, color=p-value
  go_dotplot <- reactive({
    if(is.null(go_enrichment())) {
      showModal(nosig_Modal())
    }
    if(input$GO_ont == "ALL"){
      p <- enrichplot::dotplot(go_enrichment(), split="ONTOLOGY", showCategory=input$go_number_show) + facet_grid(~ONTOLOGY)
    }
    else{
      descs <- go_similarity()[1:input$go_number_show,2]
      p <- enrichplot::dotplot(go_enrichment(), showCategory=descs)
    }
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title.x=element_blank(), axis.text.y=element_text(size=8)) + 
      scale_y_discrete(label = function(x) stringr::str_trunc(x, 55))
    return(p)
  })
}
