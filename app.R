#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinylive)
library(httpuv)
library(Seurat)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
library(zip)
library(STutility)
library(shiny)
library(bslib)
library(shinyjs)
library(shinybusy)
library(shinyvalidate)
# library(shinyalert)
library(DT)
library(biomaRt)
library("org.Hs.eg.db")
library("org.Ss.eg.db")
library("clusterProfiler")
#library(saveplot)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)
# 
# 
model_choices = c('AMI')
ami_data <- readRDS(file.path('data','AMI','Acute_Clustering.RDS'))
ami_data@meta.data$week <- unlist(lapply(ami_data@meta.data$slide_id, function(x) { unlist(strsplit(x,"_"))[1] }))
ami_data@meta.data$rep <- unlist(lapply(ami_data@meta.data$slide_id, function(x) { unlist(strsplit(x,"_"))[2] }))
ami_data@meta.data$rep[which(ami_data@meta.data$rep=="Rep1")] <- "Transplanted1"
ami_data@meta.data$rep[which(ami_data@meta.data$rep=="Rep2")] <- "Transplanted2"
ami_data@meta.data$sampleid <- 0
ami_data@meta.data$sampleid[which(ami_data@meta.data$rep=="Transplanted1" & ami_data@meta.data$week=="1wk")] <- 1
ami_data@meta.data$sampleid[which(ami_data@meta.data$rep=="Transplanted1" & ami_data@meta.data$week=="2wk")] <- 3
ami_data@meta.data$sampleid[which(ami_data@meta.data$rep=="Transplanted2" & ami_data@meta.data$week=="1wk")] <- 2
ami_data@meta.data$sampleid[which(ami_data@meta.data$rep=="Transplanted2" & ami_data@meta.data$week=="2wk")] <- 4

ami_data_avg <- read.table(file.path("data","AMI","ami_data_avg.csv"), sep="\t")
ami_genes <- sort(unique(unlist(lapply(row.names(ami_data@assays$RNA$data), function(x) { gsub(substr(x, 1,7),"", x) }))))
ami_species <- sort(unique(unlist(lapply(row.names(ami_data@assays$RNA$data), function(x) { substr(x, 1,6) }))))

ami_gene_species <- data.frame( 
                      genes = unlist(lapply(row.names(ami_data@assays$RNA$data), function(x) { gsub(substr(x, 1,7),"", x) })),
                      species = unlist(lapply(row.names(ami_data@assays$RNA$data), function(x) { substr(x, 1,6) }))
                    )
ami_gene_species$species_name <- ami_gene_species$species
ami_gene_species$species_name[which(ami_gene_species$species=='GRCh38')] <- 'GRCh38 (Human)'
ami_gene_species$species_name[which(ami_gene_species$species=='Ssus11')] <- 'Ssus11 (Pig)'


ami_clusters <- read.table(file.path("data","AMI","Acute_deMarkers_ALL_Clusters.txt"), sep="\t",row.names = 1, header = T)
#row.names(ami_clusters) <- NULL

ami_human_clusterids <- c(9,10)




gene_categories <- read.csv(file.path("data","categories_genes.csv"), header=T)

cell_timepoints <- paste(as.character(c(1,4,12)),'wk', sep = "")

human_gene2ensembl <- read.table(file.path("data","Homo_sapiens.GRCh38.102_gene_annotation_table.txt"), header=T)


pig_gene2ensembl <- read.table(file.path("data","Pig_EnsemblID.txt"), header=T)
#gene_id GeneSymbol are the column names 


# Define UI for application that draws a histogram
ui <- 
  tagList(
    navbarPage(h3(strong("MI-HEART-ST")), 
      
    tabPanel(title=h3(strong("Introduction", style="color:RoyalBlue")), 
             value="intro",
             fluidPage(
               h3(strong("BACKGROUND", style="color:RoyalBlue")),
               hr(),
                     p(style="text-align: justify; font-size:13pt; font-family: arial",
                     "The human pluripotent stem cells derived cardiovascular progenitors (CVPs) differentiated 
                     on the laminin 521+221 matrix were transplanted them into acute and chronic Myocardial Infarcted (MI) 
                     pig hearts (AMI and CMI).  The time-series spatial transcriptomes were generated to characterize these 
                     human CVP cells at AMI 1- and 2- and at CMI 1-, 4- and 12 weeks post-transplantation. 
                     The AMI and CMI transcriptome atlas is provided in this shiny app. "),
                     br(),
                     tags$div(img(src = "models-1.png", alt="AMI model", height = 300, width = 700,),img(src = "models-2.png", alt="CMI model", height = 300, width = 700,)),
                     br(), br(),

                h3(strong("H&E STAINING OF AMI AND CMI CVPs TRANSPLANTED HEART TISSUES",  style="color:RoyalBlue")),
                hr(), 
                      p(style="text-align: justify; font-size:13pt; font-family: arial", strong("I"),"--> Infarcted"),
                      p(style="text-align: justify; font-size:13pt; font-family: arial", strong("NI"),"--> Non-Infarcted"),
                      p(style="text-align: justify; font-size:13pt; font-family: arial", strong("Red dotted box"),"--> CVP Transplanted"),
                      br(), 
                      # tags$img(height = 500, width = 300,src = "acute-1.png"),
                      tags$div(img(src = "acute-1.png",height = 600, width = 550,), img(src = "chronic-1.png", height = 600, width = 550,)),
                br(), br(),
                   
                h3(strong("FUNCTIONALITIES OF EACH TAB", style="color:RoyalBlue")),
                hr(),
                   tags$ul(
                     tags$ul(style="text-align: left; font-size:13pt; font-family: arial",strong("1. Myocardial Infarction (MI) Models")), 
                     p(style="font-size:13pt; font-family: arial",
                     "This page enables users to customize their selection of MI (Myocardial Infarction) models, species (GRCh38-human, Ssus11-pig), 
                     tissue sections at various timepoints, and up to six genes of interest. A link is provided for users to view the functional spot 
                     clusters corresponding to the selected model and timepoint. Upon gene selection, the page generates spatial spot plots of the chosen 
                     genes alongside an H&E staining image, emphasizing regions of infarction (I), non-infarction (NI), and areas with CVP (Cardiovascular 
                     Progenitor Cells) transplantation, marked by a red-colored box. Furthermore, if the selected gene is involved in any signaling pathways, 
                     an extra link is available to navigate to a signaling page displaying associated ligand and receptor pairs. Lastly, a download link is 
                     provided for users to obtain the spot plot for the selected genes in PNG format."),
                     br(),
                     tags$ul(style="text-align: justify; font-size:13pt; font-family: arial",strong("2. Functional Spot Clusters")),
                     p(style="font-size:13pt; font-family: arial", "This page facilitates both functional spot cluster selection and gene selection. 
                       With cluster selection, users can choose the model and spot clusters of interest by examining the spot clustered overlayed H&E 
                       plots provided on the left panel. Upon selecting a model and cluster, functional analysis dot plots are displayed in the right panel. 
                       In the gene selection tab, users can explore differentially expressed genes across spot clusters. This tab also enables users to 
                       search for genes of interest and view their statistical significance and expression across different spot clusters."),
                     br(),
                     tags$ul(style="text-align: justify; font-size:13pt; font-family: arial",strong("3. Acute vs Chronic MI Model Change")), 
                     p(style="font-size:13pt; font-family: arial", "AMI and CMI share only one common timepoint, which is 1-week post CVP transplantation. 
                     This page allows users to investigate marker expression and differentially expressed genes between AMI and CMI at this specific timepoint. 
                     In the marker(s) expression tab, users can select the species and gene of interest to generate spatial spot plots of the chosen genes alongside 
                     H&E staining image highlighting regions of infarction (I), non-infarction (NI), and areas with CVP (Cardiovascular Progenitor Cells) transplantation, 
                     indicated by a red-colored box. In the differential gene expression tab, users can explore genes that are differentially expressed across spot clusters, 
                       as depicted in the spot clusters overlayed H&E images in the left panel. This tab also facilitates searching for genes of interest and provides information 
                       on their statistical significance and expression across different spot clusters."),
                     br(),
                     tags$ul(style="text-align: justify; font-size:13pt; font-family: arial",strong("4. Signaling")), 
                     p(style="font-size:13pt; font-family: arial", "This page allows users to choose specific timepoint from CMI revealing the enriched ligands 
                       and receptors identified using the cell-cell communication analysis. Users can then select a particular ligand receptor pair from the 
                       dropdown to view their saptial localizations."),
                     
                    ),
                br(), br(),
                
              h3(strong("CITATION",  style="color:RoyalBlue")),
              hr(),
                   p(style="text-align: justify; font-size:13pt; font-family: arial",
                     "Swarnaseetha Adusumalli, Samantha Lim, Vincent Ren, LiYen Chong, Roy Tham, YeLei, Yibin Wang, Enrico Petretto, Karl Tryggvason, Lynn Yap"),
                   p(style="text-align: justify; font-size:13pt; font-family: arial",
                     
                     "bioRxiv 2023.06.10.544480; doi:", a(href="https://doi.org/10.1101/2023.06.10.544480","https://doi.org/10.1101/2023.06.10.544480" )),
              br(), br(),
                   
              
              h3(strong("PREVIOUS STUDIES","(", em("bulk and single-cell RNA sequencing - Mouse and Pig models"),")", style="color:RoyalBlue")),
              hr(),
               p(style="text-align: justify; font-size:13pt; font-family: arial",
                "Yap,  L. et al. Pluripotent stem cell-derived committed cardiac progenitors remuscularize damaged ischemic hearts and improve their function in pigs. NPJ Regen Med 8, 26 (2023)"),
              p(style="text-align: justify; font-size:13pt; font-family: arial",
               "doi:", a(href="https://doi.org/10.1038/s41536-023-00302-6","https://doi.org/10.1038/s41536-023-00302-6" )),
              br(),
              p(style="text-align: justify; font-size:13pt; font-family: arial",
                "Yap,  L. et al. In Vivo Generation of Post-infarct Human Cardiac Muscle by Laminin- Promoted Cardiovascular Progenitors. Cell Rep 26, 3231–3245 e3239 (2019)"),
              p(style="text-align: justify; font-size:13pt; font-family: arial",
                "doi:", a(href="https://doi.org/10.1016/j.celrep.2019.02.083","https://doi.org/10.1016/j.celrep.2019.02.083" )),
                        br(), br(), br()
              
                   
              # tags$blockquote("Shiny-Box is still under continuous development. 
              #                   Please look forward to future updates!")
             )
        ),
      #First Tab
    tabPanel(title=h3(strong("Myocardial Infarction (MI) Models",  style="color:RoyalBlue")),
             theme = shinythemes::shinytheme("cerulean"),
             value="mimodels",
             useShinyjs(),
             add_busy_bar(color = "red", height = "8px"),
             column(4, 
                    column(6,
                           selectInput("model_choice",h4(strong("Select the MI Model:")), model_choices, selected=model_choices[1]),
                           selectInput("species_choice", h4(strong("Select a Species:")), choices=NULL),
                           selectInput("timepoint_choice", h4(strong("Select a Timepoint:")), choices=NULL),
                           br(),
                           #actionLink("spot-cluster-pagelink", "Functional Spot Clusters"),
                           tags$a(
                             id="spot-clusters-pagelink",
                             #shinyLink(to = "spotclusters", label = "Functional Spot Clusters"), "."
                           ),br(),
                           downloadButton("downloadGenePlots",h4(strong("Click here to download plots")))
                    ),
                    column(6,
                           selectizeInput("gene_category", h4(strong("Select a Gene Category:")), 
                                          choices = c("ALL", sort(unique(gene_categories$group )))),
                           selectizeInput("gene_choice", h4(strong("Select Gene(s) of Interest:")), choices = NULL, multiple=TRUE,
                                          options = list(maxItems = 6)),
                           p(id="gene_reference_links")
                    ),
                    fluidRow(br(),br(),br()),
                    h3(strong("H & E STAINING"), style="text-align:center;font-family: arial; display:block; font-size:20pt"),
                    br(), br(),
                    imageOutput("heplot_rep1", height="auto"),br(),
                    fluidRow(
                      p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("I"),"--> Infarcted")),
                      p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("NI"),"--> Non-Infarcted")),
                      p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("Red dotted box"),"--> CVP Transplanted")),
                    )
             ),
             column(8,
                    fluidRow(uiOutput("featureoverlayPlot")),
                    style = 'border-left: 4px solid'
                    
             )
    ),
    tabPanel(
      title = h3(strong("Functional Spot Clusters", style="color:RoyalBlue")),
      theme = shinythemes::shinytheme("cerulean"),
      value="spotclusters",
      add_busy_bar(color = "red", height = "8px"),
      column(
        4,
        selectInput("deg_model", h4(strong("Select the MI Model:")),
                    model_choices, selected = NULL),
        br(),
        selectInput("deg_cluster", HTML("<h4> <strong> Select a Spot Cluster: <br/> (as shown in H&E below)</strong> </h4>"), choices = NULL),
        br(),br(),
        #fluidRow(
          imageOutput("deg_clusterplot",width = "auto"), 
          br(),br(),br(),br(),br(),
          br(),br(),br(),
          #height = "400",
          #width = "auto"
        #),
        #fluidRow(
          h4(strong(textOutput("deg_clusterplot_text"))), br(),br()
        #)
        ,
        #fluidRow(
          downloadButton("downloadDegPlots", h4(strong("Click here to download plots")))
        #)
      ),
      
      column(
        8,
        tabsetPanel(
          tabPanel(h3(strong("Cluster Selection")),
                   fluidRow(
                     column(6,
                            h4(strong("Gene Ontology (Biological processes)"), style= "height: 50px; width: 50%; text-align:center; font-size: 20px; display: block;}")
                     ),
                     column(6,
                            fluidRow(h4(strong("KEGG Pathways")), style= "height: 50px; width: 50%; text-align:center; font-size: 20px; display: block;}")
                     )
                   ),
                   fluidRow(
                    tags$div(
                      column(6,
                            #h4(strong("Gene Ontology (Biological processes)"), style= "height: 50px; width: 50%; text-align:center; font-size: 20px; display: block;}"),
                            imageOutput("deg_goplot",inline=TRUE, width="100px",height="100px")
                      ),
                      column(6,
                            #fluidRow(h4(strong("KEGG Pathways")), style= "height: 50px; width: 50%; text-align:center; font-size: 20px; display: block;}"),
                            imageOutput("deg_pathwayplot", inline=TRUE, width="100px",height="100px"),
                      )
                    )
                  )
                #),
          ),
          tabPanel(h3(strong("Gene Selection")),
                   fluidRow(DTOutput("deg_table")),
                   fluidRow(uiOutput("deg_table_featureplot"), height = "auto")
          )
        )
      )
    ),
    tabPanel(
      title= h3(strong("Signaling", style="color:RoyalBlue")), 
      value = "signalling",
      add_busy_bar(color = "red", height = "8px"),
      fluidRow(
        column(3,
               selectInput("cell_week", "Select a CMI timepoint of interest",choices = cell_timepoints ,selected = NULL),
               selectizeInput("cell_interactions", "Select a ligand-receptor pair from this dropdown", 
                              choices = NULL, multiple=TRUE,
                              options = list(maxItems = 6)),
               br(),br(),
               downloadButton("cell_downloadBtn","Download Plots"),
               fluidRow(br(),br(),br()),
               h3(strong("H & E STAINING"), style="text-align:center;font-family: arial; display:block; font-size:20pt"),
               br(), br(),
               imageOutput("heplot_signalling", inline=TRUE, width="60px",height="60px"),br(),
               fluidRow(
                 p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("I"),"--> Infarcted")),
                 p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("NI"),"--> Non-Infarcted")),
                 p(style="text-align: justify; font-size:13pt; font-family: arial", h4(strong("Red dotted box"),"--> CVP Transplanted")),
               )
        ),
        column(9,
               fluidRow(
                 column(8,
                        imageOutput("cell_clusterPlot", inline=TRUE, width="100px",height="100px"),
                        imageOutput("cell_interactionPlot", inline=TRUE, width="100px",height="100px")
                 ),
                 column(4,
                        imageOutput("cell_heatmap", inline=TRUE, width="100px",height="100px")
                 )
               )
        )
      ),
      fluidRow(
        column(3, br()),
        column(4, 
               wellPanel(fluidRow(uiOutput("cell_featureplot"), height = "auto"))
        )
      )
    ),
    position ="static-top", # = bslib::bs_theme(bootswatch = "minty")
  ),
  tags$script(src = "customHref.js")
  )

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
    observeEvent(input$model_choice, {
      selmodel <- input$model_choice
      if(selmodel == "AMI"){
        timechoices <- unique(ami_data@meta.data$week)
        #specieschoices <- sort(unique(ami_gene_species$species))
        specieschoices <- as.list(split(unique(ami_gene_species$species), unique(ami_gene_species$species_name)))
      }
      else{
        timechoices <- character(0)
        specieschoices <- character(0)
      }
      updateSelectInput(session, "timepoint_choice", choices = timechoices)
      updateSelectInput(session, "species_choice", choices = specieschoices)
      updateSelectInput(session, "deg_model", selected = input$model_choice)
      #html("spot-clusters-pagelink", HTML("<a onclick=","customHref('spotclusters')" ,">", "SpotClusters","</a>"))
      html("spot-clusters-pagelink", "<a onclick=spotHref('spotclusters')><h4>Click here to view spot clusters</h4></a>")
      
    })
    
    observe({
      selmodel <- input$model_choice
      selspecies <- input$species_choice
      selcategory <- input$gene_category
      seltime <- input$timepoint_choice
      if(selcategory=="ALL"){
        if(selmodel=="AMI"){
          genechoices <- sort(unique(ami_gene_species$genes[ami_gene_species$species==selspecies]))
        }
        else{
          genechoices <- character(0)
        }
      }else{
        genechoices <- sort(unique(gene_categories$gene[gene_categories$group==selcategory]))
      }
      updateSelectizeInput(session, "gene_choice", choices = genechoices)
    })
    

    output$heplot_rep1 <- renderImage({
      req(input$model_choice)
      req(input$timepoint_choice)
      #tissue <- input$tissue_choice
      tissue <- "Transplanted1"
      filename <- normalizePath(file.path('data',"HandE_images",
                                          paste(input$model_choice, '-', input$timepoint_choice, '-', tissue, '.jpeg', sep='')))
      #if(file.exists(filename)){
        #plot_file_path = file.path(getwd(),"data","HandE_images", filename)
        # Return a list containing the filename and alt text
        list(src = filename, height = 280)
      #}
    },deleteFile = FALSE)
    
    
    plot_featureoverlay <- function(model_name, species_name, timepoint, genelist){

      if(model_name=="AMI"){
        sel_sids <- unique(ami_data@meta.data$sampleid[which(ami_data@meta.data$week==timepoint)])
        genes <- genelist[genelist%in%ami_gene_species$genes[ami_gene_species$species==species_name]]
        validate(need(length(genes)>0, paste0("Selected genes ",genes," are not expressed")))
        sel_genes <- paste0(species_name,"-",genes)
        overlayplot <- FeatureOverlay(ami_data,
                                      features = sel_genes,   #Pgam2, Cox6a2, and Fabp3 ACADVL ESRRA CD36,HADHA
                                      pt.size = 5,
                                      cols = c("dark blue", "cyan", "yellow", "red", "dark red"),
                                      dark.theme = F,
                                      type = "raw",
                                      sampleids = sel_sids,
                                      ncols = 2,
                                      #layout.by.feature =T,
                                      #min.cutoff = 'q1',
                                      show.sb=F
        )
      }
      else{
        sel_sids <- unique(cmi_data@meta.data$sampleid[which(cmi_data@meta.data$week==timepoint)])
        genes <- genelist[genelist%in%cmi_gene_species$genes[cmi_gene_species$species==species_name]]
        validate(need(length(genes)>0, paste0("Selected genes ",genes," are not expressed")))
        sel_genes <- paste0(species_name,"-",genes)
        # cat("CMI selected genes: ",sel_genes,"\n")
        overlayplot <- FeatureOverlay(
          cmi_data,
          features = sel_genes,
          pt.size = 5,
          cols = c("dark blue", "cyan", "yellow", "red", "dark red"),
          dark.theme = F,
          type = "raw",
          sampleids = sel_sids,
          ncols = 2,
          #layout.by.feature = T,
          #min.cutoff = 'q1',
          show.sb = F
        )
      }

      return(overlayplot)

    }
    
    # gene_choice_isolated <- isolate(input$gene_choice)
    
    ## update the gene reference links according to the selected genes and species 
    
    generate_ref_links <- function(gene_symbol, genome){
      if(genome=="GRCh38"){
       ensembl_id <- as.character(human_gene2ensembl$gene_id[human_gene2ensembl$GeneSymbol==gene_symbol])
       ensembl_link = paste0("https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",ensembl_id)
       proteinAtlas_link = paste0("https://www.proteinatlas.org/",ensembl_id)
       geneCards_link = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene_symbol)
       gene_ref_link = paste0("<font size='3'>", gene_symbol,"</font> | ","<a href=\'",ensembl_link, "\', target='_blank', style='font-size:15px;'><strong>ENSEMBL</strong></a> | ",
                              "<a href=\'",proteinAtlas_link, "\', target='_blank', style='font-size:15px;'><strong>ProteinAtlas</strong></a> | ",
                              "<a href=\'", geneCards_link, "\', target='_blank', style='font-size:15px;'><strong>GeneCards</strong></a>")
      }
      else if(genome=='Ssus11'){
        ensembl_id <- as.character(pig_gene2ensembl$gene_id[pig_gene2ensembl$GeneSymbol==gene_symbol])
        ensembl_link = paste0("https://www.ensembl.org/Sus_scrofa/Gene/Summary?db=core;g=",ensembl_id)
        gene_ref_link = paste0("<font size='3'>", gene_symbol,"</font> | ","<a href=\'",ensembl_link, "\', target='_blank', style='font-size:15px;'><strong>ENSEMBL</strong></a> ")
      }
      else{
        gene_ref_link = ''
      }
      return(gene_ref_link)
    }
    
    check_signal_link <- function(geneSymbol, modelName, timePoint){
      signal_ref_link = ''
      if(modelName=='CMI'){
        cc_file = normalizePath(
          file.path('data',"cell-cell_communications",
                    paste(timePoint, '_', 'LR', '_', 'cellchat', '.txt', sep=''))
        )
        data = read.table(cc_file, sep="\t", header=T)
        genelist <- unique(unlist(strsplit(paste0(unique(data$interaction_name), collapse="_"),"_")))
        if(any(geneSymbol%in%genelist)){
          signal_ref_link = paste0(" | <a onclick=signalHref('signalling','", geneSymbol, "')><strong>Signalling</strong></a>")
        }
      }
      return(signal_ref_link)
    }
    
    observeEvent(input$gene_choice, {
      req(input$model_choice)
      req(input$timepoint_choice)
      req(input$gene_choice)
      req(input$species_choice)
      final_html_string <- ''
      for(gene in input$gene_choice){
        gene_string = generate_ref_links(gene, input$species_choice)
        signal_string = check_signal_link (gene, input$model_choice, input$timepoint_choice)
        final_html_string = paste0(final_html_string, gene_string, signal_string, '<br><br>')
      }
      html("gene_reference_links", final_html_string)
      show("gene_reference_links")
    })
    
    observe({
      toggle(id="gene_reference_links", condition = shinyvalidate::input_provided(input$gene_choice))
    })
    
    
    observeEvent(input$timepoint_choice,{
      req(input$model_choice)
      if(input$model_choice=="CMI"){
        updateSelectInput(session, "cell_week", selected = input$timepoint_choice )
      }
    })
    
    observeEvent(input$clicked_spot_gene,{
      if(any(grepl(input$clicked_spot_gene, unique(cellchat_data()$interaction_name)))){
        selitem <- unique(cellchat_data()$interaction_name)[grepl(input$clicked_spot_gene, unique(cellchat_data()$interaction_name))][1]
        updateSelectInput(session, "cell_interactions", selected = selitem)
      }
    })
    
    output$featureoverlayPlot <- renderUI({
      req(input$model_choice)
      req(input$timepoint_choice)
      req(input$species_choice)
      req(input$gene_choice)
      #num_plots <- seq(1:length(input$gene_choice))
      plot_output_list <- lapply(input$gene_choice,
                                 function(g){
                                   plotname <- paste("overlayplot", g,sep="_")
                                   plotOutput(plotname)
                                 })
      do.call(tagList, plot_output_list)
    })
    
    plot_gene_names <- reactiveValues(genes=c())
    observe({
      req(input$gene_choice)
      for(gene_ix in 1:length(input$gene_choice)){
        local({
          icurrent <- input$gene_choice[gene_ix]
          plotname <- paste("overlayplot", icurrent, sep="_")
          output[[plotname]] <- renderPlot({
            plot_featureoverlay(input$model_choice, input$species_choice, input$timepoint_choice, icurrent)
            # plt_overlay_reactive()
          })
        })
      }
      
    })
    # 
    observe({
      toggle(id="downloadGenePlots", condition = !is.null(input$gene_choice))
    })
    
    
    plot_forDownloadGenePlots <- function(model, species, timepoint, geneid){
      tryCatch(
        {
          plt <- plot_featureoverlay(model, species, timepoint, geneid)
          return(plt)
        },
        error = function(e) {
          message(paste0("No Expression found for ",geneid))
          # print(e)
          return(NA)
        },
        warning = function(e){
          message(paste0("No Expression found for ",geneid))
          # print(e)
          return(NA)
        }
      )
    }
    
    output$downloadGenePlots <- downloadHandler(
      filename = function() { paste("gene_featureplots_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".zip", sep = "")},
      content = function(file){
        newTmpDir <- tempfile()
        if(dir.create(newTmpDir)){
          for(gene in input$gene_choice){
            # png(fileName)
            feature_geneplot <- plot_forDownloadGenePlots(input$model_choice, input$species_choice, input$timepoint_choice, gene)
            # dev.off()
            if(any(!is.na(feature_geneplot))){
              fileName <- file.path(newTmpDir, paste0("FeatureOverlayplot_", gene, ".png"))
              ggsave(filename = fileName, 
                     plot = feature_geneplot, device = "png",
                     width = 60, height = 30, units = "cm")
            }
          }
          zipr(file, files=list.files(newTmpDir, full.names = TRUE))
          # unlink(newTmpDir, recursive = TRUE)
        }
        unlink(newTmpDir, recursive = TRUE)
      },
      contentType="application/zip"
    )
    
    ######### DEG page 
    
    observe({
      model <- input$deg_model
      if(model=="AMI"){
        #cluster_options <- c("ALL", as.character(unique(ami_clusters$cluster)))
        cluster_options <- as.character(unique(ami_clusters$cluster))
        cluster_selected <- 9
        text_out <- "CVPs transplanted clusters: 9,10"
      }else{
        #cluster_options <- c("ALL", as.character(unique(cmi_clusters$cluster)))
        cluster_options <- as.character(unique(cmi_clusters$cluster))
        cluster_selected <- 7
        text_out <- 'CVPs transplanted clusters: 3,7,8,12'
      }
      updateSelectInput(session, "deg_cluster", choices = cluster_options, selected = cluster_selected)
      output$deg_clusterplot_text <- renderText({text_out})
    })
    
    output$deg_clusterplot <- renderImage({
      selmodel <-input$deg_model
      
      filename <- normalizePath(file.path('data',selmodel,
                                          paste(selmodel,'_clustering', '.png', sep='')))
      
      list(src = filename, height = 550)
      
    },deleteFile = FALSE)
    
    
    filter_table <- function(model_name, clusterid){
      if(model_name=='AMI'){
        data <- ami_clusters
      }else{
        data <- cmi_clusters
      }
      if(clusterid == 'ALL'){
        seldata <- data
      }else{
        clusterid <- as.numeric(clusterid)
        seldata <- data[data$cluster%in%clusterid ,]
      }
      resdata <- seldata %>% dplyr::select(-c('p_val','pct.1','pct.2','gene'))
      return(resdata)
    }
    
    plot_featureoverlay_deg <- function(model_name, species_name, genelist){
      
      if(model_name=="AMI"){
        # sel_sids <- unique(ami_data@meta.data$sampleid[which(ami_data@meta.data$week==timepoint)])
        sel_sids <- unique(ami_data@meta.data$sampleid)
        genes <- genelist[genelist%in%ami_gene_species$genes[ami_gene_species$species==species_name]]
        validate(need(length(genes)>0, "Selected genes are not expressed"))
        sel_genes <- paste0(species_name,"-",genes)
        #sel_genes_split <- split(sel_genes, ceiling(seq_along(sel_genes)/2))
        # overlayplot <- lapply(sel_genes, function(x){
        overlayplot <- FeatureOverlay(ami_data,
                         features = sel_genes,   #Pgam2, Cox6a2, and Fabp3 ACADVL ESRRA CD36,HADHA
                         pt.size = 5,
                         cols = c("dark blue", "cyan", "yellow", "red", "dark red"),
                         dark.theme = F,
                         type = "raw",
                         sampleids = sel_sids,
                         ncols = 4,
                         #layout.by.feature =T,
                         #min.cutoff = 'q1',
                         show.sb=F
          )
        # })
      }
      else{
        # sel_sids <- unique(cmi_data@meta.data$sampleid[which(cmi_data@meta.data$week==timepoint)])
        sel_sids <- unique(cmi_data@meta.data$sampleid)
        genes <- genelist[genelist%in%cmi_gene_species$genes[cmi_gene_species$species==species_name]]
        validate(need(length(genes)>0, "Selected genes are not expressed"))
        sel_genes <- paste0(species_name,"-",genes)
        #sel_genes_split <- split(sel_genes, ceiling(seq_along(sel_genes)/2))
        #lapply(sel_genes_split, function(x){
        # overlayplot <- lapply(sel_genes, function(x){
        overlayplot <-  FeatureOverlay(cmi_data,
                         features = sel_genes,   #Pgam2, Cox6a2, and Fabp3 ACADVL ESRRA CD36,HADHA
                         pt.size = 5,
                         cols = c("dark blue", "cyan", "yellow", "red", "dark red"),
                         dark.theme = F,
                         type = "raw",
                         sampleids = sel_sids,
                         ncols = 6,
                         #layout.by.feature = T,
                         #min.cutoff = 'q1',
                         show.sb=F
          )
        # })
      }
      
      return(overlayplot)
      
    }
    
    
    deg_tbl <- reactive({
      req(input$deg_cluster)
      req(input$deg_model)
      filter_table(input$deg_model, input$deg_cluster)
    })
    
    observe({
      req(input$deg_cluster)
      req(input$deg_model)
      #deg_tbl <- filter_table(input$deg_model, input$deg_cluster)
      output$deg_table <- renderDT(server=FALSE,{
        datatable(deg_tbl(),extensions=c('Select','Buttons'),
                    options = list(scrollX=TRUE,
                                   scrollY=TRUE,
                                   lengthMenu = c(10),
                                   pageLength = 15,
                                   paging = TRUE, 
                                   select=list(style = 'multi+shift', items = 'row'),
                                   searching = TRUE,
                                   fixedColumns = TRUE,
                                   autoWidth=TRUE, 
                                   dom = 'Blfrtip', 
                                   escape = c('gene'),#'p_val','pct.1','pct.2',
                                   buttons = c('csv','excel','selectNone')),
                  # selection = list(mode = 'multiple', target='row', selected = deg_table_reset$sel),
                  selection = 'none',
                  callback = JS("table.on( 'select', function ( e, dt, type, ix ) {
                                   var selected = dt.rows({selected: true});
                                   if ( selected.count() > 6 ) {
                                      dt.rows({selected: true}).deselect();
                                   }
                                } );")
        )
      }) 
    })
    
    #-------------------
      output$deg_table_featureplot <- renderUI({
        req(input$deg_table_rows_selected)
        num_plots <- seq(1:length(input$deg_table_rows_selected))
        plot_output_list <- lapply(num_plots,
                                   function(g){
                                     plotname <- paste("overlayplot", g,sep="_")
                                     plotOutput(plotname)
                                   })
        do.call(tagList, plot_output_list)
      }) 
  
      observe({
        req(input$deg_model)
        req(input$deg_cluster)
        # deg_table_reset$sel <- input$deg_table_rows_selected
        for(gene_ix in seq(1:length(input$deg_table_rows_selected))){
          local({
            icurrent <- input$deg_table_rows_selected[gene_ix]
            icurrent_species <- deg_tbl()[icurrent,"species"]
            icurrent_gene <- deg_tbl()[icurrent, "geneid"]
            plotname <- paste("overlayplot", gene_ix, sep="_")
            output[[plotname]] <- renderPlot({
              plot_featureoverlay_deg(input$deg_model, icurrent_species, icurrent_gene)
            })
          })
        }
      })
    
    observe({
      toggle(id="downloadDegPlots", condition = !is.null(input$deg_table_rows_selected))
    })
    
    
    plot_forDownloadDegPlots <- function(model, species, geneid){
      tryCatch(
        {
          plt <- plot_featureoverlay_deg(model, species, geneid)
          return(plt)
        },
        error = function(e) {
          message(paste0("Cannot plot for ",geneid))
          # print(e)
          return(NA)
        },
        warning = function(e){
          message(paste0("warning while plotting for ",geneid))
          # print(e)
          return(NA)
        }
      )
    }
    
    
    output$downloadDegPlots <- downloadHandler(
      filename = function() { paste("DEG_featureplots_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".zip", sep = "")},
      content = function(file){
        newTmpDir <- tempfile()
        if(dir.create(newTmpDir)){
          for(gene_ix in seq(1:length(input$deg_table_rows_selected))){
            local({
              icurrent <- input$deg_table_rows_selected[gene_ix]
              icurrent_species <- deg_tbl()[icurrent,"species"]
              icurrent_gene <- deg_tbl()[icurrent, "geneid"]
              feature_geneplot <- plot_forDownloadDegPlots(input$deg_model, icurrent_species, icurrent_gene)
              if(any(!is.na(feature_geneplot))){
                fileName <- file.path(newTmpDir, paste0("FeatureOverlayplot_", icurrent_gene, ".png"))
                ggsave(filename = fileName, 
                       plot = feature_geneplot, device = "png",
                       width = 60, height = 20, units = "cm")
              }
            })
          }
          zipr(file, files=list.files(newTmpDir, full.names = TRUE))
        }
        unlink(newTmpDir, recursive = TRUE)
      },
      contentType="application/zip"
    )
    #----------------
    # output$deg_table_featureplot <- renderText({row.names(deg_tbl())[input$deg_table_rows_selected]})
    
    go_enrich <- function(gene_ids, types, org_db){
      eid <- bitr(gene_ids, fromType = "SYMBOL", toType = types, OrgDb = org_db)
      ego_BP <- enrichGO(gene = eid$ENTREZID, #as.character(unique(inp_sub$EnsemblID)),
                         OrgDb         = org_db,
                         #keytype       = 'ENTREZID',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         #pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
      return(ego_BP)
    }
    
    kegg_enrich <- function(gene_ids, types, org_db, orgname){
      eid <- bitr(gene_ids, fromType = "SYMBOL", toType = types, OrgDb = org_db)
      ego_kegg <- enrichKEGG(gene = eid$ENTREZID, #as.character(unique(inp_sub$EnsemblID)),
                             org=orgname,
                             #keytype       = 'ENTREZID',
                             #ont           = "BP",
                             pAdjustMethod = "BH",
                             #pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05)
      return(ego_kegg)
    }
    
    
    output$deg_goplot <- renderImage({
      req(input$deg_cluster)
      req(input$deg_model)
      validate(
        need(input$deg_cluster != "ALL", "select a cluster for GO plot")
      )
      
      path_name = paste0(input$deg_model,"_GOplot_clus_", as.character(input$deg_cluster), ".png")
      
      filename <- normalizePath(file.path('data','GOplots', path_name))
      
      list(src = filename,  height = 600)
      
      #list(src = img, width = 100, height = 100)
      
    },deleteFile = FALSE)
    
    output$deg_pathwayplot <- renderImage({
      req(input$deg_cluster)
      req(input$deg_model)
      validate(
        need(input$deg_cluster != "ALL", "select a cluster for KEGG plot")
      )
      
      path_name = paste0(input$deg_model,"_KEGGplot_clus_", as.character(input$deg_cluster), ".png")
      
      filename <- normalizePath(file.path('data','KEGGplots', path_name))
      
      list(src = filename,  height = 600)
      
      #list(src = img, width = 100, height = 100)
      
    },deleteFile = FALSE)
    
    ########## server components for AMIvsCMI page
    
    
    observe({
      req(input$cmp_species)
      species_genes <- unique(ami_vs_cmi_genetable$genes[ami_vs_cmi_genetable$species==input$cmp_species])
      updateSelectizeInput(session, "cmp_gene", choices = species_genes)
      
    })
    
    
    output$cmp_ami_rep1 <- renderImage({
      req(input$cmp_species)
      filename <- normalizePath(file.path('data',"HandE_images",
                                          paste('AMI', '-', '1wk', '-', 'Transplanted1', '.jpeg', sep='')))
      
      list(src = filename, height = 350)
      
      },deleteFile = FALSE
    )

    output$cmp_cmi_rep1 <- renderImage({
        req(input$cmp_species)
        # validate(
        #   need(input$cmp_species != "", "choose a species")
        # )
        filename <- normalizePath(file.path('data',"HandE_images",
                                            paste('CMI', '-', '1wk', '-', 'Transplanted1', '.jpeg', sep='')))
        
        list(src = filename, height = 350)
        
      },deleteFile = FALSE
    )
    
    observe({
      toggle(id="downloadCmpPlots", condition = !is.null(input$cmp_gene))
    })
    
    plot_forDownloadCmpPlots <- function(geneid){
      tryCatch(
        {
          plt <- STutility::FeatureOverlay(ami_vs_cmi_data, 
                                           features = c(geneid), 
                                           pt.size = 1.8,
                                           cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                           dark.theme = F, 
                                           type = "raw",
                                           sampleids = c(1,2,3,4),
                                           ncols = 4,
                                           #min.cutoff = 'q1',
                                           show.sb=F
          )
          return(plt)
        },
        error = function(e) {
          message(paste0("Cannot plot for ",geneid))
          # print(e)
          return(NA)
        },
        warning = function(e){
          message(paste0("warning while plotting for ",geneid))
          # print(e)
          return(NA)
        }
      )
    }
    
    output$downloadCmpPlots <- downloadHandler(
      filename = function() { paste("CMP_featureplots_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".zip", sep = "")},
      content = function(file){
        newTmpDir <- tempfile()
        if(dir.create(newTmpDir)){
          for(gene_ix in seq(1:length(input$cmp_gene))){
            local({
              sel_gene <- input$cmp_gene[gene_ix]
              feature_gene <- paste0(input$cmp_species,'-',sel_gene)
              feature_geneplot <- plot_forDownloadCmpPlots(feature_gene)
              if(any(!is.na(feature_geneplot))){
                plotname <- paste("overlayplot_cmp", sel_gene, sep="_")
                fileName <- file.path(newTmpDir, paste0(plotname, ".png"))
                ggsave(filename = fileName, 
                       plot = feature_geneplot, device = "png",
                       width = 60, height = 20, units = "cm")
              }
              
            })
          }
          zipr(file, files=list.files(newTmpDir, full.names = TRUE))
        }
        unlink(newTmpDir, recursive = TRUE)
      },
      contentType="application/zip"
    )
    
    
    ## cluserplot for selected gene
  
    output$cmp_clusplot <- renderUI({
      req(input$cmp_species)
      req(input$cmp_gene)
      num_plots <- seq(1:length(input$cmp_gene))
      plot_output_list <- lapply(num_plots,
                                 function(g){
                                   plotname <- paste("overlayplot_comp", g,sep="_")
                                   plotOutput(plotname, height = 350)
                                 })
      do.call(tagList, plot_output_list)
    }) 
    
    observe({
      req(input$cmp_species)
      req(input$cmp_gene)
      # deg_table_reset$sel <- input$deg_table_rows_selected
      for(gene_ix in seq(1:length(input$cmp_gene))){
        local({
          icurrent <- input$cmp_gene[gene_ix]
          plotname <- paste("overlayplot_comp", gene_ix, sep="_")
          feature_gene <- paste0(input$cmp_species,'-',icurrent)
          output[[plotname]] <- renderPlot({
            STutility::FeatureOverlay(ami_vs_cmi_data, 
                                      features = c(feature_gene), 
                                      pt.size = 1.8,
                                      cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                      dark.theme = F, 
                                      type = "raw",
                                      sampleids = c(1,2,3,4),
                                      ncols = 4,
                                      #min.cutoff = 'q1',
                                      show.sb=F
            )
          })
        })
      }
    })
    
    output$cmp_clusterplot <- renderImage({

      filename <- normalizePath(file.path('data',"Comparison",
                                          'acute_chronic_clustering.png'))
      
      list(src = filename, height = 300)
      
      },deleteFile = FALSE
      
    )
    
    cmp_tbl <- reactive({
      ami_vs_cmi_degtable %>% dplyr::select(-c('p_val','pct.1','pct.2','gene'))
    })
    
    output$cmp_degtable <- renderDT(server=FALSE, {
      # cmp_tbl <- ami_vs_cmi_degtable %>% dplyr::select(-c('p_val','pct.1','pct.2','gene'))
      datatable(cmp_tbl(),
                extensions=c('Select','Buttons'),
                options = list(scrollX=TRUE,
                               scrollY=TRUE,
                               lengthMenu = c(10),pageLength = 15,
                               paging = TRUE, 
                               select=list(style = 'multi+shift', items = 'row'),
                               searching = TRUE,
                               fixedColumns = TRUE,
                               autoWidth=TRUE, 
                               dom = 'Blfrtip', 
                               escape = c('p_val','pct.1','pct.2','gene'),
                               buttons = c('csv','excel','selectNone')),
                # selection = list(mode = 'multiple', target='row', selected = deg_table_reset$sel),
                selection = 'none',
                callback = JS("table.on( 'select', function ( e, dt, type, ix ) {
                                   var selected = dt.rows({selected: true});
                                   if ( selected.count() > 6 ) {
                                      dt.rows({selected: true}).deselect();
                                   }
                                } );")
      )
    })
    
    output$cmp_table_featureplot <- renderUI({
      req(input$cmp_degtable_rows_selected)
      num_plots <- seq(1:length(input$cmp_degtable_rows_selected))
      plot_output_list <- lapply(num_plots,
                                 function(g){
                                   plotname <- paste("cmp_overlayplot", g,sep="_")
                                   plotOutput(plotname)
                                 })
      do.call(tagList, plot_output_list)
    }) 
    
    validate_gene_cmp <- function(gene_name){
      if(!gene_name%in%row.names(ami_vs_cmi_data@assays$RNA$data)){
          paste0("Selected gene (", gene_name, ") is not present")
      } else{
        NULL
      }
    }
    
    observe({
      req(input$cmp_degtable_rows_selected)
      for(gene_ix in seq(1:length(input$cmp_degtable_rows_selected))){
        local({
          icurrent <- input$cmp_degtable_rows_selected[gene_ix]
          icurrent_species <- cmp_tbl()[icurrent,"species"]
          icurrent_gene <- cmp_tbl()[icurrent, "geneid"]
          plotname <- paste("cmp_overlayplot", gene_ix, sep="_")
          feature_gene <- paste0(icurrent_species,'-',icurrent_gene)
          
          output[[plotname]] <- renderPlot({
            validate(validate_gene_cmp(feature_gene))
            STutility::FeatureOverlay(ami_vs_cmi_data, 
                                      features = c(feature_gene), 
                                      pt.size = 1.8,
                                      cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                      dark.theme = F, 
                                      type = "raw",
                                      sampleids = c(1,2,3,4),
                                      ncols = 4,
                                      #min.cutoff = 'q1',
                                      show.sb=F
            )
          })
        })
      }
    })
    
    observe({
      toggle(id="downloadCmpDegPlots", condition = !is.null(input$cmp_degtable_rows_selected))
    })
    
    plot_forDownloadCmpDegPlots <- function(geneid){
      tryCatch(
        {
          plt <- STutility::FeatureOverlay(ami_vs_cmi_data, 
                                           features = c(geneid), 
                                           pt.size = 1.8,
                                           cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                           dark.theme = F, 
                                           type = "raw",
                                           sampleids = c(1,2,3,4),
                                           ncols = 4,
                                           #min.cutoff = 'q1',
                                           show.sb=F
          )
          return(plt)
        },
        error = function(e) {
          message(paste0("Cannot plot for ",geneid))
          # print(e)
          return(NA)
        },
        warning = function(e){
          message(paste0("warning while plotting for ",geneid))
          # print(e)
          return(NA)
        }
      )
    }
    
    output$downloadCmpDegPlots <- downloadHandler(
      filename = function() { paste("CMP-DEG_featureplots_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".zip", sep = "")},
      content = function(file){
        newTmpDir <- tempfile()
        if(dir.create(newTmpDir)){
          for(gene_ix in seq(1:length(input$cmp_degtable_rows_selected))){
            local({
              icurrent <- input$cmp_degtable_rows_selected[gene_ix]
              icurrent_species <- cmp_tbl()[icurrent,"species"]
              icurrent_gene <- cmp_tbl()[icurrent, "geneid"]
              feature_gene <- paste0(icurrent_species,'-',icurrent_gene)
              feature_geneplot <- plot_forDownloadCmpDegPlots(feature_gene)
              if(any(!is.na(feature_geneplot))){
                plotname <- paste("cmp_overlayplot", icurrent_gene, sep="_")
                fileName <- file.path(newTmpDir, paste0(plotname, ".png"))
                ggsave(filename = fileName, 
                       plot = feature_geneplot, device = "png",
                       width = 60, height = 20, units = "cm")
              }
              
            })
          }
          zipr(file, files=list.files(newTmpDir, full.names = TRUE))
        }
        unlink(newTmpDir, recursive = TRUE)
      },
      contentType="application/zip"
    )
  
  ###########################################################################
  ###################### Cell Cell interaction tab 
  ###########################################################################
    
  cellchat_data <- reactive({
    cellchat_file = normalizePath(
      file.path('data','cell-cell_communications',
                paste0(input$cell_week, '_','LR','_','cellchat','.txt'))
    )
    read.table(cellchat_file, sep="\t", header=T)
  })
 
  observe({
    req(input$cell_week)
    # cellchat_file = normalizePath(
    #   file.path('data',"cell-cell_communications",
    #   paste(input$cell_week, '_', 'LR', '_', 'cellchat', '.txt', sep=''))
    #   )
    # cellchat_data = read.table(cellchat_file, sep="\t", header=T)
    # # ligand_receptors = paste(cellchat_data$ligand, cellchat_data$receptor, sep="-")
    updateSelectizeInput(session, "cell_interactions", choices= cellchat_data()$interaction_name)
  })
    
    
  output$heplot_signalling <- renderImage({
    req(input$cell_week)
    #tissue <- input$tissue_choice
    tissue <- "Transplanted1"
    filename <- normalizePath(file.path('data',"HandE_images",
                                        paste('CMI', '-', input$cell_week, '-', tissue, '.jpeg', sep='')))
    #if(file.exists(filename)){
    #plot_file_path = file.path(getwd(),"data","HandE_images", filename)
    # Return a list containing the filename and alt text
    list(src = filename, height = 280, width=400)
    #}
  },deleteFile = FALSE)
  
  
  output$cell_clusterPlot <- renderImage({
    req(input$cell_week)
    filename <- normalizePath(
        file.path('data',"cell-cell_communications",
                  paste(input$cell_week, '_', 'SeuratClusters', '.jpeg', sep=''))
        # tags$div(img(src = "acute-1.png",height = 600, width = 550,), img(src = "chronic-1.png", height = 600, width = 550,)),
    )
    list(src = filename,  height="auto", width = 700)
  }, deleteFile = FALSE
  )
 
  output$cell_heatmap <- renderImage({
    req(input$cell_week)
    filename <- normalizePath(
      file.path('data',"cell-cell_communications",
                paste(input$cell_week, '_', 'LR', '_', 'heatmap', '.png', sep=''))
    )
    list(src = filename, height="auto", width = 400)
  },deleteFile = FALSE
  )

  
  output$cell_featureplot <- renderUI({
    req(input$cell_interactions)
    num_plots <- seq(1:length(input$cell_interactions))
    plot_output_list <- lapply(num_plots,
                               function(g){
                                 plotname <- paste("cell_overlayplot", g,sep="_")
                                 plotOutput(plotname)
                               })
    do.call(tagList, plot_output_list)
  }) 
  # 
  observe({
    req(input$cell_interactions)
    for(gene_ix in seq(1:length(input$cell_interactions))){
      local({
        interaction <- input$cell_interactions[gene_ix]
        genes <- unlist(strsplit(interaction, "_"))
        plotname <- paste("cell_overlayplot", gene_ix, sep="_")
        feature_genes <- paste('GRCh38',genes, sep='-')
        if(input$cell_week=='1wk'){
          sids= c(1,2)
        }else if (input$cell_week=='4wk'){
          sids = c(3,4)
        }else{
          sids = c(5,6)
        }
        output[[plotname]] <- renderPlot({
          feature_genes <- feature_genes[feature_genes%in%row.names(cmi_data@assays$RNA$data)]
          validate(need(length(feature_genes)>0, 
                        paste0("Selected genes (", feature_genes,") are not expressed")))
          STutility::FeatureOverlay(cmi_data, 
                                    features = c(feature_genes), 
                                    pt.size = 1.8,
                                    cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                    dark.theme = F, 
                                    type = "raw",
                                    sampleids = sids,#unique(cmi_data@meta.data$sampleid),
                                    #ncols = num_cols,
                                    layout.by.feature = FALSE,
                                    #min.cutoff = 'q1',
                                    show.sb=F
          )
        })
      })
    }
  })
  
  observe({
    toggle(id="cell_downloadBtn", condition = !is.null(input$cell_interactions))
  })
  
  plot_forDownloadCellPlots <- function(genes){
    tryCatch(
      {
        feature_genes <- genes[genes%in%row.names(cmi_data@assays$RNA$data)]
        if(input$cell_week=='1wk'){
          sids= c(1,2)
        }else if (input$cell_week=='4wk'){
          sids = c(3,4)
        }else{
          sids = c(5,6)
        }
        plt <- STutility::FeatureOverlay(cmi_data, 
                                         features = c(feature_genes), 
                                         pt.size = 1.8,
                                         cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
                                         dark.theme = F, 
                                         type = "raw",
                                         sampleids = sids, #unique(cmi_data@meta.data$sampleid),
                                         # ncols = 6,
                                         layout.by.feature = FALSE,
                                         #min.cutoff = 'q1',
                                         show.sb=F
        )
        return(plt)
      },
      error = function(e) {
        message(paste0("Cannot plot for ",genes))
        print(e)
        return(NA)
      },
      warning = function(e){
        message(paste0("warning while plotting for ",genes))
        print(e)
        return(NA)
      }
    )
  }
  
  output$cell_downloadBtn <- downloadHandler(
    filename = function() { paste("Cell-Cell_featureplots_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".zip", sep = "")},
    content = function(file){
      newTmpDir <- tempfile()
      # cat("STDERR: Tmpdir: ", newTmpDir,)
      if(dir.create(newTmpDir)){
        for(gene_ix in seq(1:length(input$cell_interactions))){
          local({
            interaction <- input$cell_interactions[gene_ix]
            genes <- unlist(strsplit(interaction, "_"))
            plotname <- paste("cell_overlayplot", gene_ix, sep="_")
            feature_genes <- paste('GRCh38',genes, sep='-')
            feature_geneplot <- plot_forDownloadCellPlots(feature_genes)
            if(any(!is.na(feature_geneplot))){
              plotname <- paste("cell-cell_overlayplot", interaction, sep="-")
              fileName <- file.path(newTmpDir, paste0(plotname, ".png"))
              ggsave(filename = fileName, 
                     plot = feature_geneplot, device = "png",
                     width = 50, height = 30, units = "cm")
            }
            # else{
            #   cat("STDERR: Error in plotting \n")
            # }
          })
        }
        zipr(file, files=list.files(newTmpDir, full.names = TRUE))
      }
      unlink(newTmpDir, recursive = TRUE)
    },
    contentType="application/zip"
  )
}
# Run the application 
shinyApp(ui = ui, server = server)
