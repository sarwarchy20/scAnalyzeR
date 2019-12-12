

library(shiny)
library(shinyjs)
library(ggfortify)
library(ggplot2)
library("gplots")
library(tsne)
library(Rtsne)
#library(pheatmap)
library(plotly)
library("heatmaply")
library(htmlwidgets)
library(dendextend)
library(RColorBrewer)
library(viridis)
library(autoplotly)
library("ggpubr")
library(DT)
#library(Rmagic)
#setwd("E:/Shiny R/sc-vieweR1")
# **************************************** User Interafce of sc-AnalyzeR *********************************************

ui <-shinyUI(
  fluidPage(
    useShinyjs(),
    tags$head(
      tags$style(HTML("@import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
                      h1 {
                      font-style: oblique;
                      font-weight: 1000;
                      line-height:1.2;
                      color: white;
                      text-align: center;
                      font-size: 40px;
                      background-color:green;
                      }
                      h2{
                      font-style: italic;
                      color: white;
                      text-align: center;
                      font-size: 20px;
                      background-color:red;
                      }
                      h3{
                      font-style: italic;
                      color: white;
                      text-align: center;
                      font-size: 20px;
                      background-color:black;
                      }
                      h4{
                      font-style: normal;
                      color: black;
                      text-align: center;
                      font-size: 20px;
                      background-color:yellow;
                      }
                      h5{
                      font-style: italic;
                      color: white;
                      text-align: center;
                      font-size: 20px;
                      background-color:green;
                      }
                      h6{
                      font-style: italic;
                      color: orange;
                      text-align: center;
                      font-size: 20px;
                      background-color:green;
                      }
                      p{
                      font-family: verdana;
                      font-size: 15px;
                      color:#336EFF 
                      }
                      h7 {
                      font-family: verdana;
                      font-size: 15px;
                      color: white 
                      background-color: black
                      }
                        "))
    ),
    
    
    h1("scAnalyzeR"),
    
    tags$style(HTML(".tabbable > .nav > li > a                  {font-weight: bold; background-color: #D4F5F5;  color:black}
                    .tabbable > .nav > li[class=active]    > a {background-color: #00008B; color:white}")),
    #tags$style(HTML(".tabbable > .nav > li[class=inactive]    > a {background-color: black; color:white}")),
    
  #HTML("I like <u>turtles</u>"),
  # Sidebar layout with input and output definitions ----
  tags$head(tags$style(HTML('.progress-bar {background-color: red;}'))),
  sidebarLayout(
    
    tabsetPanel( id = "tabset",
                 
      tabPanel("Instructions", 
               
               img(src = "line_font.png"),
               tags$br(),
               p('  Welcome to scAnalyzeR!'),
               tags$br(),
               "  Please follow the instructions to perform your task successfully:"
               #verbatimTextOutput("inst")
               ),
      tabPanel("Upload Dataset", 
               img(src = "line_font.png"),
               sidebarPanel(id = "data_source", width = 4,
                            tags$style("#data_source{background-color:#54F8EC;}"),
                            radioButtons("species_name", "Species",
                                         choices = c('Human' = "species_human",'Mouse' = "species_mouse"),selected = "species_human"),
                            radioButtons("sor_data", "Source Dataset",
                                         choices = c('Use sample dataset(Human-PBMC)' = "sam_data",'Upload dataset' = "upd_data"),selected = ""),
                                                      
                                          conditionalPanel(
                                                condition = "input.sor_data == 'upd_data'", 
                                                fileInput("file1", "Choose a File(.csv/.txt)"),
                                                multiple = FALSE,
                                                accept = ".csv",
                                                #accept = c("text/csv",
                                                #"text/comma-separated-values,text/plain",
                                                #".csv"),
                                                # Horizontal line ----
                                                tags$hr(),
                                                
                                                # Input: Checkbox if file has header ----
                                                checkboxInput("header", "Header", TRUE),
                                                
                                                # Input: Select separator ----
                                                radioButtons("sep", "Separator",
                                                             choices = c(Comma = ",",
                                                                         Semicolon = ";",
                                                                         Tab = "\t"),
                                                             selected = ","),
                                                
                                                # Input: Select quotes ----
                                                radioButtons("quote", "Quote",
                                                             choices = c(None = "",
                                                                         "Double Quote" = '"',
                                                                         "Single Quote" = "'"),
                                                             selected = '"'),
                                                # Horizontal line ----
                                                #tags$hr(),
                                                #checkboxInput('upload_summary',p('Uploaded Data Summary')),
                                                tags$hr()
                                                )
                                                      #uiOutput("source_UI"),
                                                   
                                                      # Input: Select number of rows to display ----
                                                       #, checkboxInput('show',p('Show Loaded Dataset')),
                                                              # radioButtons("disp", "Display",
                                                                            #choices = c(Head = "head",
                                                                            #All = "all"),
                                                                            #selected = "head")
      
      )# end of sidebarPanel for Uplading Dataset
      , column(width =6,  br(),h2(htmlOutput("text"))),
      DT::dataTableOutput("testdata")
      ), # end of tabPanel for Uplading Dataset
      tabPanel("Pre-processing",
               img(src = "line_font.png"),
               sidebarPanel(id= "prep_pro" ,width = 4,
                            tags$style("#prep_pro{background-color:#EED4F9;}"),
                                        h4("Create Dataset "),
                                        " Keep all genes expressed in >= MC cells and Keep all cells with at
                                        least MFs detected genes",
                                        numericInput("cells", "Minimum no. of cells(MC):", 3, min = 1, max = Inf),
                                        numericInput("genes", "Minimum no. of features(MFs):", 200, min = 1, max = Inf),
                            
                                         tags$head(
                                          tags$style(HTML('#ok{background-color:#00008B;color:white;}'))
                                          ),
                            
                                        actionButton('ok', 'Submit'),
                                        tags$hr(),
                                        radioButtons("disp1", "Display Genes",
                                                     choices = c(Mitochondrial= "mgene",
                                                                 None = "none"),
                                                     selected = "none"),
                                        checkboxInput('mdatack',p('Plot Meta Data')),
                                        checkboxInput('mitogene',p('Download Mitochondrial Genes(as a .csv file)')),
                                        conditionalPanel('input.mitogene==1',
                                                         downloadButton("mtData", "Download")),
                                        tags$hr(),
                                        checkboxInput('gene_summary',p('Plot Gene Summary')),
                                        tags$hr(),
                                        checkboxInput('fth',h6("  Filtering ")),
                                        conditionalPanel('input.fth==1',
                                        " filter cells that have unique feature counts over nFeature_RNA(>) and 
                                        less than nFeature_RNA(<)",
                                        numericInput("lt", "nFeature_RNA(>):", 200, min = -Inf, max = Inf),
                                        numericInput("ut", "nFeature_RNA(<):", 2500, min = 1, max = Inf),
                                        numericInput("pmt", "MT gene % :", 2, min = 0 , max = Inf),
                                        tags$head(
                                          tags$style(HTML('#th{background-color:#8B0000;color:white;}'))
                                        ),
                                        actionButton('th', 'Submit')
                                        )
                                            
                                        
      ) # end of sidebarPanel for Filtering
      ,column(width =6,br(), br(),h4(htmlOutput("text1")))
      ,column(width =6,br(), br(),h6(htmlOutput("thtext")))
      #,column(width =12,br(), br(),plotOutput("fdata_umi_gene"))
      ,column(width =6,br(), br(),DT::dataTableOutput("fgene"))
      ,column(width =12,br(), br(),plotOutput("mdataplot"))
      ,column(width =12,br(), br(),
              splitLayout(
                style = "border: 1px solid silver;",
                cellWidths = 300,
                cellArgs = list(style = "padding: 6px"),
                plotOutput("umi_mito"),
                plotOutput("umi_gene"),
                plotOutput("fdata_umi_gene")
              ))
      ),# end of tabPanel for Filtering
      
      tabPanel("Normalization",
               img(src = "line_font.png"),
               sidebarPanel(id  = "nor_slid" , width = 4,
                            tags$style("#nor_slid{background-color:#EFF38F;}"),
                             h3("Set Normalization Parameters"),
                             selectizeInput("nor_method", "Normalization method", 
                                                           c('Log Normalization'="LogNormalize",
                                                             'Centered log ratio transformation(CLR)' = "CLR",
                                                             'Relative counts(log-transformation is NOT applied)' = "RC"),
                                                                selected = 'LogNormalize'),
                             numericInput("nor_factor", "Scaling factor",10000,min = 1, max = Inf),
                             tags$head(
                               tags$style(HTML('#nor_ok{background-color:#00008B;color:white;}'))
                             ),
                            actionButton('nor_ok', 'Submit'),
                tags$hr(),  
                  h5("Find Highly Variable Features"),
                  numericInput("top_vgenes", "Features", 1000, min = 1, max = Inf),
                tags$head(
                  tags$style(HTML('#hvg_ok{background-color:#8B0000;color:white;}'))
                ),
                  actionButton('hvg_ok', 'Submit'),
                tags$hr(),  
                checkboxInput('hfp_disp',p('Variable Features Plot')),
                conditionalPanel('input.hfp_disp==1',"",
                                 numericInput("top_hfp", "Top Features:", 10, min = 1, max = Inf)
                ),
                tags$hr(),
                  checkboxInput('hvg_disp',p('Show Highly Variable Genes')),
                  checkboxInput('high_hvg',p('Download Highly Variable Genes')),
                  conditionalPanel('input.high_hvg==1',
                                downloadButton("hvg_download", "Download"))
    ), # # end of sidebarPanel for Normalization
    column(width =6,br(), br(),h6(htmlOutput("nor_text"))),
    column(width =6,br(), br(), verbatimTextOutput("hvg_list")),
    column(width = 6, br(), br(),DT::dataTableOutput("hvg_table")),
    column(width =12,br(), br(),plotOutput("hfp_plot"))
    ),# end of tabPanel for Normalization
      
    
    tabPanel("Dimensionality Reduction", 
             img(src = "line_font.png"),
             sidebarPanel(id = "dim_slid", width = 4, 
                          tags$style("#dim_slid{background-color:#A9F19C;}"),
                                 radioButtons("pc_genes", "Gene use",
                                              choices = c('High Variable Genes' = "pc_hvg",'All Genes' = "pc_all"),selected = "pc_hvg"),
                                 column(width = 8, numericInput("pc_cells", "Number of PCs:", 20, min = 1, max = Inf)),
                                 
                                 #checkboxInput('pc_impute',h3('Apply imputation')),
                                 
                                  
                                 tags$head(
                                   tags$style(HTML('#com_pca{background-color:#8B0000;color:white;}'))
                                 ),
                                 actionButton('com_pca', 'Compute PCA'),
                                               
                                  tags$hr(),
                                 checkboxInput('ck_pca',p('PCA Plot')),
                                               #tags$hr(),
                                 checkboxInput('pc_print',p('Print PCA')),
                                              # tags$hr(),
                                 checkboxInput('pc_elbow',p('Elbow Plot')),
                                               # tags$hr(),
                                 checkboxInput('pc_jac_plot',p('Jack Straw Plot(process can take a long time for big datasets)')),
                                                #tags$hr(),
                                 checkboxInput('pc_hmap',p('PC Heatmap')),
                                 conditionalPanel('input.pc_hmap==1', 
                                                  "Please set the PCs to plot the Heatmap:",
                                                  br(),
                                                  column(width=4,numericInput("pc_low", "Lower limit",1,min = 1, max = Inf)),
                                                  column(width=4,numericInput("pc_upper", "Upper limit",1,min = 1, max = Inf)),
                                                  br(),br(),br(),br(),br(),
                                                  column(width = 8, numericInput("pc_hmap_cells", "Number of Cells:", 500, min = 2, max = Inf)),
                                                  br(), br(), 
                                                  actionButton('pc_hmap_ok', 'Show plot')
                                ),
                             checkboxInput('ck_tsne',p('t-SNE Plot')),
                             conditionalPanel('input.ck_tsne == 1',
                                              radioButtons("tsn_genes", " ",
                                                           choices = c('High Variable Genes' = "tsn_hvg",'All Genes' = "tsn_all"),selected = "tsn_hvg")
                               )
                                 
    ),# end of sidebarPanel for PCA 
    column(width =6,br(), br(),h4(htmlOutput("pc_text"))),
    column(width =8, br(), br(), plotOutput("sc_pca_plot")),
    column(width = 6, br(),verbatimTextOutput("pca_list")),
    column(width = 8, br(), plotOutput("pc_jac_show")),
    column(width = 6, br(), plotOutput("pc_elbow_plot")),
    column(width =10, br(), br(), plotOutput("pc_hplot")),
    column(width = 8, br(), plotOutput("tsnep"))
    ),# end of tabPanel for PCA
      tabPanel("Clustering", 
               img(src = "line_font.png"),
               sidebarPanel(id  = "cls_slid", width = 4,
                            tags$style("#cls_slid{background-color:#A3CFF5;}"),
                            "Please select the PCs for clustering:",
                                                             br(),
                                                            column(width=4,numericInput("clus_low", "Lower limit",1,min = 1, max = Inf)),
                                                            column(width=4,numericInput("clus_upper", "Upper limit",2,min = 1, max = Inf)),
                                                            column(width =8, numericInput("clus_restn","Resolution",0.8,min=0,max = Inf)),
                                                            column(width =8, numericInput("clus_itr","Iterations",10,min=1,max = Inf)),
							                                              selectInput('clus_algo','Select Algorithm',choices = c("Louvain Algorithm(LA)", 
							                                                                                                     "Multilevel LA", "SLM algorithm",
							                                                                                                     "Leiden algorithm"),selected = 'Louvain Algorithm(LA)'),
                                                            tags$head(
                                                            tags$style(HTML('#clus_ok{background-color:black;color:white;}'))
                                                            ),
                                                            actionButton('clus_ok', 'Do cluster'),
                                                            
                                                            #tags$hr(),
                                                            #checkboxInput('clus_umap',p('UMAP plot')),
                                          
                                                            tags$hr(),
                                                            checkboxInput('clus_tsne',p('t-SNE plot')),
                                          
                                                            tags$hr(),
                                                            checkboxInput('clus_barplot',p('Show Cluster Bar Plot'))
                              #selectInput('Clusm','Clustering Method',choices = c("K-Means", "Hierarchical Clustering"),selected = 'K-Means')
      ), # end of sidebarPanel for Clustering
      column(width =6,br(), br(),h6(htmlOutput("clus_text"))),
      column(width =6, br(), br(), plotOutput("clus_tsne_plot")),
      #column(width =6, br(), br(), plotOutput("clus_umap_plot")),
      column(width =12, br(), br(), plotOutput("clus_barplot_show"))
      ), # end of tabpanel of Clustering
    
    
      tabPanel("Differential Expression",
               img(src = "line_font.png"),
               mainPanel(
                                                 
                                          uiOutput("DE_clus_UI")
                                          
                                          
                    )  #  DE MainPanel
      ), # # end of tabpanel of Differential Expression 
    
    tabPanel("Plots", 
             img(src = "line_font.png"),
             mainPanel(
    tabsetPanel(id = "plots_tabset",
    
    tabPanel("Violin Plot", 
             img(src = "line_font1.png"),
             sidebarPanel(id = "vplot_slid" , width = 6,
                          tags$style("#vplot_slid{background-color:#A6F0F4;}"),
                          textAreaInput("vln_genes","Write Gene List", 
                                                                   placeholder = "Please write gene sybmols: one gene per line.",
                                                                   width = "200px", height = "250px")),
             
             column(width = 12,br(), br(), plotOutput("vln_plot"))
    ), # end taPanel for violin plot
    
    # ------Dynamic Heatmap
    tabPanel("Dynamic Heatmap",
             img(src = "line_font1.png"),
             sidebarPanel(id  = "dhm_slid", width = 5, 
                          tags$style("#dhm_slid{background-color:#DEF3C9;}"),
                                            #h4("Choose Your Input Option!"),
                                            # Input: Select the options----
                                            radioButtons("dhmap_input_opt", h3("Choose the Input Option"),
                                                         choices = c('Input Gene list as a text file' = "dhmap_file_inp",
                                                                     'Write Gene List' = "dhmap_write_inp"),
                                                         selected = "dhmap_file_inp"),
                                            
                                            uiOutput("dhmap_input_UI"),
                                            #submitButton("Submit", icon("refresh")),
                                            #verbatimTextOutput("value"),
                                            #textInput("text", label = h3("Text input"), value = "Enter text...")
                                            #sliderInput("controller", "Controller", 1, 4341, 10)
                                            
                                            br(),h5('Row dendrogram'),
                                            column(width=6,selectizeInput("distFun_row", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                            column(width=6,selectizeInput("hclustFun_row", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                            column(width=12,sliderInput("r", "Number of Clusters", min = 1, max = 15, value = 2)) ,
                                            
                                            br(),br(),br(),br(),br(),br(),br(),br(),br(),hr(),h3('Column dendrogram'),
                                            column(width=6,selectizeInput("distFun_col", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                            column(width=6,selectizeInput("hclustFun_col", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                            column(width=12,sliderInput("c", "Number of Clusters", min = 1, max = 15, value = 2)),
                                            
                                            br(),br(),br(),br(),br(),br(),br(),br(),br(),hr(),   #h4('Color Manipulation'),
                                            h6(uiOutput('colUI')),
                                            #column(3,checkboxInput('showColor','Color'))
                                            #br(),#
                                            hr(), 
                                            h2('Dendrogram Manipulation'),
                                            selectInput('dtype','Dendrogram Type',choices = c("both", "row", "column", "none"),selected = 'both'),
                                            sliderInput('branches_lwd','Dendrogram Branch Width',value = 0.6,min=0,max=5,step = 0.1)
                                            
    ) # end of sidebarPanel for Dynamic Heatmap
    
    , column(width =7,  br(),plotlyOutput("dh_plot", height='600px'))
    ),   # end of the Dynamic heatmap 
    
    
    # --------------------------------------------------- Static Heatmap --------------------------------------------------
    
    tabPanel("Static Heatamp", 
             img(src = "line_font1.png"),
             sidebarPanel(id  = "stmp_slid" ,width = 5,
                          tags$style("#stmp_slid{background-color:#EAC9F3;}"),
                                            #h4("Choose Your Input Optio:!"),
                                            radioButtons("stmap_input_opt", h4("Choose Your Input Option"),
                                                         choices = c('Input Gene list as a text file' = "stmap_file_inp",
                                                                     'Write Gene List' = "stmap_write_inp"),
                                                         selected = "stmap_file_inp"),
                                            
                                            uiOutput("stmap_input_UI"),
                                            radioButtons("stmap_optn", "Select Plot",
                                                         choices = c('Heatmap with clusters' = "stmap_optn_1",'Heatmap without Clusters' = "stmap_optn_2"),
                                                         selected = "stmap_optn_1")
                                            
                                            #numericInput("stmap_cells", "Cells",200,min = 2, max = Inf)
                                            
                                            
                                            #fileInput("file", "Input Gene list as a text file:")
                                            
                                            
    ), # end of sidebarPanel for Static Heatmap
    
    
    column(width =8,  br(),br(), plotOutput("sta_hmap", height='800px'))
    
    ), # End of tabPanel for Static Heatmap
    
    #------------------- Dynamic PCA -------------------
    
    #tabPanel("Data Summary", br(), br(), DT::dataTableOutput("mytable")),
    tabPanel("Dynamic PCA",
             img(src = "line_font1.png"),
             sidebarPanel(id  = "dypca_slid", width = 4, 
                          tags$style("#dypca_slid{background-color:#8FF065;}"),
                                        #h4("Choose Your Input Option!"),
                                        fileInput("dpca_file", "Input Gene list as a text file (one gene per line)")
                                        
    ),
    column(width = 12,br(), br(), plotlyOutput("dpcaplot", height='600px'))
    
    ), # End of tabPanel for Dynamic PCA
    
    
    
    #tabPanel("T-SNE Plot",
            # img(src = "line_font1.png"),
            # sidebarPanel(width = 4, 
                                        #h4("Choose Your Input Option!"),
                                       # fileInput("tsne_file", "Input Gene list as a text file:")
                                        
    #),
    #column(width = 12,br(), br(), plotOutput("tsnep", height='600px'))
    
    #),
    
    tabPanel("Correlation Plot",
             img(src = "line_font1.png"),
             sidebarPanel(id = "corr_slid", width = 4,
                          tags$style("#corr_slid{background-color:#DCDEF0;}"),
                                              textInput("gene1","Gene1(X)", value = "GAPDH"),
                                              textInput("gene2","Gene1(Y)", value = "TP53"),
                                              selectInput('corm','Correlation Method',choices = c("pearson", "spearman", "kendall"),selected = 'pearson'),
                                              selectInput('regline','Regression Line',choices = c('adding linear regression line'="reg.line", 
                                                                                                  'adding local regression fitting'="loess", "none"),
                                                          selected = 'reg.line'),
                                              checkboxInput('cor_log',p('Apply Log2'))
                                              #, checkboxInput('zeroexp',p('Excluding ZERO Expression value'))
                                             # , checkboxInput('cor_impute',p('Impute'))
                                              
    ) # end of sidebarPanel for Correlation plot
    ,column(width =6, br(), br(), plotOutput("corplot"))
    ) # End of the Correlation plot
    
     ) #End of tabsetPanel for plots
      )
    ),
    
    
    
    # -------------------------------------------------- Start Pathways Analysis------------------------------------------
    tabPanel("Pathways Analysis",
             img(src = "line_font.png"),
             sidebarPanel(id = "path_slid" ,width = 4,
                          tags$style("#path_slid{background-color:#EEA6F4;}"),
                                              radioButtons("path_input_opt", h4("Select DE Genes"),
                                                           
                                                           choices = c('Cluster Specific DE Genes' = "path_clus",
                                                                       'Cluster(s) vs Cluster(s) DE Genes' = "path_clus_vs_clus"),
                                                           selected = "path_clus"),
                                              
                                              #choices = c('All DE Genes' = "path_all_de",
                                              #           'Cluster Specific DE Genes' = "path_clus",
                                              #           'Clusters vs Clusters DE Genes' = "path_clus_vs_clus"),
                                              #selected = "path_all_de"),
                                              
                                              #h2('Select Pathway'),
                                              selectInput('path_source',h3('Select Pathway Source'),choices = c('KEGG Pathways' = "path_KEGG",
                                                                                                                'Gene Ontology' =  "path_GO"),
                                                          selected = "Path_KEGG"),
                                              #br(), br(), 
                                              tags$head(
                                                tags$style(HTML('#path_submit{background-color:red;color:white;}'))
                                              ),
                                              actionButton('path_submit', 'Compute Pathways'),
                                              br(),br(),
                                              
                                              uiOutput("path_download_UI"),
                                              br(),
                                              br(),
                                              #conditionalPanel('input.path_submit==1',
                                                               
                                                               
                                                               radioButtons("path_select", p("Select Pathways"),
                                                                            choices = c('Positive Pathways' = "path_pos",
                                                                                        'Negative pathways' = "path_neg"),
                                                                            selected = "path_pos"),
                                              
                                                               numericInput("path_qval","Cutoff q-value(<= cutoff):", 0.1, min = 0,max = Inf),                 
                                                                
                                                              uiOutput("path_ftr_download_UI"),
                                              
                                                               checkboxInput('path_top_plot',h4('Show Pathways(plot)')),
                                                               conditionalPanel('input.path_top_plot==1',
                                                                                numericInput("path_top","Top Pathways:", 10, min = 1,max = Inf)
                                                               ),
                                                               
                                                               checkboxInput('path_list_show',h4('Show Pathways list')),
                                              
                                                               checkboxInput("path_genesets",p('Show Gene sets:')),
                                                                      conditionalPanel('input.path_genesets==1',
                                                                                       numericInput("path_no","Pathway No.:", 1, min = 1,max = Inf),
                                                                                       radioButtons("path_hmap_genes", p("Select Gene Set"),
                                                                                                    choices = c('Significant Genes Only' = "path_hmap_genes_sig",
                                                                                                                'All Genes ' = "path_hmap_genes_all"),
                                                                                                    selected = "path_hmap_genes_sig"),
                                                                                       
                                                                                       checkboxInput("path_hmap",p('Show Gene set Heatmap:'))
                                                                                       
                                                                      )
                                                                
                                                               
                                                               
                                                               #downloadButton("path_all_download", "Save Pathways")
                                              #)
                                              
                                              
    ), # end of sidebarPanel for Pathways
    column(width =6,br(), br(),DT::dataTableOutput("path_list_table")),
    column(width =10, br(), br(), plotOutput("path_show")),
    column(width = 8, br(),verbatimTextOutput("path_gene_list")),
    column(width =10, br(), br(), plotOutput("path_hmap_show"))
    ),
   
    # ------------ End Pathways Analysis
    
    # ----------------------------------------------------------- Trajectory Analysis --------------------------
      tabPanel("Trajectory Analysis",
               img(src = "line_font.png"),
               sidebarPanel(id = "tray" , width = 4, 
                                                  tags$style("#tray{background-color:#A3F5D1;}"),
                                        #h4("Choose Your Input Option!"),
                                        h2("Select Analysis Options"),
                                        radioButtons("tray_genes", "Gene use",
                                                     choices = c('High Variable Genes Only' = "tray_hvg",'All Genes' = "tray_all"),
                                                     selected = "tray_hvg"),
                                        
                                        checkboxInput("tray_plot_set",p('Cell trajectory plot')),
                                        conditionalPanel('input.tray_plot_set==1',
                                                         
                                                         radioButtons("tray_plot_set_opts", p("Select plot style"),
                                                                      choices = c('Clusters' = "tray_clst_ck",
                                                                                  'Pseudotime' = "pseudotime_ck",
                                                                                  'State' = 'tray_stat_ck'),
                                                                      selected = "tray_clst_ck")
                                                         
                                        ),
                                        
                                        checkboxInput("tray_pseu_set",p('Pseudotemporal plot')),
                                                     
                                                    
                                        conditionalPanel('input.tray_pseu_set==1',
                                                         radioButtons("tray_input_opt", h4("Choose the input option"),
                                                                      choices = c('Input Gene list as a text file' = "tray_file_inp",
                                                                                  'Write Gene List' = "tray_write_inp"),
                                                                      selected = "tray_file_inp"),
                                                         
                                                         uiOutput("tray_input_UI"),
                                                         
                                                         numericInput("tray_pseu_clst","Clusters:", 5, min = 1,max = Inf),
                                                         tags$head(
                                                           tags$style(HTML('#tray_submit{background-color:red;color:white;}'))
                                                         ),
                                                         actionButton('tray_submit', 'OK')
                                                         
                                                         
                                                         
                                        )
                                        
    ),
    column(width = 8,br(), br(), plotOutput("tray_plot", height='500px')),
    column(width = 8,br(), br(), 
           plotOutput("tray_pseu_plot")
           )
    
    ) # End of tabPanel for Trajectory Analysis
        ),
    
    mainPanel(
      #tabsetPanel(type = "tabs",
       # tabPanel("Instructions", verbatimTextOutput("inst")), 
       # tabPanel("Dynamic Heatmap", br(), br(), plotlyOutput("plot", height='600px')),
        #tabPanel("Static Heatamp", br(), br(), plotOutput("plot1")), 
       # tabPanel("Data Summary", br(), br(), tableOutput("table")),
        #tabPanel("PCA Plot", br(), br(), plotlyOutput("pcaplot", height='600px'))
      # )
      
      ) # end main panel
  )
 )    
)
