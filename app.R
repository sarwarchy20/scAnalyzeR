
#setwd("E:/seurat/Seurat_3 doc/")
#use_python("C:/Program Files/Python37/python.exe")


source("install_dependencies_scAnalyzer.R")
source("load_packages.R")
source("ui.R")

#use_python("C:/Users/GSC/AppData/Local/conda/conda/envs/r-reticulate/python.exe")
# ************************************************ Global section ***************************************************

#library(survival)
#library("leiden")
#library("igraph")
#new2_data<-read.csv(unz("pbmc4k.zip", "pbmc4k.csv")) 
#new2_data<-read.csv("summing_data_1.csv") 
#colnames(new2_data) <- c("Gene_Symbol", colnames(new2_data)[-1])
#row.names(new2_data) <-new2_data$Gene_Symbol 


#hg<-new2_data[filenames,]
#By default, Shiny limits file uploads to 5MB per file.
# Increase the uploading file size up to 30 MB
#options(shiny.maxRequestSize = 30*1024^2)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10240*1024^2)
#That way the max limit is 10GB when the app is run locally and 5MB when run from the server.
f_total_col = 0 # total number of column of filtered dataset
mito.genes<-c()
path_title<- NULL
new2_data<- NULL

#-----------------------End of Global Section----------------------------------------------

#******************************************* Start Server Section ************************************************
server <- function(input, output,session) {
  #--------------------------------- Instructions -------------------------------------
  # output$inst<-renderText({
  # paste("Please follow the instructions to perform your task successfully:")
  #print(input$Tgenelist)
  
  #})
  
  # --------------------Dataset uploading--------------------------------------------------------------------------
  new2_data<-reactive({
    #matrix(1, ncol = 1, nrow = 1)
    req(input$sor_data)
    if(input$sor_data == "upd_data")
    {
      req(input$file1)
      # withProgress(message = 'Filtering', value = 0, {
      
      #n <- 10
      
      #for (i in 1:n) {
      df <- read.csv(input$file1$datapath,
                     header = input$header,
                     sep = input$sep,
                     quote = input$quote,
                     row.names = 1)
      #colnames(df) <- c("Gene_Symbol",colnames(df)[-1])
      #row.names(df) <-df$Gene_Symbol
      # Increment the progress bar, and update the detail text.
      #incProgress(1/n, detail = paste("Completed ", i*10,"%"))
      
      # Pause for 0.1 seconds to simulate a long computation.
      #Sys.sleep(0)
      #  }    
      #})
      df
    }
    else if(input$sor_data == "sam_data")
    {
      
      #inp_data = read.csv(file = ".csv", # load genes with Z-scores
      # header = TRUE,
      # row.names = 1)
      #setwd("E:/seurat/Seurat_3 doc/www/")
      withProgress(message = 'Loading...', value = 0, {
        
        n <- 1
        
        for (i in 1:n) {
          pbmc.data_2 <- read.csv("sample_data/mouse_ovarian.csv", header = TRUE,
                                  row.names = 1)
          
          # Increment the progress bar, and update the detail text.
          incProgress(10, detail = paste("Completed ", i*100,"%"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(4)
        }
      })
      pbmc.data_2
      
    }
    
    else if(input$sor_data == "data_10x")
    {
      req(input$file_10x)
      #print(input$file_10x$datapath)
      #print(input$file_10x$datapath[1])
      #print(str_sub(input$file_10x$datapath[1], end=-6))
      
      #temp1<- read.delim("filtered_gene_bc_matrices/hg19/matrix.mtx")
      
      #write.table(temp1,"filtered_gene_bc_matrices/hg19/ttt_3.tsv",row.names =F,quote = F)
      
      
      p1<- str_sub(input$file_10x$datapath[1], end=-6)
      p11<- paste(p1, sep = "", "barcodes.tsv")
      file.rename(from =input$file_10x$datapath[1], to = p11)
      
      q1<- str_sub(input$file_10x$datapath[2], end=-6)
      q11<- paste(q1, sep = "", "genes.tsv")
      file.rename(from =input$file_10x$datapath[2], to = q11)
      
      r1<- str_sub(input$file_10x$datapath[3], end=-6)
      r11<- paste(r1, sep = "", "matrix.mtx")
      file.rename(from =input$file_10x$datapath[3], to = r11)
      
      pbmc.data_10x<- Read10X(data.dir = str_sub(input$file_10x$datapath[1], end=-6))
      data.matrix(pbmc.data_10x)
      
    }
    
    
    
  })
  
  #..................... Test Uploaded dataset ........................
  output$text <- renderUI({
    if(!is.null(new2_data())){
      str1<-paste("Summary of The Uploaded Dataset:")
      str2 <- paste("Total number of features(genes): ", nrow(new2_data()))
      str3 <- paste("Total number of samples(cells): ", ncol(new2_data()))
      HTML(paste(str1, str2 ,str3, sep = '<br/>'))
    }
    
  })
  
  #observeEvent(input$upload_summary,{
  # show("text")
  #})
  # Display uploaded data
  #output$testdata = DT::renderDataTable({
  #req(input$file1)
  #if(input$show == 1){
  # hide("text")
  #if(input$disp == "head") {
  # return(head(new2_data()[1:5,2:3]))
  # }
  #else {
  #  return(new2_data()[1:5,2:3])
  #  }
  #}
  #})
  # ------------------------------------------------Pre-processing--------------------------------------------
  ################### Create Dataset
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$ok, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$ok
  })
  
  #observeEvent(input$tabset, {
  #v$doPlot <- FALSE
  #})
  
  pbmc<- reactive({
    #if(input$ok == TRUE){
    if (v$doPlot == FALSE) return()
    #withProgress(message = 'Filtering', value = 0, {
    
    # n <- 100
    #
    # for (i in 1:n) {
    isolate({
      total_col<-ncol(new2_data())
      pbmct<- CreateSeuratObject(counts = new2_data(), 
                                 min.cells = input$cells, min.feature = input$genes, 
                                 project = "Cells")
      # Increment the progress bar, and update the detail text.
      #incProgress(1/n, detail = paste("Completed ", i,"%"))
      
      # Pause for 0.1 seconds to simulate a long computation.
      # Sys.sleep(0)
      # }    
      #})
      tcol<- ncol(GetAssayData(object = pbmct)) # test seurat object , return if object conatins zero columns(cells)
      if(tcol<=0)
      {
        pbmct
        
      }
      else
      {
        pattern_code<- switch(input$species_name, species_human = '^MT-', species_mouse = '^mt-')
        pbmct[["percent.mt"]] <- PercentageFeatureSet(object = pbmct, pattern = pattern_code)
        pbmct
      }
    })
  })
  # ............. summary of dataset........
  output$text1<-renderUI({
    
    if(!is.null(pbmc())){ 
      str1<-paste("Dataset Created with:")
      str2 <- paste("Total number of features(genes): ", nrow(GetAssayData(object = pbmc())))
      str3 <- paste("Total number of samples(cells): ", f_total_col=ncol(GetAssayData(object = pbmc())))
      HTML(paste(str1, str2 ,str3, sep = '<br/>'))
    }
  })
  # ....................Show Filtered Gene(mito genes)..............
  output$fgene = DT::renderDataTable({
    #if(input$ok == TRUE){
    if(input$disp1 == "mgene") {
      if(input$species_name == "species_human"){
        mito.genes<- grep(pattern = "^MT-", x = rownames(GetAssayData(object = pbmc())), value = TRUE)
      }
      else
      {
        mito.genes<- grep(pattern = "^mt-", x = rownames(GetAssayData(object = pbmc())), value = TRUE)
      }
      
      mt_data<- data.matrix(GetAssayData(object = pbmc()[mito.genes,]))
      
      mtd2<- data.frame(rowSums(mt_data))
      
      mtd2<- mtd2 %>% 
        dplyr:: rename(
          Total_Count = rowSums.mt_data.
        )
      mtd2
      
    }
    else 
      return()
    #}
  })
  # ............plot meta data
  output$mdataplot<- renderPlot({
    if(input$mdatack == 1){
      vln_plot<- VlnPlot(object = pbmc(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      vln_plot + ggtitle("Before Filtering") +
        theme(plot.title = element_text(hjust = 0.5, size = 25,face = "bold"))
    }
    else 
      return ()
  })
  # ...............Download Mitochondrial  genes data set.........
  output$mtData<-downloadHandler(
    
    if(input$mitogene == FALSE) return(),
    
    filename = function() {
      paste("mitochondrial_genes", ".csv", sep ="")
    },
    content = function(file) {
      if(input$species_name == "species_human"){
        mito.genes<- grep(pattern = "^MT-", x = rownames(GetAssayData(object = pbmc())), value = TRUE)
      }
      else if(input$species_name == "species_mouse")
      {
        mito.genes<- grep(pattern = "^mt-", x = rownames(GetAssayData(object = pbmc())), value = TRUE)
      }
      md_data<- GetAssayData(object = pbmc())[mito.genes,]
      write.csv(as.matrix(md_data), file, row.names = TRUE, na = "")
    }
  )
  # plot gene summary
  output$umi_mito<-renderPlot({
    if(input$gene_summary == FALSE) return()
    FeatureScatter(object = pbmc(), feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("Before filtering")
  })
  
  output$umi_gene<-renderPlot({
    if(input$gene_summary == FALSE) return()
    FeatureScatter(object = pbmc(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("Before filtering")
  })
  # ------------------------------------------------ Filtering -------------------
  pbmcf<-reactive({
    if(input$th == FALSE) return()
    
    isolate({
      tpbmc<-pbmc()
      
      tpbmc[["MTg"]]<- input$pmt
      tpbmc[["ltv"]]<- input$lt
      tpbmc[["utv"]]<- input$ut
      te<-0.5
      
      print(head(x = tpbmc@meta.data,5))
      print(input$pmt)
      pbmc1<-subset(x = tpbmc,(subset =  nFeature_RNA >ltv  & nFeature_RNA < utv  & percent.mt < MTg & ltv< Inf))
      
      pbmc1
    })
  })
  
  output$thtext<-renderUI({
    
    if(is.null(pbmcf())) return()
    else if(input$th == FALSE) return()
    str1<-paste("Filtered Dataset:")
    str2 <- paste("Total number of features(genes): ", nrow(GetAssayData(object = pbmcf())))
    str3 <- paste("Total number of samples(cells): ", ncol(GetAssayData(object = pbmcf())))
    
    HTML(paste(str1, str2 ,str3, sep = '<br/>'))
    
  })
  
  # ............plot meta data after filtering
  output$fth_mdataplot<- renderPlot({
    if(input$fth_mdatack == 1){
      vln_plot<- VlnPlot(object = pbmcf(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      vln_plot + ggtitle("After Filtering") +
        theme(plot.title = element_text(hjust = 0.5, size = 25,face = "bold"))
    }
    else 
      return ()
  })
  
  # plot gene summary after filtering
  output$fth_umi_mito<-renderPlot({
    if(input$fth_gene_summary == FALSE) return()
    FeatureScatter(object = pbmcf(), feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("After filtering")
  })
  
  output$fth_umi_gene<-renderPlot({
    if(input$fth_gene_summary == FALSE) return()
    FeatureScatter(object = pbmcf(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("After filtering")
  })
  
  # ------------------------------------------------------------- Normalization -----------------------------------
  
  nor_pbmc<-reactive({
    #if(is.null(pbmcf())) return()
    if(input$nor_ok ==FALSE) return()
    isolate({
      nor_data<-NormalizeData(object = pbmcf(), normalization.method = input$nor_method, 
                              scale.factor = input$nor_factor,display.progress = TRUE)
      nor_data
      #sc_pbmc <- ScaleData(object = nor_data, vars.to.regress = c("nUMI", "percent.mito"))
      #sc_pbmc
    })
  })
  
  #........... Normalized data displsy.........
  output$nor_text<-renderUI({
    if(is.null(nor_pbmc())) return()
    else if(input$nor_ok == FALSE) return()
    str1<-paste("Normalized Dataset:")
    str2 <- paste("Total number of features(genes): ", nrow(GetAssayData(object = nor_pbmc())))
    str3 <- paste("Total number of samples(cells): ", ncol(GetAssayData(object = nor_pbmc())))
    HTML(paste(str1, str2 ,str3, sep = '<br/>'))
  })
  #.................... Highly Variable Genes.....................................................................................
  hv_pbmc<-reactive({
    if(is.null(nor_pbmc())) return()
    else if(input$hvg_ok == FALSE) return()
    isolate({
      # temp_hvg<-FindVariableGenes(object = nor_pbmc(), mean.function = ExpMean, dispersion.function = LogVMR, 
      #x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,top.genes = input$top_vgenes)
      temp_hvg<-FindVariableFeatures(object = nor_pbmc(), selection.method = "vst", 
                                     nfeatures = input$top_vgenes)
      temp_hvg
    })
  })
  
  # Find Highly Variable Features
  output$hvg_list<-renderText({
    if(is.null(hv_pbmc())) return()
    else if(input$hvg_ok == FALSE) return()
    #temp<-hv_pbmc()@var.genes
    ##temp
    print(paste("Highly Variable Features Detected: ",length(x= VariableFeatures(object = hv_pbmc()))))
    
  })
  
  # Variable Features Plot
  output$hfp_plot<-renderPlot({
    if(is.null(hv_pbmc())) return()
    else if(input$hfp_disp == FALSE) return()
    top10 <- head(x = VariableFeatures(object = hv_pbmc()),input$top_hfp)
    plot1 <- VariableFeaturePlot(object = hv_pbmc())
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
  })
  # Show highly variable genes 
  output$hvg_table = DT::renderDataTable({
    if(is.null(hv_pbmc())) return()
    else if(input$hvg_ok == FALSE) return()
    else if(input$hvg_disp == FALSE) return() 
    #h_info<-as.matrix(HVFInfo(object = hv_pbmc()[VariableFeatures(object = hv_pbmc()),]))
    
    h_info<- as.data.frame(HVFInfo(object = hv_pbmc()[VariableFeatures(object = hv_pbmc()),]))
    h_info<- h_info[order(-h_info$variance.standardized),]
    
    h_info<- h_info %>% 
      dplyr:: rename(
        Standardized.variance = variance.standardized
      )
    h_info
    #temp_data<-GetAssayData(object = hv_pbmc())[VariableFeatures(object = hv_pbmc()),1:2]
    #data.matrix(temp_data)
  })
  # Download High variable genes data set.........
  output$hvg_download<- downloadHandler(
    if(is.null(hv_pbmc())) return()
    else if(input$hvg_ok == FALSE) return()
    else if(input$high_hvg == FALSE) return(),
    filename = function() {
      paste("highly_variable_genes", ".csv", sep = "")
    },
    content = function(file) {
      hv_data<-GetAssayData(object = hv_pbmc())[VariableFeatures(object = hv_pbmc()),]
      write.csv(as.matrix(hv_data), file, row.names = TRUE, na = "")
    }
  )
  
  #... scaling the Normalized data
  scdata <- reactive({
    if(is.null(hv_pbmc())) return()
    else if(input$hvg_ok == FALSE) return()
    #isolate({
      #print(hv_pbmc())
      all_genes <- rownames(x = hv_pbmc())
      ScaleData(object = hv_pbmc(), features = all_genes)
    #})
    
  })
  #----------------------------------------Seurat PCA--------------------------
  pc_pbmc<-reactive({
    #if(is.null(hv_pbmc())) return()
    if(input$com_pca==FALSE) return()
    isolate({
      
      if(input$pc_genes == "pc_hvg")
      {
        
        #if(input$pc_impute == 1)
        #{
        #MAGIC_pbmc_obj <- magic(hv_pbmc(), genes=hv_pbmc()@var.genes)
        # MAGIC_pbmc_obj <- RunPCA(MAGIC_pbmc_obj,pcs.compute = input$pc_cells)
        #MAGIC_pbmc_obj
        #}
        #else{
        
        
        #scdata<-ScaleData(object = hv_pbmc())
        RunPCA(object = scdata(), features = VariableFeatures(object = scdata()),
               npcs = input$pc_cells, verbose = FALSE)
        #}
      } # pc_hvg
      
      else if(input$pc_genes == "pc_all")
      {
        
        all.genes <- rownames(x = scdata())
        #pbmc_all <- ScaleData(object = hv_pbmc(), features = all.genes)
        
        RunPCA(object = scdata(), features = all.genes,
               npcs = input$pc_cells, verbose = FALSE)
      }
    }) # end isolate
  })
  # after PCA notification 
  output$pc_text<-renderUI({
    if(is.null(pc_pbmc())) return()
    #else if(input$clus_ok == FALSE) return()
    str1<-paste("PCA Completed!")
    HTML(str1)
  })
  #PCAplot
  output$sc_pca_plot<-renderPlot({
    if(is.null(pc_pbmc())) return() 
    else if(input$ck_pca==FALSE) return() #hide("sc_pca_plot")
    #PCAPlot(object = pc_pbmc(), dim.1 = 1, dim.2 = 2)
    DimPlot(object = pc_pbmc(), reduction = "pca")
  })
  
  #PC t-SNE plot
  output$pc_tsne_plot<-renderPlot({
    if(is.null(pc_pbmc())) return() 
    else if(input$ck_tsne_plot==FALSE) return() #hide("sc_pca_plot")
    t_tsne<- RunTSNE(object = pc_pbmc(), reduction = "pca", 
                     dims.use = 1:input$pc_cells, do.fast = TRUE)
    DimPlot(object = t_tsne, reduction = "tsne")
    
  })
  
  # PCA Elbow plot
  output$pc_elbow_plot<-renderPlot({
    if(input$pc_elbow == FALSE) return()
    ElbowPlot(object = pc_pbmc(), ndims = input$pc_cells)
  })
  
  # print PCA
  
  output$pca_list<-renderPrint({
    if(is.null(pc_pbmc())) return()
    else if(input$pc_print==FALSE) return()
    #tp<-PrintPCA(object = pc_pbmc(), pcs.print = 1:input$pc_cells, genes.print = NULL, use.full = FALSE)
    tp<-print(x = pc_pbmc()[["pca"]], dims = 1:input$pc_cells, nfeatures = NULL)
    tp
  })
  
  # PCA Jack Straw Plot 
  output$pc_jac_show<-renderPlot({
    if(is.null(pc_pbmc())) return()
    else if(input$pc_jac_plot == FALSE) return()
    # pdata <- ProjectPCA(object = pc_pbmc(), do.print = T)
    
    jplot <- JackStraw(object = pc_pbmc(), num.replicate = 100, dims = input$pc_cells)
    jplot <- ScoreJackStraw(object = jplot, dims = 1:input$pc_cells)
    JackStrawPlot(object = jplot, dims = 1:input$pc_cells)
  })
  
  #output$tsnep<- renderPlot({
  # if(is.null(scdata())) return()
  #else if(input$ck_tsne == FALSE) return()
  
  #if(input$tsn_genes == "tsn_hvg"){
  # var_gene<- VariableFeatures(object = hv_pbmc())
  # m_sub<- GetAssayData(object = scdata(), slot = "scale.data")[var_gene,] # used scaled data
  #}
  #else 
  # {
  #all_genes <- rownames(x = hv_pbmc())
  #  m_sub<- GetAssayData(object = scdata(), slot = "scale.data") #[all_genes,] # used scaled data
  
  #}
  
  #tm_sub<- t(m_sub)
  #colour1 = rainbow(length(rownames(tm_sub)))
  #pt<-Rtsne(tm_sub,dims=2, perplexity=1,check_duplicates = F)
  #tsne_plot<- data.frame(x = pt$Y[,1], y = pt$Y[,2], col = colour1)
  #ggplot(tsne_plot) + geom_point(aes(x=x, y=y), col=colour1, size = 3)
  # })
  
  # PCA Heatmap
  output$pc_hplot<-renderPlot({
    if(input$pc_hmap_ok == FALSE) return()
    #isolate({
      #DimHeatmap(object = pbmc, dims = 1:15, cells = 500, balanced = TRUE)
      DimHeatmap(object = pc_pbmc(), dims = input$pc_low:input$pc_upper, nfeatures = input$pc_hmap_genes,
                 balanced = TRUE, fast = TRUE)
    #})
  })
  
  # ----------------------------------------------------- Clustering ------------------------------
  clus_pbmc<-reactive({
    if(is.null(pc_pbmc())) return()
    else if (input$clus_ok == FALSE) return()
    isolate({
      if(input$clus_algo == "Louvain Algorithm(LA)")
        algo_type = 1
      else if(input$clus_algo == "Multilevel LA")
        algo_type = 2
      else if(input$clus_algo == "SLM algorithm")
        algo_type = 3
      else if(input$clus_algo == "Leiden algorithm")
        algo_type = 4
      print(pc_pbmc())
      ct_data<- FindNeighbors(object = pc_pbmc(), dims = input$clus_low:input$clus_upper)
      
      ct_data<- FindClusters(object = ct_data, graph.name = NULL,
                             modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                             node.sizes = NULL, resolution = input$clus_restn, algorithm = algo_type, n.start = 10,
                             n.iter = 10, random.seed = 0, temp.file.location = NULL,
                             edge.file.name = NULL, verbose = TRUE)
      
      #ct_data<- RunUMAP(object = ct_data, dims = input$clus_low:input$clus_upper)
      
      ct_data<- RunTSNE(object = ct_data, reduction = "pca", 
                        dims.use = input$clus_low:input$clus_upper, do.fast = TRUE)
      ct_data
    })
  })
  
  # after clustering notification 
  output$clus_text<-renderUI({
    if(is.null(clus_pbmc())) return()
    #else if(input$clus_ok == FALSE) return()
    str1<-paste("Clustering Completed")
    HTML(str1)
  })
  # t-SNE plot
  output$clus_tsne_plot<-renderPlot({
    if(is.null(clus_pbmc())) return()
    else if (input$clus_tsne == FALSE) return()
    isolate({
      
      # Calculate number of cells per cluster from object@ident
      cell.num <- table(Idents(clus_pbmc()))
      
      # Add cell number per cluster to cluster labels
      ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))
      
      # Order legend labels in plot in the same order as 'ClusterLabels'
      ClusterBreaks = names(cell.num)
      
      # Plot tSNE with new legend labels for clusters
      #TSNEPlot(object = clus_pbmc(), do.return = T) +
      #   scale_colour_discrete(breaks = ClusterBreaks, 
      #                        labels = ClusterLabels) +
      # labs(x = "tSNE_1",
      #      y = "tSNE_2")
      
      DimPlot(object = clus_pbmc(), reduction = "tsne") + 
        scale_colour_discrete(breaks = ClusterBreaks, 
                              labels = ClusterLabels) +
        labs(x = "tSNE_1",
             y = "tSNE_2")
    })
  })
  
  # UMAP plot
  
  output$clus_umap_plot<-renderPlot({
    if(is.null(clus_pbmc())) return()
    else if (input$clus_umap == FALSE) return()
    isolate({
      #DimPlot(object = clus_pbmc(), reduction = "umap")
      UMAPPlot(clus_pbmc())
    })
  })
  
  # ..........................Cluster Bar Plot  
  output$clus_barplot_show<-renderPlot({
    if(is.null(clus_pbmc())) return() 
    else if(input$clus_barplot==FALSE) return() 
    
    #How many cells are in each cluster
    bar_data<-table(Idents(clus_pbmc()))
    #dim(x  = bar_data)
    bar_data<-data.frame(bar_data)
    
    bar_data<- bar_data %>% 
      dplyr:: rename(
        Clusters = Var1
      )
    # bar plor show
    clus_bar_plot<- ggplot(bar_data, aes( x = Clusters, y = Freq, fill = Clusters )) + geom_bar(stat="identity") +
      scale_fill_hue(c = 100) +
      labs(x = "Clusters", 
           y = "Number of Cells") +
      ggtitle("Number of cells in the each cluster") +
      geom_text(aes(label=bar_data$Freq), position=position_dodge(width=1), vjust=-0.25, col= "black")
    
    barp1<-clus_bar_plot+ theme(plot.title = element_text(hjust = 0.5), 
                                # Hide panel borders and remove grid lines
                                panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                # Change axis line
                                axis.line = element_line(colour = "black"))
    barp1 + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(color="#993333", size=14, face="bold"),
      axis.title.y = element_text(color="blue", size=14, face="bold")
    )
    
    
  })
  # -------------------------------------------------------- Differential Expression Analysis ---------------
  output$DE_clus_UI <- renderUI({
    
    tabsetPanel(id = "de_tabset",
                 
                 tabPanel("Find all markers",
                          img(src = "line_font1.png"),
                          sidebarPanel(id = "allclus_slid", width = 6, height = 1000, 
                                       tags$style("#allclus_slid{background-color:#F6F9CD;}"),
                                       br(),
                                       
                                       
                                       
                                       column(width = 12,selectizeInput("all_de_test", "Select Test Method", 
                                                                        c(Wilcox = "wilcox",
                                                                          Bimod='bimod',
                                                                          t_test='t',
                                                                          Logistic_Regression= 'LR',
                                                                          MAST = 'MAST'),
                                                                        selected = 'wilcox')),        
                                       
                                       tags$head(
                                         tags$style(HTML('#all_de_ok{background-color:red;color:white;}'))
                                       ),
                                       actionButton('all_de_ok', 'Find All markers'),
                                       
                                       #tags$hr(),
                                       uiOutput("DE_download_UI"),
                                       
                                       #tags$hr(),
                                       checkboxInput('all_de_show',p('Show All Markers')),
                                       
                                       #tags$hr(),
                                       
                                       
                                       #.....................
                                       checkboxInput('de_filter',p('Filtering Markers')),
                                       conditionalPanel('input.de_filter == 1',
                                                        column(width = 6, numericInput("all_clst_thres","Avg.logFC threshold",0.25,min = 0,max = Inf)),
                                                        column(width = 6, numericInput("all_clst_pct","Min % (min.pct)",0.1,min = 0,max = Inf)),
                                                        column(width = 12, numericInput("all_clst_pvalue","Adjust p-value",0.05,min = 0,max = Inf)),
                                                        #actionButton('de_ok_filter','Search Markers'),
                                                        #tags$hr(),
                                                        radioButtons("de_selection", "Markers Selection",
                                                                     choices = c('Positive only' = "de_pos",'Negative only' = "de_neg"),selected = "de_pos"),
                                                        checkboxInput('de_filter_show',p('Show filtered list')),
                                                        
                                                        checkboxInput('all_de_show_top',p('Show Top genes')),
                                                        conditionalPanel('input.all_de_show_top == 1',
                                                                         column(width = 6, numericInput("all_de_top_gene", "Top DE genes:", 5, min = 1, max = Inf)),
                                                                         tags$head(
                                                                           tags$style(HTML('#top_genes_download{background-color:#800080;color:white;}'))
                                                                         ),
                                                                         downloadButton("top_genes_download", "Save All Top genes as a text file"),
                                                                         #tags$hr(),
                                                                         downloadButton("top_pos_genes_download", "Save +ve Top genes as a text file"),
                                                                         downloadButton("top_neg_genes_download", "Save -ve Top genes as a text file"),
                                                                         #tags$hr(),
                                                                         checkboxInput('de_filter_hmap',p('Show Heatmap(top genes)'))
                                                        )
                                                        
                                                        
                                                        
                                       ), 
                                       
                                       # tags$hr(),
                                       
                                       radioButtons("de_ck_barplot", "Markers Barplot",
                                                    choices = c('All DE genes' = "de_all_bar",'Filtered DE genes' = "de_filt_bar", 
                                                                'Hide Bar Plot' = 'hide_de_bar'),
                                                    selected = "de_all_bar")
                                       
                                       
                                       
                                       #,
                                       
                                       
                                       #conditionalPanel('input.all_de_show == 1',
                                       #downloadButton("all_de_download", "Download"))
                          ),
                          column(width =6,br(), br(),h6(htmlOutput("de_text"))),
                          #column( width = 6 , downloadButton("all_de_download", "Download")),
                          column(width = 6, br(), br(),DT::dataTableOutput("all_de_list")), 
                          column(width = 6, br(), br(),DT::dataTableOutput("filt_de_list")),
                          column(width =12, br(), br(), plotOutput("all_de_barplot")),
                          column(width = 6, br(), br(),DT::dataTableOutput("all_de_show_top_list")),
                          column(width =12, br(), br(), plotOutput("de_filter_hmap_plot"))
                 ), # tabPanel of Find all markers 
                 
                 
                 tabPanel("Find markers by cluster",
                          img(src = "line_font1.png"),
                          sidebarPanel(id = "fclus_slid",width = 6,
                                       tags$style("#fclus_slid{background-color:#FAE4E0;}"),
                                       "Please select parameters:",
                                       #tags$hr(),
                                       br(),
                                       tags$head(
                                         tags$style(HTML('#sclst1{background-color:#0000FF;color:white;font-size: 15px;}'))
                                       ),
                                       
                                       column( width = 5,selectInput('sclst1','Select cluster',choices = levels(x = clus_pbmc()),selected = '0')),
                                       #column( width = 5,numericInput("sclst1", "Select cluster:", 0, min = 0, max = Inf)),
                                       column(width = 6, numericInput("sclst_thres","logfc threshold",0.25,min = 0,max = Inf)),
                                       column(width = 12, numericInput("sclst_pct","Min % (min.pct)",0.1,min = 0,max = Inf)),
                                       
                                       column(width = 12,selectizeInput("sde_test", "Select Test Method", 
                                                                        c(Wilcox = "wilcox",
                                                                          Bimod='bimod',
                                                                          t_test='t',
                                                                          Logistic_Regression= 'LR',
                                                                          MAST = 'MAST'),
                                                                        selected = 'wilcox')),
                                       tags$head(
                                         tags$style(HTML('#sde_ok{background-color:black;color:white;}'))
                                       ),
                                       actionButton('sde_ok', 'Find Markers'),
                                       tags$hr(),
                                       
                                       
                                       uiOutput("sDE_download_UI"),
                                       
                                       
                                       tags$hr(),
                                       checkboxInput('sde_show',p('Show All Markers(both positive and negative)')),
                                       
                                       
                                       checkboxInput('sde_filter',p('Filtering Markes')),
                                       conditionalPanel('input.sde_filter == 1',
                                                        column(width = 12, numericInput("sde_clst_pvalue","Adjust p-value",0.05,min = 0,max = Inf)),
                                                        tags$hr(),
                                                        radioButtons("sde_selection", "Markers Selection",
                                                                     choices = c('Positive only' = "sde_pos",'Negative only' = "sde_neg"),selected = "sde_pos"),
                                                        
                                                        checkboxInput('sde_filter_hmap',p('Show Heatmap')),
                                                        conditionalPanel('input.sde_filter_hmap == 1',
                                                                         column(width = 6, numericInput("sde_top","Top Genes:", 10, min = 2,max = Inf)),
                                                                         checkboxInput('sde_filter_hmap_cluster',p('Zoom the cluster'))
                                                        ),
                                                        
                                                        checkboxInput('sde_filter_show',p('Show Filtered Markers'))
                                                        
                                       )
                                       
                          ),
                          
                          column(width = 6, br(), br(),DT::dataTableOutput("sde_list")),
                          column(width = 6, br(), br(),DT::dataTableOutput("sde_filtered_list")),
                          column(width =12, br(), br(), plotOutput("sde_filter_hmap_plot")),
                          column(width =12, br(), br(), plotOutput("sde_filter_hmap_cluster_plot"))
                          
                          
                 ), # End tabPanel for markers by cluster
                 
                 # ---------------------------------------------------------   
                 tabPanel("Find markers by clusters vs other clusters",
                          img(src = "line_font1.png"),
                          sidebarPanel(id = "fclus_slid2",width = 6, 
                                       tags$style("#fclus_slid2{background-color:#DAF3DA;}"),
                                       "Please select parameters:",
                                       tags$hr(),
                                       #br(),
                                       tags$head(
                                         tags$style(HTML('#clst1{background-color:#F0FFF0;color:black;}'))
                                       ),
                                       
                                       #column( width = 5,selectInput('clst1','Select cluster',choices = levels(x = clus_pbmc()),selected = '0')),
                                       
                                       column(width = 12, textAreaInput("clst1","Select Cluster/s", placeholder =  "The cluster number must be separated by a commma")),
                                       
                                       tags$head(
                                         tags$style(HTML('#clst2{background-color:#800080;color:white;}'))
                                       ),
                                       #column( width = 6,numericInput("clst2", "Select complementary cluster:", 0, min = 0, max = Inf)),
                                       #column( width = 6,selectInput('clst2','Select complementary cluster:',choices = levels(x = clus_pbmc()),selected = '1')),
                                       br(),
                                       tags$hr(),
                                       br(),
                                       column(width = 12, textAreaInput("clst2","Select complementary cluster/s", placeholder = "The cluster number must be separated by a commma")),
                                       
                                       column(width = 6, numericInput("clst_thres","logfc threshold",0.25,min = 0,max = Inf)),
                                       column(width = 6, numericInput("clst_pct","Min % (min.pct)",0.1,min = 0,max = Inf)),
                                       
                                       column(width = 12,selectizeInput("de_test", "Select Test Method", 
                                                                        c(Wilcox = "wilcox",
                                                                          Bimod='bimod',
                                                                          t_test='t',
                                                                          Logistic_Regression= 'LR',
                                                                          MAST = 'MAST'),
                                                                        selected = 'wilcox')),        
                                       
                                       tags$head(
                                         tags$style(HTML('#de_ok{background-color:#00008B;color:white;}'))
                                       ),
                                       actionButton('de_ok', 'Find markers'),
                                       #conditionalPanel('input.de_ok == 1',
                                       tags$hr(),
                                       
                                       tags$head(
                                         tags$style(HTML('#fde_download{background-color:blue;color:white;}'))
                                       ),
                                       
                                       checkboxInput('de_show',p('Show All Markers(Both Positive and Negative)')),
                                       checkboxInput('fde_filter',p('Filtering Markers')),
                                       conditionalPanel('input.fde_filter == 1',
                                                        #tags$hr(),
                                                        column(width = 12, numericInput("fde_clst_pvalue","Adjust p-value",0.05,min = 0,max = Inf)),
                                                        radioButtons("fde_selection", "Markers Selection",
                                                                     choices = c('Positive only' = "fde_pos",'Negative only' = "fde_neg"
                                                                     ),selected = "fde_pos"),
                                                        
                                                        downloadButton("fde_download", "Save Filtered Markers as csv"),
                                                        
                                                        checkboxInput('fde_filter_hmap',p('Show Heatmap')),
                                                        conditionalPanel('input.fde_filter_hmap == 1',
                                                                         column(width = 6, numericInput("fde_top","Top Genes:", 10, min = 2,max = Inf))              
                                                        ),
                                                        
                                                        checkboxInput('fde_filter_show',p('Show Filtered Markers'))
                                       )
                                       
                                       #)
                          ),
                          
                          column(width = 6, br(), br(),DT::dataTableOutput("de_list")),
                          column(width = 6, br(), br(),DT::dataTableOutput("fde_filtered_list")),
                          column(width =12, br(), br(), plotOutput("fde_filter_hmap_plot"))
                 )
                 
                 
    ) # De tabset
    
  })
  
  # .................................................................... find all markers  ............................
  all_markers<-reactive({
    if(is.null(clus_pbmc())) return()
    else if (input$all_de_ok == FALSE) return()
    isolate({
      
      FindAllMarkers(clus_pbmc(), assay = NULL, features = NULL,
                     logfc.threshold = 0, test.use = input$all_de_test, slot = "data",
                     min.pct = 0, min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                     only.pos = FALSE, max.cells.per.ident = Inf, random.seed = 1,
                     latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                     pseudocount.use = 1, return.thresh = 0.01)
      
    })
  })
  
  # after de analysis notification 
  output$de_text<-renderUI({
    if(is.null(all_markers())) return()
    str1<-paste("DE Analysis Completed!!")
    HTML(str1)
  })
  
  # Show All DE genes # Show List
  output$all_de_list = DT::renderDataTable({
    if(is.null(all_markers())) return()
    else if (input$all_de_show == FALSE) return()
    #all_de_genes<-data.matrix(all_markers())
    #all_de_genes
    dplyr::select(all_markers(), -gene)
    
  })
  
  
  # Download all DE gene list .........
  output$all_de_download<- downloadHandler(
    if(is.null(all_markers())) return(),
    filename = function() {
      paste("all_de_genes", ".csv", sep = "")
    },
    content = function(file) {
      tgs<-dplyr::select(all_markers(), -gene)
      all_de_data<-data.matrix(tgs)
      write.csv(all_de_data, file, row.names = TRUE)
    }
  )
  
  output$DE_download_UI <- renderUI({
    
    if(!is.null(all_markers())){
      tags$head(
        tags$style(HTML('#all_de_download{background-color:#800080;color:white;}'))
      )
      downloadButton("all_de_download", "Save All markers as a csv")
    }
  })
  
  
  #.............. Markers filtering
  fmarkers<-reactive({
    if(is.null(all_markers())) return()
    else if(input$de_selection == "de_pos")
      filter(all_markers(),avg_logFC>=input$all_clst_thres, pct.1>=input$all_clst_pct , p_val_adj<= input$all_clst_pvalue)
    
    else if(input$de_selection == "de_neg")
      filter(all_markers(),avg_logFC<=(input$all_clst_thres)*(-1), pct.1>=input$all_clst_pct, p_val_adj<= input$all_clst_pvalue)
    
    #else {
    #filter(all_markers(),abs(avg_logFC)>=input$all_clst_thres, pct.1>=input$all_clst_pct, p_val_adj<= input$all_clst_pvalue) 
    #}
  })
  
  # .....................Show Fitered Markers
  output$filt_de_list = DT::renderDataTable({
    if(is.null(fmarkers())) return()
    else if (input$de_filter_show == FALSE) return()
    
    #else if (input$de_ok_filter == FALSE) return()
    #dplyr::select(fmarkers(), -gene)
    #filt_de_genes<-data.matrix(fmarkers())
    #filt_de_genes
    #isolate({
    fmarkers()
    
    #})
    
  })
  
  # .........Show Top DE genes
  output$all_de_show_top_list = DT::renderDataTable({
    if(is.null(fmarkers())) return()
    else if (input$all_de_show_top == FALSE) return()
    
    if(input$de_selection == "de_pos") {
      pos_top_genes_hmap<- fmarkers() %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene , wt = avg_logFC)
      sort_pos_top_genes_hmap<- pos_top_genes_hmap[order(-pos_top_genes_hmap$avg_logFC),] # descending order 
      sort_pos_top_genes_hmap <- sort_pos_top_genes_hmap[order(sort_pos_top_genes_hmap$cluster),] # descending order by cluster
      sort_pos_top_genes_hmap
    }
    
    else if(input$de_selection == "de_neg") {
      neg_top<- fmarkers()
      neg_top$avg_logFC <- neg_top$avg_logFC * (-1) # make positive to find top genes
      neg_tgenes<- neg_top %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene, wt = avg_logFC)
      #dim(neg_tgenes)
      #print.data.frame(neg_tgenes) # before sort
      sort_neg_tgenes<-neg_tgenes[order(-neg_tgenes$avg_logFC),]
      sort_neg_tgenes <- sort_neg_tgenes[order(sort_neg_tgenes$cluster),] # sort by cluster 
      #print.data.frame(sort_neg_tgenes) # after sort
      sort_neg_tgenes$avg_logFC <- sort_neg_tgenes$avg_logFC * (-1)
      sort_neg_tgenes
    }
    
  })
  
  # ...... All Top DE genes download 
  output$top_genes_download<- downloadHandler(
    if(is.null(all_markers())) return(),
    filename = function() {
      paste("top_de_genes", ".txt", sep = "")
    },
    
    content = function(file) {
      
      tgenesd<- fmarkers() %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene, wt = avg_logFC)
      
      # gene name
      tgene_name<-dplyr::select(tgenesd,gene)
      #dim(x = tgene_name)
      genes_list<- tgene_name[,'gene']
      genes_list<- data.frame(genes_list)
      
      ugenes<-unique(genes_list[,1:1])
      
      #fileConn<-file("Top_genes.txt")
      writeLines(c(ugenes), file)
      #close(file)
      
    }
  )
  
  # ............ Positive Top DE genes download
  output$top_pos_genes_download<- downloadHandler(
    if(is.null(all_markers())) return(),
    filename = function() {
      paste("positive_top_de_genes", ".txt", sep = "")
    },
    
    content = function(file) {
      
      tgenesd<- all_markers() %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene, wt = avg_logFC)
      
      # gene name
      tgene_name<-dplyr::select(tgenesd,gene)
      #dim(x = tgene_name)
      genes_list<- tgene_name[,'gene']
      genes_list<- data.frame(genes_list)
      
      ugenes<-unique(genes_list[,1:1])
      
      #fileConn<-file("Top_genes.txt")
      writeLines(c(ugenes), file)
      #close(file)
      
    }
  )
  
  # ............ Negative Top DE genes download
  output$top_neg_genes_download<- downloadHandler(
    if(is.null(all_markers())) return(),
    filename = function() {
      paste("negative_top_de_genes", ".txt", sep = "")
    },
    
    content = function(file) {
      
      tgenesd<- all_markers() %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene, wt = avg_logFC)
      
      # gene name
      tgene_name<-dplyr::select(tgenesd,gene)
      #dim(x = tgene_name)
      genes_list<- tgene_name[,'gene']
      genes_list<- data.frame(genes_list)
      
      ugenes<-unique(genes_list[,1:1])
      
      #fileConn<-file("Top_genes.txt")
      writeLines(c(ugenes), file)
      #close(file)
      
    }
  )
  
  # DE filtered Heatmap 
  
  output$de_filter_hmap_plot <- renderPlot({
    if(is.null(fmarkers())) return()
    else if (input$de_filter_hmap == FALSE) return()
    
    if(input$de_selection == "de_pos") {
      pos_top_genes_hmap<- fmarkers() %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene , wt = avg_logFC)
      sort_pos_top_genes_hmap<- pos_top_genes_hmap[order(-pos_top_genes_hmap$avg_logFC),] # descending order 
      sort_pos_top_genes_hmap <- sort_pos_top_genes_hmap[order(sort_pos_top_genes_hmap$cluster),] # descending order by cluster
      DoHeatmap(clus_pbmc(), features = sort_pos_top_genes_hmap$gene) + NoLegend() # heatmap
    }
    
    else if(input$de_selection == "de_neg") {
      neg_top<- fmarkers()
      neg_top$avg_logFC <- neg_top$avg_logFC * (-1)
      neg_tgenes<- neg_top %>% group_by(cluster) %>% top_n(n = input$all_de_top_gene, wt = avg_logFC)
      #dim(neg_tgenes)
      #print.data.frame(neg_tgenes) # before sort
      sort_neg_tgenes<-neg_tgenes[order(-neg_tgenes$avg_logFC),]
      sort_neg_tgenes <- sort_neg_tgenes[order(sort_neg_tgenes$cluster),] # sort by cluster 
      #print.data.frame(sort_neg_tgenes) # after sort
      sort_neg_tgenes$avg_logFC <- sort_neg_tgenes$avg_logFC * (-1) # add negative sign 
      DoHeatmap(clus_pbmc(), features = sort_neg_tgenes$gene, angle = 45) # heatmap
    }
    
  })
  
  #..................... DE bar plot 
  
  output$all_de_barplot <- renderPlot({
    
    if(is.null(all_markers())) return()
    else if (input$all_de_ok == FALSE) return()
    
    #How many DE genes in  the each cluster
    else if (input$de_ck_barplot == "de_all_bar"){
      degenes<- data.frame(filter(dplyr::select(all_markers() ,gene, cluster)))
    }
    
    else if (input$de_ck_barplot == "de_filt_bar"){
      degenes<- data.frame(filter(dplyr::select(fmarkers() ,gene, cluster)))
    }
    
    else if (input$de_ck_barplot == "hide_de_bar") { return ()}
    
    
    #dim(x  = degenes)
    
    #degenes[1:5,]
    
    fr<-table(degenes$cluster)
    #fr
    frd<-data.frame(fr)
    
    # frd
    
    frd<- frd %>% 
      dplyr:: rename(
        Clusters = Var1
      )
    #frd
    
    # DE genes cluster bar plot
    
    de_clus_bar_plot<- ggplot(frd, aes( x = Clusters, y = Freq, fill = Clusters )) + geom_bar(stat="identity") +
      scale_fill_hue(c = 200) +
      labs(x = "Clusters", 
           y = "Number of Genes") +
      ggtitle("Number of Genes in the each cluster") +
      geom_text(aes(label=frd$Freq), position=position_dodge(width=1), vjust=-0.25, col= "black")
    
    de_barp1<-de_clus_bar_plot+ theme(plot.title = element_text(hjust = 0.5), 
                                      # Hide panel borders and remove grid lines
                                      panel.border = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      # Change axis line
                                      axis.line = element_line(colour = "black"))
    de_barp1 + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(color="#993333", size=14, face="bold"),
      axis.title.y = element_text(color="blue", size=14, face="bold")
    )
    
    
  })
  # ................................................................ find markers by a specific cluster .........................
  
  scluster_markers<-reactive({
    if(is.null(clus_pbmc())) return()
    else if (input$sde_ok == FALSE) return()
    isolate({
      resd<-FindMarkers(clus_pbmc(), ident.1 = input$sclst1,
                        group.by = NULL, subset.ident = NULL, assay = NULL,
                        slot = "data", reduction = NULL, features = NULL,
                        logfc.threshold = input$sclst_thres, test.use = input$sde_test, min.pct = input$sclst_pct,
                        min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE,
                        max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
                        min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1)
      
      resd[,'Gene']<-rownames(resd)
      resd
    })
  })
  
  
  
  # Download DE gene list .........
  output$all_sde_download<- downloadHandler(
    if(is.null(scluster_markers())) return()
    else if(input$sde_ok == FALSE) return(),
    filename = function() {
      paste("sde_genes", ".csv", sep = "")
    },
    content = function(file) {
      tgs<-dplyr::select(scluster_markers(), -Gene)
      all_de_data<-data.frame(tgs)
      write.csv(all_de_data, file, row.names = TRUE)
    }
  )
  
  output$sDE_download_UI <- renderUI({
    
    if(!is.null(scluster_markers())){
      tags$head(
        tags$style(HTML('#all_sde_download{background-color:#800080;color:white;}'))
      )
      downloadButton("all_sde_download", "Save All markers as a csv")
    }
  })
  
  
  # Show All Markers
  output$sde_list = DT::renderDataTable({
    
    if(is.null(scluster_markers())) return()
    
    else if (input$sde_show == FALSE) return()
    
    #dplyr::select(scluster_markers(),Gene, p_val,avg_logFC,pct.1,pct.2,p_val_adj)
    scluster_markers()
  })
  
  
  
  #.............. Markers filtering
  sde_fmarkers<-reactive({
    if(is.null(scluster_markers())) return()
    
    else if(input$sde_selection == "sde_pos"){
      sfd<- filter(scluster_markers(),avg_logFC>=0, p_val_adj<= input$sde_clst_pvalue)
      sfd <- sfd[order(-sfd$avg_logFC),]
    }
    
    else if(input$sde_selection == "sde_neg"){
      sfd<- filter(scluster_markers(),avg_logFC<0, p_val_adj<= input$sde_clst_pvalue)
      sfd<- sfd[order(sfd$avg_logFC),]
    }
    
    dplyr::select(sfd,Gene, p_val,avg_logFC,pct.1,pct.2,p_val_adj)
  })
  
  # .....................Show Fitered Markers
  output$sde_filtered_list = DT::renderDataTable({
    if(is.null(sde_fmarkers())) return()
    else if (input$sde_filter_show == FALSE) return()
    
    sde_fmarkers()
    
  })
  
  # DE filtered Heatmap for a specific cluster 
  
  output$sde_filter_hmap_plot <- renderPlot({
    if(is.null(sde_fmarkers())) return()
    else if (input$sde_filter_hmap == FALSE) return()
    
    sde_top<- sde_fmarkers()[1:input$sde_top,]
    
    DoHeatmap(clus_pbmc(), features = sde_top$Gene, group.bar = T) # heatmap
    
    
  })
  
  # DE filtered Heatmap for a specific cluster : only for the particular cluster 
  
  output$sde_filter_hmap_cluster_plot <- renderPlot({
    if(is.null(sde_fmarkers())) return()
    else if (input$sde_filter_hmap_cluster == FALSE) return()
    
    sde_top<- sde_fmarkers()[1:input$sde_top,]
    hmap_cells<- WhichCells(clus_pbmc(), idents = input$sclst1) # identify cells for a specific cluster
    DoHeatmap(clus_pbmc()[,hmap_cells], features = sde_top$Gene, group.bar = T) # heatmap
    
    
  })
  # ................................................. find markers: cluster vs other cluster..........
  
  cluster_markers<-reactive({
    if(is.null(clus_pbmc())) return()
    else if (input$de_ok == FALSE) return()
    
    isolate({
      cls_Num1<-unlist(strsplit(input$clst1, ","))
      
      cls_Numi1<-as.integer(cls_Num1)
      
      cls_Numi1<-cls_Numi1+1
      #print(fd_clus)
      clus_NUmf1<-levels( x= clus_pbmc())[cls_Numi1] 
      #clus_NUmf
      
      #print(clus_NUmf)
      #.............
      cls_Num2<-unlist(strsplit(input$clst2, ","))
      
      cls_Numi2<-as.integer(cls_Num2)
      
      cls_Numi2<-cls_Numi2+1
      
      clus_NUmf2<-levels( x= clus_pbmc())[cls_Numi2] 
      
      resd<- FindMarkers(clus_pbmc(), ident.1 =  clus_NUmf1 , ident.2 = clus_NUmf2,
                         group.by = NULL, subset.ident = NULL, assay = NULL,
                         slot = "data", reduction = NULL, features = NULL,
                         logfc.threshold = input$clst_thres, test.use = input$de_test, min.pct = input$clst_pct,
                         min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE,
                         max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
                         min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1)
      
      resd[,'Gene']<-rownames(resd)
      resd
    })
  })
  
  # Show  All DE genes: clusters vs clusters
  output$de_list = DT::renderDataTable({
    if(is.null(cluster_markers())) return()
    else if (input$de_show == FALSE) return()
    #det<-data.matrix(fdegenes())
    #det
    cluster_markers()
    
  })
  #  Fitering Markers clusters vs clusters 
  fdegenes<- reactive({
    if(is.null(cluster_markers())) return()
    else if(input$fde_selection == "fde_pos")
    {
      sfd<- filter(cluster_markers(),avg_logFC>=0, p_val_adj<= input$fde_clst_pvalue)
      sfd <- sfd[order(-sfd$avg_logFC),]
    }
    
    else if(input$fde_selection == "fde_neg"){
      sfd<- filter(cluster_markers(),avg_logFC<0, p_val_adj<= input$fde_clst_pvalue)
      sfd<- sfd[order(sfd$avg_logFC),]
    }
    
    dplyr::select(sfd,Gene, p_val,avg_logFC,pct.1,pct.2,p_val_adj)
    
  })
  
  # Show  Filtered DE genes: clusters vs clusters
  output$fde_filtered_list = DT::renderDataTable({
    if(is.null(fdegenes())) return()
    else if (input$fde_filter_show == FALSE) return()
    #det<-data.matrix(fdegenes())
    #det
    fdegenes()
    
  })
  
  # Download filtered DE gene list .........
  output$fde_download<- downloadHandler(
    if(is.null(fdegenes())) return(),
    
    filename = function() {paste("Filtered_Markers_Clusters_vs_Clusters", ".csv", sep = "")},
    content = function(file) {
      #de_data<-data.matrix(fdegenes())
      write.csv(fdegenes(), file, row.names = TRUE)}
  )
  
  # DE filtered Heatmap for  clusters vs Clusters 
  
  output$fde_filter_hmap_plot <- renderPlot({
    if(is.null(fdegenes())) return()
    else if (input$fde_filter_hmap == FALSE) return()
    
    cls_Num1<-unlist(strsplit(input$clst1, ","))
    
    cls_Numi1<-as.integer(cls_Num1) # cluster number
    
    fde_top<- fdegenes()[1:input$fde_top,]
    print(cls_Numi1)
    fd_hmap_cells<- WhichCells(clus_pbmc(), idents = cls_Numi1)
    DoHeatmap(clus_pbmc()[,fd_hmap_cells], features = fde_top$Gene, group.bar = T) # heatmap
    
  })
  
  # -------------------------------------------------------------------- Violin Plot --------------------------------------------------
  
  vln_gene_list<- reactive({
    
    gene_list<-  unlist(strsplit(input$vln_genes, "\n"))
    if(is.null(gene_list)) return()
    rev(gene_list)
    
  })
  
  output$vln_plot<- renderPlot({
    
    if(is.null(clus_pbmc())) return() 
    
    #vlngne<- rownames(GetAssayData(object = clus_pbmc(), slot = "scale.data"))
    gene_list<- vln_gene_list()
    
    #print(gene_list)
    
    #VlnPlot(clus_pbmc(), slot = "scale.data", features = gene_list) 
    VlnPlot_2 <- function(object, features.plot,xlab) {
      
      # Main function
      main_function <- function(object = object, features.plot = features.plot, xlab = xlab) {
        VlnPlot(object = object, features= features.plot) + 
          labs(x = xlab) + theme(legend.position = 'none')
      }
      
      # Apply main function on all features
      p <- lapply(X = features.plot, object = object, xlab = xlab,
                  FUN = main_function)
      
      cowplot::plot_grid(plotlist = p, ncol = 2)
    }
    
    VlnPlot_2(object = clus_pbmc(), features.plot = gene_list, xlab = "Clusters")
    
  })
  
  # ------------------------------------------------------------------ Feature Plot ---------------------
  fetr_gene_list<- reactive({
    
    gene_list<-  unlist(strsplit(input$fetr_genes, "\n"))
    if(is.null(gene_list)) return()
    rev(gene_list)
    
  })
  
  output$fetr_plot<- renderPlot({
    
    if(is.null(clus_pbmc())) return() 
    
    gene_list<- fetr_gene_list()
    
    #print(gene_list)
    
    FeaturePlot(clus_pbmc(), features = gene_list, ncol = 2, label  = T) 
    
  })
  # ----------------------------------------------------------------- Pathways Analysis --------------------- 
  
  input_paths<- reactive({
    
    #if(is.null(cluster_markers())) return() 
    #else if(is.null(scluster_markers())) return()
    if(input$path_submit == FALSE) return()
    spc_code<- switch(input$species_name, species_human = 'hsa', species_mouse = 'mmu')
    isolate({
      if(input$path_source == "path_KEGG")
      {
        
        kegg_pathways_t <- mappedkeys(KEGGPATHID2EXTID) %>% str_subset(spc_code) # maps KEGG pathway identifiers to Entrez Gene; 'hsa' for human and 'mmu' for mouse
        #kegg_pathways_t
        
        kegg_pathways<- kegg_pathways_t%>% str_replace(sprintf('^%s(.*)$', spc_code), '\\1')
        #kegg_pathways
        
        #............Error: not a graph object
        # It may be that you have both igraph and sna packages loaded at the same time.  If you detach the igraph package your code should run fine.
        
        #detach("package:igraph")
        
        kegg_names <- KEGGPATHID2NAME[kegg_pathways] %>% as.list %>% simplify # maps KEGG pathway identifiers to KEGG pathway names
        #kegg_names
        
        kegg_to_kegg_path <- KEGGPATHID2EXTID[kegg_pathways %>% {
          sprintf('%s%s', spc_code, .)
        }] %>% as.list
        
        #kegg_to_kegg_path
        
        names(kegg_to_kegg_path) <- kegg_names # replce KEGG pathway identifiers by KEGG pathways names
        #kegg_to_kegg_path
        #length(x = kegg_to_kegg_path) # 229
        
        # .... convert entrez gene id to gene Symbol
        # kegg_to_symbol_path <- lapply(kegg_to_kegg_path,function(x){ 
        #  syms=eg2id(x, org="Hs", category ="symbol")
        # return(syms[,2]) })
        
        #kegg_to_symbol_path
        
        kegg_to_kegg_path
        
      }
      
      else if(input$path_source == "path_GO")
      {
        
        #data(go.gs)
        spc_code_go<- switch(input$species_name, species_human = 'human', species_mouse = 'mouse')
        go.hs=go.gsets(species=spc_code_go)
        
        go.bp=go.hs$go.sets[go.hs$go.subs$BP] # Biological Process
        #length(x = go.bp) # 15,927
        #go.bp
        
        # .... convert entrez gene id to gene Symbol
        #go.gs.sym <- lapply(go.bp,function(x){ 
        # syms=eg2id(x, org="Hs", category ="symbol")
        # return(syms[,2]) })
        #go.gs.sym
        
        go.bp
        
        
      }
    })
  })
  
  
  output_paths <-  reactive({
    
    if(is.null(input_paths())) return()
    spc_code_sym<- switch(input$species_name, species_human = 'Hs', species_mouse = 'Mm')
    isolate({
      # get gene symbols; mouse = 'Mm', human = 'Hs'
      SYMBOL2EG <-
        eval(parse(text = sprintf(
          'org.%s.egSYMBOL2EG', spc_code_sym
        )))
      
      if(input$path_input_opt == "path_clus")
      {
        tde<- scluster_markers()
      }
      
      else if(input$path_input_opt == "path_clus_vs_clus")
      {
        tde<- cluster_markers()
        
        
      }
      
      
      logFC_score <- tde$avg_logFC
      #logFC_score
      de_genes<- rownames(tde)
      #length(x = de_genes)
      names(logFC_score) <- de_genes
      #logFC_score
      
      #if(input$path_source == "path_KEGG")
      #....{
      genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
      #length(x = genes) 
      
      
      logFC_score<- logFC_score[genes] # access logFC score 
      #length(x = logFC_score) 
      
      gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% simplify
      #length(x = gene_entrez)
      
      names(logFC_score) <- gene_entrez
      #barplot(sort(logFC_score, decreasing = T))
      #length(x = logFC_score) 
      
      # dput(genes)
      #....}
      
      
      
      kegg_to_symbol_path<- input_paths() # input pathways
      
      fgseaRes1 <- fgsea(kegg_to_symbol_path, logFC_score,
                         minSize=1,
                         maxSize=Inf,
                         nperm = 10000) # run fgsea
      fgseaRes1
      
    })
    
  })
  
  # show plthways plot
  output$path_show<- renderPlot({
    if(is.null(output_paths())) return(NULL)
    else if(input$path_top_plot == FALSE) return(NULL) 
    
    #isolate({
    if(input$path_select == "path_pos"){
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES>=0, padj<=input$path_qval)
      tt<- ": Upregulated"
    }
    
    else if(input$path_select == "path_neg")  
    {
      
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES<0,padj<=input$path_qval)
      tt<- ": Downregulated"
    }
    #}) # end isolation
    
    # select top paths
    ggdat<- ftr_res %>% head(input$path_top) %>% mutate(pathway =fct_inorder(pathway))
    
    isolate({
      if(input$path_source == "path_GO")
      {
        tt0<- "GO (Biological Process)Pathways"
        
      }
      
      else if(input$path_source == "path_KEGG")
      {
        
        tt0<- "KEGG Pathways"
        
      }
      
    }) # End isolate 
    
    ggplot(ggdat) +
      geom_point(aes(
        x = pathway,
        y = -log10(padj),
        #size = size,
        colour = size
      ), size = 8) + #geom_text(aes(x = pathway, y = NES, label = size, colour = size),
      # vjust = "inward", hjust = "inward",
      # show.legend = FALSE, size = 5, angle = 45)+
      labs(title = paste(tt0, tt),
           x = 'Pathways',
           y = '-log10(padj)') +
      #scale_size_continuous(name = 'Size of\nthe pathway',range = c(2,8)) +
      theme_grey(base_size =14 ) +
      theme(axis.text.x = element_text(angle = -23, hjust = 0, size = 11, color  = "Black"),
            plot.margin = margin(10,10,5,5)) + scale_color_gradient("Size of\nthe pathway" , low="blue", high="red")
    
    
  })
  
  # Show Path list
  output$path_list_table = DT::renderDataTable({
    if(is.null(output_paths())) return()
    else if(input$path_list_show == FALSE) return() 
    
    if(input$path_select == "path_pos"){
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES>=0, padj<=input$path_qval)
    }
    
    else if(input$path_select == "path_neg")  
    {
      
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES<0, padj<=input$path_qval)
      
    }
    
    t_path<- length(rownames(ftr_res))
    slist<- c()
    for (i in 1:t_path) {
      slist[i] = length(ftr_res$leadingEdge[[i]])
    }
    
    ftr_res$leadingEdge <- slist
    
    ftr_res <- ftr_res %>% 
      dplyr::  rename(
        Sig_genesets = leadingEdge
      )
    
    ftr_res
  })
  
  
  output$path_download_UI <- renderUI({
    
    if(!is.null(output_paths())){
      tags$head(
        tags$style(HTML('#path_all_download{background-color:#800080;color:white;}'))
      )
      downloadButton("path_all_download", "Download All Pathways as a csv")
    }
  })
  
  # Download All Pathways list .........
  output$path_all_download <- downloadHandler(
    if(is.null(output_paths())) return(),
    filename = function() {paste("All_Pathways_list", ".csv", sep = "")},
    content = function(file) {
      tm<- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(-NES)
      write.csv(tm[,1:7], file, row.names = TRUE)}
  )
  
  #  download filtered datasets .......
  output$path_ftr_download_UI <- renderUI({
    
    if(!is.null(output_paths())){
      tags$head(
        tags$style(HTML('#path_ftr_download{background-color:#800080;color:white;}'))
      )
      downloadButton("path_ftr_download", "Download Filtered Pathways as a csv")
    }
  })
  
  
  output$path_ftr_download <- downloadHandler(
    if(is.null(output_paths())) return(),
    filename = function() {paste("Filtered_Pathways_list", ".csv", sep = "")},
    content = function(file) {
      
      if(input$path_select == "path_pos"){
        ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
        ftr_res<- dplyr::filter(ggdat_t,NES>=0, padj<=input$path_qval)
      }
      
      else if(input$path_select == "path_neg")  
      {
        
        ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
        ftr_res<- dplyr::filter(ggdat_t,NES<0, padj<=input$path_qval)
        
      }
      
      write.csv(ftr_res[,1:7], file, row.names = TRUE)
    }
    
    
  )
  
  # .....................................access specific Pathway gene set
  
  path_genes<-reactive({
    if(is.null(output_paths())) return()
    else if(input$path_genesets==FALSE) return()
    
    
    if(input$path_select == "path_pos"){
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES>=0, padj<=input$path_qval)
    }
    
    else if(input$path_select == "path_neg")  
    {
      
      ggdat_t <- output_paths() %>% as.data.frame %>% as_data_frame %>% arrange(padj)
      ftr_res<- dplyr::filter(ggdat_t,NES<0, padj<=input$path_qval)
      
    }
    
    
    if(input$path_input_opt == "path_clus")
    {
      tde_data<- scluster_markers()
    }
    
    else if(input$path_input_opt == "path_clus_vs_clus")
    {
      tde_data<- cluster_markers()
      
    }
    
    tmp<- ftr_res[input$path_no,]
    tmp<- data.frame(tmp)
    row.names(tmp) <- tmp$pathway
    
    path_title<<- rownames(tmp)
    
    spc_code_sym<- switch(input$species_name, species_human = 'Hs', species_mouse = 'Mm')
    
    if(input$path_hmap_genes == "path_hmap_genes_sig")
    {
      
      tlp<- ftr_res$leadingEdge[input$path_no]
      
      #length(tlp[[1]])
      
      syms_sig=eg2id(tlp[[1]], org= spc_code_sym, category ="symbol") # convert gene entrez ID to gene symbol 
      #length(syms_sig[,2])
      
      #print(syms_sig[,2])
      syms_sig[,2]
    }
    
    
    else if(input$path_hmap_genes == "path_hmap_genes_all")
    {
      
      #print(path_title)
      # convert gene entrez ID to gene symbol 
      id_to_sym <- lapply(input_paths()[rownames(tmp)],function(x){ 
        syms=eg2id(x, org= spc_code_sym, category ="symbol")
        return(syms[,2])})
      
      
      t_genes<- intersect(rownames(tde_data), id_to_sym[[1]])
      t_genes
      #path_genes<- tde_data[id_to_sym[[1]],]
      #dim(path_genes)
      #ompleterecords <- na.omit(path_genes) # remove NA rows
      # dim(ompleterecords)
      #t_genes<- rownames(ompleterecords)
      #t_genes<- rownames(GetAssayData(object = clus_pbmc()[path_genes], slot = "scale.data"))
      #length(x  = t_genes)
      
      #print(t_genes)
      
    }
    
  })
  
  # .....................................print specific Pathway gene set
  output$path_gene_list<-renderPrint({
    
    if(is.null(output_paths())) return()
    else if(input$path_genesets==FALSE) return()
    if(is.null(path_genes())) return()
    
    print(path_genes())
    
  })
  
  output$path_hmap_show<- renderPlot({
    if(is.null(path_genes())) return()
    else if(input$path_hmap==FALSE) return()
    
    DoHeatmap(clus_pbmc(), features = path_genes() , group.bar = T) + 
      ggtitle(paste("Pathway Name:", path_title))
  })
  
  
  #--------------------------------------------------------------------- Dynamic Heatmap -------------------------------------
  #output$dhmap_input_UI<- renderUI({
    
   # if(input$dhmap_input_opt == "dhmap_file_inp")
    #{
      
     # fileInput("dhmap_file", "Input Gene list as a text file:")
      
    #}
    
    #else if (input$dhmap_input_opt == "dhmap_write_inp")
    #{
     # textAreaInput("dh_map_genelist","Write Gene List", width = "200px", height = "150px", placeholder = "Please write gene sybmols: one gene per line")
      
    #}
    
  #})  
  
  # upload gene list processing
  #dhmap_data<-reactive({
   # if (input$dhmap_input_opt == "dhmap_file_inp")
    #{
     # inFile<-input$dhmap_file
      #if(is.null(inFile)) return()
      #conn2 <- file(inFile$datapath,open="r")
      #linn <-readLines(conn2)
      #close(conn2)
      #linn
    #}
    
    #else if(input$dhmap_input_opt == "dhmap_write_inp")
    #{
      
     # gene_list<-  unlist(strsplit(input$dh_map_genelist, "\n"))
      #if(is.null(gene_list)) return()
      #gene_list
    #}
    
  #})
  
  
  #output$dh_plot <- renderPlotly({
    
   # if(is.null(dhmap_data())) return()
    #else if(dhmap_data() == 10) return()
    #isolate({
    #withProgress(message = 'Calculation in progress',
    # detail = 'This may take a while...', value = 0, {
    #for (i in 1:15) {
    # incProgress(1/15)
    # Sys.sleep(0.25)
    # }
    #})
    
    #distfun_row = function(x) dist(x, method = input$distFun_row)
    #distfun_col =  function(x) dist(x, method = input$distFun_col)
    
    #hclustfun_row = function(x) hclust(x, method = input$hclustFun_row)
    #hclustfun_col = function(x) hclust(x, method = input$hclustFun_col)
    
    #sub_genes <- dhmap_data() # gene names 
    #rev_gene_list<-rev(sub_genes) # reversed gene list 
    #m_sub  <- GetAssayData(object = scdata(), slot = "scale.data")[rev_gene_list,] # used scaled data
    
    #heatmaply(m_sub,seriate = "mean", 
    #          row_dend_left = T, 
     #         plot_method = "plotly",showticklabels = c(F,T),labCol = colnames('NA'),
      #        dendrogram = input$dtype,
       #       k_col = input$c, 
        #      k_row = input$r,
         #     distfun_row =  distfun_row,
          #    distfun_col =  distfun_col,
           #   branches_lwd = input$branches_lwd, 
            #  hclustfun_col = hclustfun_col,
             # hclustfun_row = hclustfun_row)
    
    #sub_genes <- dhmap_data()
    #})
    
#  })
  
  #Color Pallete UI ----
  #output$colUI<-renderUI({
   # colSel=ifelse(input$transform_fun=='cor','RdBu','Vidiris')
    #selectizeInput(inputId ="pal", label ="Select Color Palette",
     #              choices = c('Vidiris (Sequential)'="viridis",
      #                         'Magma (Sequential)'="magma",
       #                        'Plasma (Sequential)'="plasma",
        #                       'Inferno (Sequential)'="inferno",
         #                      'Magma (Sequential)'="magma",
          #                     'Magma (Sequential)'="magma",
                               
           #                    'RdBu (Diverging)'="RdBu",
            #                   'RdYlBu (Diverging)'="RdYlBu",
             #                  'RdYlGn (Diverging)'="RdYlGn",
              #                 'BrBG (Diverging)'="BrBG",
               #                'Spectral (Diverging)'="Spectral",
                               
                 #              'BuGn (Sequential)'='BuGn',
                #               'PuBuGn (Sequential)'='PuBuGn',
                  #             'YlOrRd (Sequential)'='YlOrRd',
                    #           'Heat (Sequential)'='heat.colors',
                   #            'Grey (Sequential)'='grey.colors'),
                   #selected=colSel)
  #})
  
  
  
  
  # ----------------------------------------------------- Static Heatmap ---------------------------
  
  output$stmap_input_UI <- renderUI({
    
    if(input$stmap_input_opt == "stmap_file_inp")
    {
      
      fileInput("stmap_file", "Input Gene list as a text file:")
      
    }
    
    else if (input$stmap_input_opt == "stmap_write_inp")
    {
      #textAreaInput("st_map_genelist","Gene List", "Please write gene sybmols: one gene per line")
      textAreaInput("st_map_genelist","Write Gene List", width = "200px", height = "150px", 
                    placeholder = "Please write gene sybmols: one gene per line")
    }
    
  })    
  
  
  # upload gene list processing
  stmap_data<-reactive({
    #isolate({
    #gene_name<- c()
    if (input$stmap_input_opt == "stmap_file_inp")
    {
      inFile<-input$stmap_file
      if(is.null(inFile)) return()
      conn2 <- file(inFile$datapath,open="r")
      linn <-readLines(conn2)
      close(conn2)
      linn
      #rev_linn<-rev(linn)
      #rev_linn
      #sub_genes<- GetAssayData(object = scdata(), slot = "scale.data")[rev_linn,1:input$stmap_cells]
      #sub_genes
      
      
    }
    
    else if(input$stmap_input_opt == "stmap_write_inp")
    {
      
      gene_list<-  unlist(strsplit(input$st_map_genelist, "\n"))
      if(is.null(gene_list)) return()
      gene_list
      #rev_gene_list<-rev(gene_list)
      #rev_gene_list
      #sub_genes<- GetAssayData(object = scdata(), slot = "scale.data")[rev_gene_list,1:input$stmap_cells]
      #sub_genes
      
      
    }
    
    #sub_genes<- GetAssayData(object = scdata(), slot = "scale.data")[gene_name,1:input$stmap_cells]
    #sub_genes<- data.matrix(as.matrix(sub_genes))
    #sub_genes
    
    #}) end isolate 
    
  })
  
  #.... draw static heatmap plot    
  output$sta_hmap<- renderPlot({
    if(is.null(stmap_data())) return() 
    
    
    if(input$stmap_optn == "stmap_optn_1")
      DoHeatmap(clus_pbmc(), features = stmap_data(), group.bar = T) 
    
    else if(input$stmap_optn == 'stmap_optn_2')
    {
      
      sub_genes <- stmap_data() # gene names 
      rev_gene_list<-rev(sub_genes) # reversed gene list 
      m_sub  <- GetAssayData(object = scdata(), slot = "scale.data")[rev_gene_list,] # used scaled data
      
      hplot3<-heatmap3(m_sub,Rowv = NA, Colv = m_sub, showRowDendro = F, showColDendro = F, labCol  = NA,
                       lasRow = 1, margins = c(5,5) )
      hplot3
    }
    
    
    
    
    
    
    #.........................heatmap.2 ....................
    #my_palette <- colorRampPalette(c("forestgreen", "yellow", "red"))(n = 299)
    #col_breaks = c(seq(-1,-0.5,length=100),  # forestgreen
    # seq(-0.5,0.5,length=100), # yellow
    #seq(0.5,1,length=100))
    # distance= dist(m_sub, method ="euclidean")    
    #hcluster = hclust(distance, method ="ward.D")
    #dend1 <- as.dendrogram(hcluster)
    #cols_branches <- c("darkred", "forestgreen", "orange", "blue")
    #dend1 <- color_branches(dend1, k = 4, col = cols_branches)
    #col_labels <- get_leaves_branches_col(dend1)
    #col_labels <- col_labels[order(order.dendrogram(dend1))]
    
    #heatmap.2(m_sub, trace = "none", col = viridis(100),
    #         RowSideColors=col_labels, key = T, 
    #          colRow = col_labels,Rowv = dend1, labCol = TRUE,cexRow = 1.3)
    
    
  }) # end static heatmap  plot
  
  # ---------------------------------------------------- Dynamic PCA Plot ---------------------------------------------
  
  output$dpcaplot <- renderPlotly({
    
    inFile<-input$dpca_file
    if(is.null(inFile)) return()
    conn2 <- file(inFile$datapath,open="r")
    linn <-readLines(conn2)
    close(conn2)
    
    m_sub  <- GetAssayData(object = scdata(), slot = "scale.data")[linn,] # used scaled data
    tm_sub<- t(m_sub)
    colour1 = rainbow(length(rownames(tm_sub)))
    pre<-prcomp(tm_sub)
    autoplotly(pre, data = tm_sub
               , colour = colour1, size = 3)
  })
  
  # T-sne plot
  #output$tsnep<- renderPlot({
  
  # inFile<-input$tsne_file
  #if(is.null(inFile)) return()
  #conn2 <- file(inFile$datapath,open="r")
  #gene_lst <-readLines(conn2)
  #close(conn2)
  
  #m_sub<- GetAssayData(object = scdata(), slot = "scale.data")[gene_lst,] # used scaled data
  #tm_sub<- t(m_sub)
  #colour1 = rainbow(length(rownames(tm_sub)))
  #pt<-Rtsne(tm_sub,dims=2, perplexity=1,check_duplicates = F)
  #tsne_plot<- data.frame(x = pt$Y[,1], y = pt$Y[,2], col = colour1)
  #ggplot(tsne_plot) + geom_point(aes(x=x, y=y), col=colour1, size = 3)
  #})
  
  
  
  # load TSNE image ...................
  #output$tsnep <- renderImage({
  
  #list(src = filenames)
  
  #})
  
  # --------------------------------------------------- Correlation Plot ------------------------------------------
  output$corplot<- renderPlot({
    if(is.null(scdata())) return()
    else if(is.null(input$gene2)) return()
    
    #cord<-hv_pbmc()@scale.data[c(input$gene1,input$gene2),]
    cord<-GetAssayData(object = scdata(), slot = "scale.data")[c(input$gene1,input$gene2),] # used scaled data
    cordt<-t(as.matrix(cord))
    mat_cord<-data.frame(cordt)
    
    #if(input$cor_impute == 1){
    
    #cor_imputed_data<-AddImputedScore(object = hv_pbmc(), genes.use = rownames(hv_pbmc()@scale.data),
    #genes.fit = c(input$gene1,input$gene2), do.print = TRUE, gram=FALSE)
    
    #cordt<-t(as.matrix(cor_imputed_data@imputed[c(input$gene1,input$gene2),]))
    #mat_cord<-NULL
    #mat_cord<-data.frame(cordt)
    # MAGIC_data <- magic(hv_pbmc(), genes=c(input$gene1,input$gene2))
    # cordt<-t(as.matrix(MAGIC_data@data))
    # mat_cord<-data.frame(cordt)
    
    #}
    
    
    if(input$cor_log == 1)   # apply log 2
    {
      
      mat_cord<-log2(mat_cord+0.01)
    }
    
    # if(input$zeroexp == 1) # excluding cells which zero values in both genes
    # {
    # p=rowSums(cordt>0)>=1
    # nv=cordt[p,]
    # mat_cord<-data.frame(nv)
    # }
    
    ggscatter(mat_cord, x = input$gene1, y = input$gene2, 
              add =input$regline, conf.int = TRUE, 
              cor.coef = TRUE, cor.method = input$corm,
              xlab = input$gene1, ylab = input$gene2, color = "black",
              cor.coef.size = 4,shape = 19, size  = 2,
              add.params = list(color = "blue", fill = "lightgray"))
    
  })
  # ------------------------------------------------------ Trajectory Analysis ------------------------------------------------------------------
  
  output$tray_input_UI <- renderUI({
    
    if(input$tray_input_opt == "tray_file_inp")
    {
      
      fileInput("tray_file", "Input Gene list as a text file:")
      
    }
    
    else if (input$tray_input_opt == "tray_write_inp")
    {
      #textAreaInput("tray_genelist","Gene List", "Please write gene sybmols: one gene per line")
      textAreaInput("tray_genelist","Write Gene List", width = "200px", height = "150px", 
                    placeholder = "Please write gene sybmols: one gene per line")
      
    }
    
  })    
  
  
  tray_data<-reactive({
    if(is.null(clus_pbmc())) return()
    #else if (input$de_ok == FALSE) return()
    
    #isolate({
    data <- as(as.matrix(clus_pbmc()@assays$RNA@scale.data), 'sparseMatrix')
    #dim(x = data)
    #data[1:2,1:3]
    pd <- new('AnnotatedDataFrame', data = clus_pbmc()@meta.data)
    #pd[1:5,]
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    #fData
    
    fd <- new('AnnotatedDataFrame', data = fData)
    #fd[1:2,]
    
    #Construct monocle cds
    HSMM2 <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            #lowerDetectionLimit = 0.5,
                            expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
    
    #View data
    #pData(HSMM)
    #fData(HSMM)
    
    # Select genes and Run ordering algorithm
    if(input$tray_genes == "tray_hvg"){
      ordering_genes<- clus_pbmc()[["RNA"]]@var.features
    }
    
    else if(input$tray_genes == "tray_all")
    {
      ordering_genes<- rownames(x = clus_pbmc()[["RNA"]])
    }
    
    HSMM2 <- setOrderingFilter(HSMM2, ordering_genes)
    #print(dim(exprs(HSMM)))
    
    ## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
    HSMM2 <- reduceDimension(HSMM2,norm_method="none", 
                             reduction_method="DDRTree",
                             max_components=4,
                             scaling=TRUE,
                             verbose=TRUE,
                             pseudo_expr=0)
    
    HSMM2 <- orderCells(HSMM2)
    HSMM2
    #})
  })
  
  # .....trajectory plot
  output$tray_plot<-renderPlot({
    
    if(is.null(tray_data())) return()
    else if(input$tray_plot_set == F) return()
    
    #isolate({
    try2_data<- tray_data()
    if(input$tray_plot_set_opts == "tray_clst_ck")
      plot_opt<- "seurat_clusters"
    else if(input$tray_plot_set_opts == "pseudotime_ck")
      plot_opt<- "Pseudotime"
    else if(input$tray_plot_set_opts == "tray_stat_ck")
      plot_opt<- "State"
    
    plot_cell_trajectory(try2_data, 
                         color_by = plot_opt,
                         theta = -15,
                         show_branch_points = T,
                         show_tree = TRUE,
                         cell_size = 1.5)
    
    #})
    
  })
  # upload gene list processing for pseudotemporal plot
  tray_input_genes<-reactive({
    if(input$tray_pseu_set == F) return()
    
    if (input$tray_input_opt == "tray_file_inp")
    {
      inFile<-input$tray_file
      if(is.null(inFile)) return()
      conn2 <- file(inFile$datapath,open="r")
      linn <-readLines(conn2)
      close(conn2)
      linn
    }
    
    else if(input$tray_input_opt == "tray_write_inp")
    {
      
      gene_list<-  unlist(strsplit(input$tray_genelist, "\n"))
      if(is.null(gene_list)) return()
      gene_list
      #print(gene_list)
      
    }
    
  })
  
  # .....pseudotemporal plot
  output$tray_pseu_plot<-renderPlot({
    
    if(input$tray_pseu_set== F) return()
    #else if(is.null(tray_input_genes())) return()
    if(input$tray_submit == F) return()
    
    #print(tray_input_genes())
    ptdata<- tray_data()
    #genex1<- c("GAPDH", "TP53")
    isolate({
      
      withProgress(message = 'Generating plot...', value = 0, {
        
        n <- 1
        
        for (i in 1:n) {
          plot_pseudotime_heatmap(ptdata[tray_input_genes(),],
                                  num_clusters = input$tray_pseu_clst, 
                                  cores = 4,
                                  hmcols = NULL,
                                  show_rownames = T)
          # Increment the progress bar, and update the detail text.
          incProgress(10, detail = paste("Completed ", i*100,"%"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(4)
        }
      })
    })
    
  })
  
  # --------- End Trajectory Analysis
  
  
} # End Server function



# ------------------- ----------------CAll App-------------------------

shinyApp(ui = ui, server = server)
