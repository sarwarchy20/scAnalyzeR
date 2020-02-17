
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

installIfNeeded = function (packages, installFn = install.packages) {
  newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(newPackages)) installFn(newPackages)
}

cranPackages = c(
  "devtools", "Rtsne", "igraph", "visNetwork", "shiny", "htmltools", "shinyjs", "htmlwidgets", "plotly", "ggpubr", "autoplotly",
  "dplyr","viridis", "RColorBrewer", "dendextend", "plotly","tsne", "gplots","ggplot2" ,"ggpubr", "ggfortify", "dplyr","heatmap3","DT","heatmaply",
  "tidyverse","Matrix", "cowplot")

installIfNeeded(cranPackages)

if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")

if (!requireNamespace("fgsea", quietly = TRUE))
  BiocManager::install("fgsea")

if (!requireNamespace("KEGG.db", quietly = TRUE))
  BiocManager::install("KEGG.db")

if (!requireNamespace("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")

if (!requireNamespace("GO.db", quietly = TRUE))
BiocManager::install("GO.db")

if (!requireNamespace("pathview", quietly = TRUE))
  BiocManager::install("pathview")

if (!requireNamespace("gage", quietly = TRUE))
  BiocManager::install("gage")

if (!requireNamespace("monocle", quietly = TRUE))
  BiocManager::install("monocle")

if (!requireNamespace("MAST", quietly = TRUE))
  BiocManager::install("MAST")
