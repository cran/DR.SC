## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library("DR.SC")

## ----eval=FALSE---------------------------------------------------------------
#  seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
#  head(seu@meta.data)

## ----eval=FALSE---------------------------------------------------------------
#  ### Given K
#  library(Seurat)
#  seu <- NormalizeData(seu)
#  # choose 2000 variable features using Seurat
#  seu <- FindVariableFeatures(seu, nfeatures = 2000)
#  seu2 <- DR.SC(seu, K=4, platform = 'ST', verbose=F)

## ----eval=FALSE---------------------------------------------------------------
#  mclust::adjustedRandIndex(seu2$spatial.drsc.cluster, seu$true_clusters)

## ----eval=FALSE---------------------------------------------------------------
#  spatialPlotClusters(seu2)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seu2)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seu2, visu.method = 'UMAP')

## ----eval=FALSE---------------------------------------------------------------
#  seu2 <- DR.SC(seu, q=10, K=NULL, K_set =2:6, platform = 'ST', verbose=F)
#  mbicPlot(seu2)

## ----eval=FALSE---------------------------------------------------------------
#  genes <- c("gene-24","gene-68", "gene-95","gene-55")
#  RidgePlot(seu2, features = genes, ncol = 2)

## ----eval=FALSE---------------------------------------------------------------
#  
#  VlnPlot(seu2, features = genes, ncol=2)

## ----eval=FALSE---------------------------------------------------------------
#  seu2 <- RunTSNE(seu2, reduction="dr-sc", reduction.key='drsc_tSNE_')
#  FeaturePlot(seu2, features = genes, reduction = 'tsne' ,ncol=2)
#  

## ----eval=FALSE---------------------------------------------------------------
#  DotPlot(seu2, features = genes)

## ----eval=FALSE---------------------------------------------------------------
#  # standard scaling (no regression)
#  seu2 <- ScaleData(seu2)
#  DoHeatmap(subset(seu2, downsample = 500), features = genes, size = 5)

## ----eval=FALSE---------------------------------------------------------------
#  sessionInfo()

