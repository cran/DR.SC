## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library("DR.SC")

## ----eval=FALSE,message=FALSE, warning=FALSE----------------------------------
#  data("HCC1", package = 'DR.SC')

## ----eval=FALSE---------------------------------------------------------------
#  # standard log-normalization
#  HCC1 <- NormalizeData(HCC1)
#  # choose 2000 highly variable features
#  seu <- FindVariableFeatures(HCC1, nfeatures = 1000)
#  

## ----eval=FALSE---------------------------------------------------------------
#  ### Given K
#  seu <- DR.SC(seu, K=6, platform = 'Visium', variable.type = 'HVGs',verbose=F)

## ----eval=FALSE---------------------------------------------------------------
#  spatialPlotClusters(seu)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seu)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seu, visu.method = 'UMAP')

## ----eval=FALSE---------------------------------------------------------------
#  # choose 2000 spatially variable features
#  seus <- FindSVGs(seu, nfeatures = 1000)

## ----eval=FALSE---------------------------------------------------------------
#  ### Given K
#  seus <- DR.SC(seus, K=6, platform = 'Visium', variable.type='SVGs', verbose=T)

## ----eval=FALSE---------------------------------------------------------------
#  spatialPlotClusters(seus)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seus)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seus, visu.method = 'UMAP')

## ----eval=FALSE---------------------------------------------------------------
#  dat <- FindAllMarkers(seus)
#  head(dat)
#  library(dplyr, verbose=F)
#  top2 <-  dat %>%
#    group_by(cluster) %>%
#    top_n(n = 2, wt = avg_log2FC)
#  top2

## ----eval=FALSE---------------------------------------------------------------
#  genes <- top2$gene[seq(1, 12, by=2)]
#  RidgePlot(seus, features = genes, ncol = 2)

## ----eval=FALSE---------------------------------------------------------------
#  
#  VlnPlot(seus, features = genes, ncol=2)

## ----eval=FALSE---------------------------------------------------------------
#  seus <- RunTSNE(seus, reduction="dr-sc", reduction.key='drsc_tSNE_')
#  FeaturePlot(seus, features = genes, reduction = 'tsne' ,ncol=2)
#  

## ----eval=FALSE---------------------------------------------------------------
#  DotPlot(seus, features = genes)

## ----eval=FALSE---------------------------------------------------------------
#  # standard scaling (no regression)
#  seus <- ScaleData(seus)
#  DoHeatmap(subset(seus, downsample = 500), features = genes, size = 5)

## ----eval=FALSE---------------------------------------------------------------
#  # choose 2000 spatially variable features
#  seus <- FindSVGs(seu, nfeatures = 2000, verbose = F)

## ----eval=FALSE---------------------------------------------------------------
#  ### Given K
#  seus <- DR.SC(seus, K_set=3:8, platform = 'Visium', variable.type='SVGs', verbose=F)

## ----eval=FALSE---------------------------------------------------------------
#  mbicPlot(seus)

## ----eval=FALSE---------------------------------------------------------------
#  spatialPlotClusters(seus)

## ----eval=FALSE---------------------------------------------------------------
#  drscPlot(seus, dims=1:10)

## ----eval=FALSE---------------------------------------------------------------
#  sessionInfo()

