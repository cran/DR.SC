## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  library(DR.SC)
#  seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
#  head(seu@meta.data)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  library(Seurat)
#  seu <- NormalizeData(seu)
#  # choose 2000 variable features using Seurat
#  seu <- FindVariableFeatures(seu, nfeatures = 400)
#  seu2 <- DR.SC(seu, K=4, platform = 'ST',  verbose=F)

## ----eval = FALSE-------------------------------------------------------------
#  mclust::adjustedRandIndex(seu2$spatial.drsc.cluster, seu$true_clusters)

## ----eval = FALSE-------------------------------------------------------------
#  spatialPlotClusters(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu2, visu.method = 'UMAP')

## ----eval = FALSE-------------------------------------------------------------
#  seu2 <- DR.SC(seu, q=10, K=2:6, platform = 'ST', verbose=F)
#  mbicPlot(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  seu <- NormalizeData(seu, verbose=F)
#  # choose 400 spatially variable features using FindSVGs
#  seus <- FindSVGs(seu, nfeatures = 400, verbose = F)
#  seu2 <- DR.SC(seus, K=4, platform = 'ST', verbose=F)

## ----eval = FALSE-------------------------------------------------------------
#  mclust::adjustedRandIndex(seu2$spatial.drsc.cluster, seu$true_clusters)

## ----eval = FALSE-------------------------------------------------------------
#  spatialPlotClusters(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu2, visu.method = 'UMAP')

## ----eval = FALSE-------------------------------------------------------------
#  seu2 <- DR.SC(seus, q=10, K=2:6, platform = 'ST',  verbose=F)
#  mbicPlot(seu2)
#  # or plot BIC or AIC
#  # mbicPlot(seu2, criteria = 'BIC')
#  # mbicPlot(seu2, criteria = 'AIC')
#  # tune pen.const
#  seu2 <- selectModel(seu2, pen.const = 0.7)
#  mbicPlot(seu2)

## ----eval = FALSE-------------------------------------------------------------
#  dat <- FindAllMarkers(seu2)
#  suppressPackageStartupMessages(library(dplyr) )
#  # Find the top 1 marker genes, user can change n to access more marker genes
#  dat %>%group_by(cluster) %>%
#      top_n(n = 1, wt = avg_log2FC) -> top
#  genes <- top$gene
#  RidgePlot(seu2, features = genes, ncol = 2)

## ----eval = FALSE-------------------------------------------------------------
#  VlnPlot(seu2, features = genes, ncol=2)

## ----eval = FALSE-------------------------------------------------------------
#  seu2 <- RunTSNE(seu2, reduction="dr-sc", reduction.key='drsctSNE_')
#  FeaturePlot(seu2, features = genes, reduction = 'tsne' ,ncol=2)
#  

## ----eval = FALSE-------------------------------------------------------------
#  DotPlot(seu2, features = genes)

## ----eval = FALSE-------------------------------------------------------------
#  # standard scaling (no regression)
#  dat %>%group_by(cluster) %>%
#      top_n(n = 30, wt = avg_log2FC) -> top
#  ### select the marker genes that are also the variable genes.
#  
#  genes <- intersect(top$gene, seu2[['RNA']]@var.features)
#  ## Change the HVGs to SVGs
#  #  <- topSVGs(seu2, 400)
#  seu2 <- ScaleData(seu2, verbose = F)
#  DoHeatmap(subset(seu2, downsample = 500),features = genes, size = 5)

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()

