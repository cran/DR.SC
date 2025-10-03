## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2024)
#  library("DR.SC")

## ----eval = FALSE,message=FALSE, warning=FALSE--------------------------------
#  data("dlpfc151510", package = 'DR.SC')
#  

## ----eval = FALSE,message=FALSE, warning=FALSE--------------------------------
#  library(Seurat)
#  count <- dlpfc151510@assays$RNA@counts
#  meta_data <- data.frame(row=dlpfc151510@meta.data$row, col=dlpfc151510@meta.data$col, annotation=dlpfc151510$annotation)
#  row.names(meta_data) <- colnames(count)
#  ## create Seurat object
#  dlpfc151510 <- CreateSeuratObject(counts=count, meta.data = meta_data)
#  head(dlpfc151510)

## ----eval = FALSE-------------------------------------------------------------
#  
#  # standard log-normalization
#  dlpfc151510 <- NormalizeData(dlpfc151510, verbose = F)
#  # choose 500 highly variable features
#  seu <- FindVariableFeatures(dlpfc151510, nfeatures = 500, verbose = F)
#  VariableFeatures(seu)[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  seu <- DR.SC(seu, K=7, platform = 'Visium', verbose=T, approxPCA=T)

## ----eval = FALSE-------------------------------------------------------------
#  mclust::adjustedRandIndex(seu$spatial.drsc.cluster, seu$annotation)
#  spatialPlotClusters(seu)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seu, visu.method = 'UMAP')

## ----eval = FALSE-------------------------------------------------------------
#  # choose 480 spatially variable features
#  seus <- FindSVGs(seu, nfeatures = 480)
#  VariableFeatures(seu)[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  seus <- DR.SC(seus,  K=7, platform = 'Visium', verbose=F, approxPCA=T)

## ----eval = FALSE-------------------------------------------------------------
#  spatialPlotClusters(seus)
#  mclust::adjustedRandIndex(seus$spatial.drsc.cluster, seus$annotation)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seus)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seus, visu.method = 'UMAP')

## ----eval = FALSE-------------------------------------------------------------
#  SVGs <- topSVGs(seus, ntop = 400)
#  dat <- FindAllMarkers(seus, features = SVGs)
#  head(dat)
#  library(dplyr, verbose=F)
#  top2 <-  dat %>%
#    group_by(cluster) %>%
#    top_n(n = 2, wt = avg_log2FC)
#  top2

## ----eval = FALSE, fig.height=6, fig.width=10---------------------------------
#  genes <- top2$gene[seq(1, 12, by=2)]
#  RidgePlot(seus, features = genes, ncol = 2)

## ----eval = FALSE, fig.height=8, fig.width=10---------------------------------
#  
#  VlnPlot(seus, features = genes, ncol=2)

## ----eval = FALSE, fig.height=8, fig.width=10---------------------------------
#  seus <- RunTSNE(seus, reduction="dr-sc", reduction.key='drsc_tSNE_')
#  FeaturePlot(seus, features = genes, reduction = 'tsne' ,ncol=2)
#  

## ----eval = FALSE-------------------------------------------------------------
#  DotPlot(seus, features = genes)

## ----eval = FALSE, fig.height=8, fig.width=10---------------------------------
#  top20 <-  dat %>%
#    group_by(cluster) %>%
#    top_n(n = 20, wt = avg_log2FC)
#  genes <- top20$gene
#  # standard scaling (no regression)
#  seus <- ScaleData(seus)
#  DoHeatmap(subset(seus, downsample = 500), features = genes, size = 5)

## ----eval = FALSE-------------------------------------------------------------
#  # choose spatially variable features
#  seus <- FindSVGs(seu, nfeatures = 480, verbose = F)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  seus <- DR.SC(seus, K=3:9, platform = 'Visium', verbose=F)

## ----eval = FALSE-------------------------------------------------------------
#  seus <- selectModel(seus, pen.const = 0.8)
#  mbicPlot(seus)

## ----eval = FALSE-------------------------------------------------------------
#  spatialPlotClusters(seus)

## ----eval = FALSE-------------------------------------------------------------
#  drscPlot(seus, dims=1:10)

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()

