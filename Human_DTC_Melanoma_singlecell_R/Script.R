#  Installing and loading packages
###################################
install.packages("Seurat")
install.packages("ggpubr")
install.packages("Matrix")
install.packages("tidyverse")
install.packages("hdf5r")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

install.packages("devtools")
devtools_version <- packageDescription("devtools", fields = "Version")
if (package_version(devtools_version) < "1.6") {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval", force = TRUE)
devtools::install_github("hadley/dplyr", force = TRUE)

library(devtools)
library("ggplot2")
library("Seurat")
library(ggpubr)
library(Matrix)
library(tidyverse)
library(dplyr)
library(GEOquery)
library(hdf5r)
##################################

path <- getwd()
list.files(path)
set.seed(42)

### Loading data ==================
hdf5_obj <- Read10X_h5(filename = "10k_Human_DTC_Melanoma_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)

# Peak the "feature barcode (sparse) matrix"
hdf5_obj[1:10,1:10]
# 10 x 10 sparse Matrix of class "dgCMatrix"
#  [[ suppressing 10 column names ‘AAACCAAAGAACCAGG-1’, ‘AAACCAAAGAACCATT-1’, ‘AAACCAAAGAACCTCA-1’ ... ]]
#
#  DDX11L2         . . . . . . . . . .
#  MIR1302-2HG     . . . . . . . . . .
#  FAM138A         . . . . . . . . . .
#  ENSG00000290826 . . . . . . . . . .
#  OR4F5           . . . . . . . . . .
#  ENSG00000238009 . . . . . . . . . .
#  ENSG00000239945 . . . . . . . . . .
#  ENSG00000239906 . . . . . . . . . .
#  ENSG00000241860 . . . . . . . . . .
#  ENSG00000241599 . . . . . . . . . .

### Create Seaurat Object ===============================
seuratobj <- CreateSeuratObject(counts = hdf5_obj, 
                                min.cells = 3,
                                min.features = 200)

### QC (Quality Control) =================================
##########################################################
# % MT reads
seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern ="^MT-")
View(seuratobj@meta.data)

# Violin plot visual
VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Scatter plot visual
FeatureScatter(seuratobj, feature1="nFeature_RNA", feature2="nCount_RNA") +
                 geom_smooth(method = 'lm')

### Filtering =============================================
seuratobj <- subset(seuratobj, 
                    subset = nFeature_RNA > 200 
                    & nFeature_RNA < 2500
                    & percent.mt < 5)


### Normalizing ========================================
seuratobj <- NormalizeData(seuratobj)


### Identifying "highly variable features"===============
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000, span = 0.3)

## Top 10 most highly var genes
top10 <- head(VariableFeatures(seuratobj), 10)
#top10
#[1] "FAAH2"          
#[2] "GNLY"           
#[3] "ANXA1"          
#[4] "CCL5"           
#[5] "HOMER1"         
#[6] "TOX"            
#[7] "LHFPL6"         
#[8] "ENSG00000259033"
#[9] "CCSER1"         
#[10] "S100B" 


##  plotting the vairable features
plot_var <- VariableFeaturePlot(seuratobj)
LabelPoints(plot = plot_var, points= top10)


### Scaling  ========================================
all.genes <- rownames(seuratobj)
seuratobj <- ScaleData(seuratobj, features = all.genes)


### Linear Dimension reduction  ========================
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))

# visual. PCA results
print(seuratobj[["pca"]], dims = 1:5, nfeatures =5)
DimHeatmap(seuratobj, dims =1, cells =500, balanced = TRUE)

# determine dim of data with elbow-plot
ElbowPlot(seuratobj)


### Clustering  ========================================
seuratobj <- FindNeighbors(seuratobj, dims = 1:15)

# result understanding
seuratobj <- FindClusters(seuratobj, resolution = c(0.3, 0.7, 1))
"""
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 713
Number of edges: 27836

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
  [----|----|----|----|----|----|----|----|----|----|
     **************************************************|
     Maximum modularity in 10 random starts: 0.8914
   Number of communities: 4
   Elapsed time: 0 seconds
   Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
   
   Number of nodes: 713
   Number of edges: 27836
   
   Running Louvain algorithm...
   0%   10   20   30   40   50   60   70   80   90   100%
     [----|----|----|----|----|----|----|----|----|----|
        **************************************************|
        Maximum modularity in 10 random starts: 0.7655
      Number of communities: 5
      Elapsed time: 0 seconds
      Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
      
      Number of nodes: 713
      Number of edges: 27836
      
      Running Louvain algorithm...
      0%   10   20   30   40   50   60   70   80   90   100%
        [----|----|----|----|----|----|----|----|----|----|
           **************************************************|
           Maximum modularity in 10 random starts: 0.6907
         Number of communities: 6
         Elapsed time: 0 seconds
"""
View(seuratobj@meta.data)

DimPlot(seuratobj, group.by = "RNA_snn_res.0.3", label=TRUE)
#DimPlot(seuratobj, group.by = "RNA_snn_res.0.7", label=TRUE)

# Set identity of clusters
Idents(seuratobj)
Idents(seuratobj) <- "RNA_snn_res.0.7"
Idents(seuratobj)

### UMAP ===============================
seuratobj <- RunUMAP(seuratobj, dims = 1:15)

# visualize umap
DimPlot(seuratobj, reduction = 'umap')
















