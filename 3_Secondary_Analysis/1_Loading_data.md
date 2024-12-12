# 1. Loading library 
Make sure that these libraries are already installed in your conda environment:
```
library(Seurat)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)
library(SummarizedExperiment) #to use SingleR
library(scuttle)  #to use SingleR
```
We can also store the output directory path to facilitate the saving process further   
```
OUT_DIR <- "~/CITEseq_analysis/3_Secondary_Analysis"
```

# 2. Create Seurat object 
## 2.1 Load multiomic data (ADT&RNA) 

In this analysis, I will start directly with the filtered matrix provided by Cell Ranger. I will explain how to start from raw data in another script.
```
wk4 <- Read10X_h5("filtered_feature_bc_matrix.h5")
```
```
Genome matrix has multiple modalities, returning a list of matrices for this genome
```
The message "Genome matrix has multiple modalities, returning a list of matrices for this genome" indicates that the list comprises a Gene Expression matrix and an Antibody Capture matrix, so we can separate these:
```
wk4_gex <- wk4$`Gene Expression` 
wk4_ADT <- wk4$`Antibody Capture`
```
The wk4_gex is a matrix where each column represents a cell and each row represents a gene. It contains 27,998 rows in total. The wk4_ADT is a matrix where each column represents a cell and each row represents an antibody/hashtag. It contains 33 rows in total.
```
> head(wk4_gex)
6 x 13198 sparse Matrix of class "dgCMatrix" #only first 6 rows

> head(wk4_ADT)
6 x 13198 sparse Matrix of class "dgCMatrix" #only first 6 rows
```
Since **wk4_ADT** contains both HTO and ADT information, we can split it:
```
adt_features <- rownames(wk4_ADT)
hto_features <- grep("^HTO", adt_features, value = TRUE)
adt_features_only <- setdiff(adt_features, hto_features)

wk4_HTO <- wk4_ADT[hto_features, ]
wk4_ADT <- wk4_ADT[adt_features_only, ]
```
## 2.2 Setup Seurat object
We can first create a Seurat object for gene expression then add ADT&HTO assays:
```
so_wk4 <- CreateSeuratObject(counts = wk4_gex)
so_wk4[["ADT"]] <- CreateAssayObject(counts = wk4_ADT)
so_wk4[["HTO"]] <- CreateAssayObject(counts = wk4_HTO)
```
We can now start the standard [pre-processing workflow for scRNA](https://satijalab.org/seurat/articles/pbmc3k_tutorial) from Seurat 
