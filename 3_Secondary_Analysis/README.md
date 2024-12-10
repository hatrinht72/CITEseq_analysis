# 1. Loading library 
Make sure that these libraries already implemented in your conda environment 

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

# 2. Create Seurat object 
## 2.1 Load multiomic data (ADT&RNA) 

In this analysis, I will start directly with filtered filtered matrix by cellranger, I will explain the way I would do if I start from raw in another script

```
wk4 <- Read10X_h5("filtered_feature_bc_matrix.h5")
```
```
Genome matrix has multiple modalities, returning a list of matrices for this genome
```
The list compromised Gene Expression matrix and Antibody Capture matrix, so we can separate these
```
wk4_gex <- wk4$`Gene Expression` 
wk4_ADT <- wk4$`Antibody Capture`
```
The **wk4_gex** is a matrix with each column  a cell and each row an gene, we have 27998 rows in overall
The **wk4_ADT** is a matrix with each column a cell and each row an antibody/hashtag, we have 33 rows in overall 
```
> head(wk4_gex)
6 x 13198 sparse Matrix of class "dgCMatrix" #only first 6 rows

> head(wk4_ADT)
6 x 13198 sparse Matrix of class "dgCMatrix" #only first 6 rows
```
Since in **wk4_ADT**, we have both HTO and ADT informations, we can split it
```
adt_features <- rownames(wk4_ADT)
hto_features <- grep("^HTO", adt_features, value = TRUE)
adt_features_only <- setdiff(adt_features, hto_features)

wk4_HTO <- wk4_ADT[hto_features, ]
wk4_ADT <- wk4_ADT[adt_features_only, ]
```
## 2.2 Setup Seurat object
We can first create Seurat object for gene expression then add ADT&HTO assays
```
so_wk4 <- CreateSeuratObject(counts = wk4_gex)
so_wk4[["ADT"]] <- CreateAssayObject(counts = wk4_ADT)
so_wk4[["HTO"]] <- CreateAssayObject(counts = wk4_HTO)
```
# 3. Demultiplexing
Since CITE-sequencing pool all of cells from different samples, we need to demultiplex them based on the HTO tag. They re a function of Seurat to do this. 
```
so_wk4 <- HTODemux(so_wk4, assay = "HTO", positive.quantile = 0.99)
```
Based on their expression in HTO, cells will be classified to 3 categories :
- Doublet : expression of HTO is more than for one specific cell
- Singlet : significant expression of HTO corresponding to only one specific sample
- Negative : show anything significant expression of HTO
These data will be stored in HTO_classification.global
```
table(so_wk4$HTO_classification.global)
#Doublet Negative  Singlet 
# 962     2879     9357 
```
We can also visualize the distribution of each HTO across cells, showing if we are distinct population of cells correspond to each sample 
```
png(paste0(OUT_DIR, "3_1_HTO_demultiplexing.png"), width = 1600, height = 1100, res = 150)
Idents(so_wk4) <- "HTO_maxID"
RidgePlot(so_wk4, assay = "HTO", features = rownames(so_wk4[["HTO"]])[1:6], ncol = 3)
dev.off()
```
