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

# 3. Standard pre-processing workflow
## 3.1 Quality control 
To select cells for further analysis, we can base on several metrics : 
The number of unique genes detected in each cell :
  - Cells with too low counts : empty droplet / low quality cells
  - Cells with too high counts : Doublets / multiplets cells 
Contamination with mitochondrial 

We will detect first mitochondrial genes by a function in Seurat:
```
# Calculate the percentage of mitochondrial genes : "mt" for mice, "MT" for human
so_wk4[["percent.mt"]] <- PercentageFeatureSet(so_wk4, pattern = "^mt-")
head(so_wk4@meta.data)
```
We can then visualize the QC metrics and use these to filter cells
```
png(paste0(OUT_DIR, "3_1_QC_RNA.png"), width = 1200, height = 1000, res = 150)
VlnPlot(so_wk4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
```
![3_1_QC_RNA](https://github.com/user-attachments/assets/2f088b34-c8c2-4cf5-b28d-c90647adfd66)

We can also visualize feature-feature relationship 
```
png(paste0(OUT_DIR, "4_2_QC_RNA.png"), width = 1200, height = 1000, res = 150)
plot1 <- FeatureScatter(so_wk4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so_wk4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
```
![3_2_QC_RNA](https://github.com/user-attachments/assets/79f986a4-8e14-42b8-93f9-60c4052279ac)

Based on the observation, we can filter cell by these criterias:
 - nFeature_RNA: Minimum 200, Maximum 7,500
 - nCount_RNA: Minimum 1,000
 - percent.mt: Maximum 5%
```
so_wk4 <- subset(so_wk4, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & 
                   nCount_RNA > 1000 & percent.mt < 5)
```
## 3.2 Normalizing the data
Normalization is a data preprocessing technique used to transform features in a dataset to a common scale.
This is a process can be considered as 2 steps :
 - Scaling : multiply each UMI count by a specific factor to get all cells the same UMI counts, then we can comparing the concentration of counts across cells, not the absolute counts.
 - Transformation : there are two type of transformation :
     - Simple transformation : using NormalizeData() function, assuming that each cell originally contains the same number of RNA molecules
     - Pearson residuals for transformation : using SCTransform() function (Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures), without the assumption that we mentionned ahead 

In our context, the aim of this step is to normalize the feature expression measurements for each cell by the total expression, then multiplies this by a scale factor.
Here I will use the NormalizeData()
```
so_wk4 <- NormalizeData(so_wk4, normalization.method = "LogNormalize", scale.factor = 10000)
```

## 3.3 Scaling the data
This is a standard pre-processing step prior to any reduction technique
with the ScaleData() function in Seurat, we can :
- Shifts the expression of each gene, so that the mean expression across cells is 0
- Scales the expression of each gene, so that the variance across cells is 1
- By default, only variable features are scaled
I will use all genes, not only variable features

```
all.genes <- rownames(so_wk4)
so_wk4 <- ScaleData(so_wk4, features = all.genes)
```
Just a reminder that these steps could replaced by SCTransform

# 3.4 Demultiplexing
In standard scRNA analysis, we can jump directly into part **3.5**, but since CITE-sequencing pools all cells from different samples, we need to demultiplex them based on the HTO tag in prior. 
We can start by normalize HTO data the using a function of Seurat to demultiplex:
```

so_wk4 <- HTODemux(so_wk4, assay = "HTO", positive.quantile = 0.99)
```
Based on their expression in HTO, cells will be classified into 3 categories :
- Doublet : expression of HTO is more than for one specific cell
- Singlet : significant expression of HTO corresponding to only one specific sample
- Negative : show anything significant expression of HTO
These data will be stored in HTO_classification.global
```
table(so_wk4$HTO_classification.global)
#Doublet Negative  Singlet 
# 1275      115    10004  
```
Visualize the classification 
```
png(paste0(OUT_DIR, "3_4_HTO_classification.png"), width = 1200, height = 1000, res = 150)
Idents(so_wk4) <- "HTO_classification.global"
VlnPlot(so_wk4, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()
```
![3_4_HTO_classification](https://github.com/user-attachments/assets/9f416a88-0c9a-4646-adc3-a0f95793542a)

We can also visualize the distribution of each HTO across cells, showing if we are distinct population of cells correspond to each sample 
```
png(paste0(OUT_DIR, "3_2_HTO_enrichment.png"), width = 1600, height = 1100, res = 150)
Idents(so_wk4) <- "HTO_maxID"
RidgePlot(so_wk4, assay = "HTO", features = rownames(so_wk4[["HTO"]])[1:6], ncol = 3)
dev.off()
```
![3_4_HTO_enrichment](https://github.com/user-attachments/assets/df53410a-9f52-4321-bd78-a937dcc1e26d)

Then we will only keep the singlet cells for further analysis 
```
so_wk4_doublet <- subset(so_wk4_doublet, subset = HTO_classification.global != "Negative")
so_wk4_doublet <- subset(so_wk4_doublet, subset = HTO_classification.global != "Doublet")
```
Also with the HTO classification, we can now determine the group of sample, whether is from old or young mice:
  - Young mice : HTO 1, 2, 3
  - Old mice : HTO 4, 5, 6
```
so_wk4@meta.data$condition <- NA
so_wk4@meta.data$condition[so_wk4@meta.data$HTO_classification %in% c("HTO1-TotalA", "HTO2-TotalA", "HTO3-TotalA")] <- "young"
so_wk4@meta.data$condition[is.na(so_wk4@meta.data$condition)] <- "old"
table(so_wk4@meta.data$condition)
#old young 
#5006  4998 
```
We can also save our data singlet filtered
```
saveRDS(so_wk4, file = "so_wk4.rds")
```
## 3.5 Dimension reduction
One important step in scRNA analysis is reduce the number of input features while retaining maximum information.[We can perform it by 2 techniques: feature selection and feature extraction (PCA, t-SNE)](https://medium.com/@aastha.code/dimensionality-reduction-pca-t-sne-and-umap-41d499da2df2).
- Principal Component Analysis (PCA) : linear dimension reduction algorithm
- t-Distributed Stochastic Neighbor Embedding (t-SNE) : non linear dimensionality reduction algorithm
- Unform Manifold Approximation and Projection (UMAP) : non linear dimensionality reduction algorithm

### 3.5.1 Features selection
Then highly variable features across cells will be choosen, basically they are highly expressed in some cells and lowly expressed in others. The default of features number is 2000 
```
so_wk4 <- FindVariableFeatures(so_wk4, selection.method = "vst", nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(so_wk4), 10) 
#Plot 
png(paste0(OUT_DIR, "3_3_Variable_features.png"), width = 2000, height = 2000, res = 150)
plot1 <- VariableFeaturePlot(so_wk4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()
```
![3_3_Variable_features](https://github.com/user-attachments/assets/4fa57ce1-6b0f-4bbc-8254-9977744d880d)

### 3.3.3 Linear dimensional reduction 
We can now apply PCA to our seurat object, only variable features will be used as imput 
```
so_wk4 <- RunPCA(so_wk4, features = VariableFeatures(object = so_wk4))
```
To visualize the result, we can use VizDimReduction(), DimPlot(), and DimHeatmap() 
```
print(so_wk4[["pca"]], dims = 1:5, nfeatures = 5)
PC_ 1 
Positive:  Ctsg, Mpo, Hp, Ptma, Anxa3 
Negative:  Angpt1, Meis1, Prkg1, Car2, Mecom 
PC_ 2 
Positive:  Ermap, Rhd, Slc25a21, Klf1, Car1 
Negative:  Prtn3, Arpc1b, Pkm, Tmsb4x, Napsa 
PC_ 3 
Positive:  Nusap1, Ckap2l, Lockd, Prc1, Mki67 
Negative:  Car1, Ermap, Slc25a21, Tspo2, Atp1b2 
PC_ 4 
Positive:  Il1rapl2, Gda, Selm, Tox, Plcl1 
Negative:  Slc22a3, Csrp3, Sox4, Tmsb10, St8sia6 
PC_ 5 
Positive:  Vim, Kif23, Ube2c, Nusap1, Ect2 
Negative:  Pf4, Pdcd4, Treml1, Itga2b, Serpine2 
```

or 
```
png(paste0(OUT_DIR, "3_3_PCA.png"), width = 2000, height = 1000, res = 150)
VizDimLoadings(so_wk4, dims = 1:2, reduction = "pca")
dev.off()

```
![3_3_PCA_viz](https://github.com/user-attachments/assets/b73c65a9-5f71-4ebc-b99d-5025c91a7128)

With DimHeatmap(), we can explore the primary sources of heterogeneity in our dataset 
```
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```
![DimHeatmap_PC1](https://github.com/user-attachments/assets/6d6d772b-0289-4cff-a952-c2dd97499c4d)
To determine the right dimension to use, we can generate an elbow plot, ranking of the principal components based on the percentage of variance explained by each one

ElbowPlot(pbmc)




