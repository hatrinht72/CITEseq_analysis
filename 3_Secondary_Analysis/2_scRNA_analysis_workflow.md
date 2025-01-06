# 3. Standard workflow for RNA data set
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
## 3.3 Features selection
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
![3_3_Variable_features](https://github.com/user-attachments/assets/40112665-db00-40d9-8274-2e762296b007)

## 3.4 Scaling the data
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


# 3.5 Demultiplexing
In standard scRNA analysis, we can jump directly into part **3.6**, but since CITE-sequencing pools all cells from different samples, we need to demultiplex them based on the HTO tag in prior. 
We can start by normalize HTO data the using a function of [Seurat](https://satijalab.org/seurat/articles/hashing_vignette) to demultiplex:
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
We can first identify which cell is come from which sample
```
so_wk4@meta.data <- so_wk4@meta.data %>%
  mutate(sample = case_when(
    HTO_classification == "HTO1-TotalA" ~ "young 1",
    HTO_classification == "HTO2-TotalA" ~ "young 2",
    HTO_classification == "HTO3-TotalA" ~ "young 3",
    HTO_classification == "HTO4-TotalA" ~ "old 1",
    HTO_classification == "HTO5-TotalA" ~ "old 2",
    HTO_classification == "HTO6-TotalA" ~ "old 3",
    TRUE ~ NA_character_  # If none of the conditions match, assign NA
  ))

#old 1   old 2   old 3 young 1 young 2 young 3 
#1874    1726    1406    1821    1712    1465 
```
Then we can set up the condition 
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
## 3.6 Dimension reduction
One important step in scRNA analysis is reduce the number of input features while retaining maximum information.[We can perform it by 2 techniques: feature selection and feature extraction (PCA, t-SNE)](https://medium.com/@aastha.code/dimensionality-reduction-pca-t-sne-and-umap-41d499da2df2).
- Principal Component Analysis (PCA) : linear dimension reduction algorithm
- t-Distributed Stochastic Neighbor Embedding (t-SNE) : non linear dimensionality reduction algorithm
- Unform Manifold Approximation and Projection (UMAP) : non linear dimensionality reduction algorithm

### 3.6.1 Linear dimensional reduction with PCA
####  3.6.1.1 Dimensional reduction
We can now apply PCA to our seurat object, only variable features will be used as imput 
```
so_wk4 <- RunPCA(so_wk4, features = VariableFeatures(object = so_wk4))
```
To visualize the result, we can use VizDimReduction(), DimPlot(), and DimHeatmap() 
```
print(so_wk4[["pca"]], dims = 1:5, nfeatures = 5)

PC_ 1 
Positive:  Ptma, Ctsg, Mpo, Hp, Anxa3 
Negative:  Angpt1, Meis1, Prkg1, Mecom, Car2 
PC_ 2 
Positive:  Ermap, Klf1, Slc25a21, Car1, Rhd 
Negative:  Prtn3, Arpc1b, Pkm, Napsa, Limd2 
PC_ 3 
Positive:  Car1, Slc25a21, Ermap, Tspo2, Atp1b2 
Negative:  Nusap1, Ckap2l, Lockd, Prc1, Fam64a 
PC_ 4 
Positive:  Slc22a3, Tmsb10, Sox4, Csrp3, St8sia6 
Negative:  Il1rapl2, Gda, Selm, Plcl1, Tox 
PC_ 5 
Positive:  Pf4, Pdcd4, Treml1, Itga2b, Serpine2 
Negative:  Kif23, Nusap1, Ube2c, Vim, Ect2 

```
or 
```
png(paste0(OUT_DIR, "3_3_PCA.png"), width = 2000, height = 1000, res = 150)
VizDimLoadings(so_wk4, dims = 1:2, reduction = "pca")
dev.off()

```
![3_6_PCA_vizdim](https://github.com/user-attachments/assets/df30c9ee-5f60-43d3-a1f0-38bb8602e4fa)

With DimHeatmap(), we can explore the primary sources of heterogeneity in our dataset 
```
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```
![DimHeatmap_PC1](https://github.com/user-attachments/assets/f7552fb0-a05f-40a9-a9f7-8bd28cd5a7bd)

To determine the right dimension to use, we can generate an elbow plot, ranking of the principal components based on the percentage of variance explained by each one

```
png(paste0(OUT_DIR, "3_6_PCA_elbow.png"), width = 2000, height = 2000, res = 150)
ElbowPlot(so_wk4)
dev.off()
```
![3_6_PCA_elbow](https://github.com/user-attachments/assets/2b6cb259-5de5-47a7-a24c-3b386983164a)

So according to the plot that we have, we can chose 15 as a cutoff

####  3.6.1.2 Cluster the cells
Based on the cutoff that we choosen, now we can cluster the cell 
```
so_wk4 <- FindNeighbors(so_wk4, dims = 1:15)
so_wk4 <- FindClusters(so_wk4, resolution = 0.5)
```
We can look at cluster IDs of the first 5 cells :
```
head(Idents(so_wk4), 5)

AAACCCAAGCAGGCAT-1 AAACCCAAGTATGATG-1 AAACCCAAGTCGAATA-1 AAACCCACAACAGCTT-1 
                 2                  0                  0                  7 
AAACCCACAGCTCTGG-1 
                 4 
Levels: 0 1 2 3 4 5 6 7 8 9 10
```

### 3.6.2 Non-linear dimensional reduction (UMAP/tSNE)
We can now try to create UMAP plot or t-SNE, based on the PCA dimensions
```
so_wk4 <- RunUMAP(so_wk4, reduction = "pca", dims = 1:15, assay = "RNA", reduction.name = "rna.umap")
png(paste0(OUT_DIR, "3_6_UMAP_PCA_RNA.png"), width = 2000, height = 1000, res = 150)
plot1 <- DimPlot(so_wk4, reduction = "rna.umap", group.by = "seurat_clusters")
plot2 <- DimPlot(so_wk4, reduction = "rna.umap", group.by = "condition")
plot1 + plot2
dev.off()
```
![3_6_UMAP_PCA_RNA](https://github.com/user-attachments/assets/a12c3933-e37c-4c03-ae10-143d08cb7e8f)

```
so_wk4 <- RunTSNE(so_wk4, reduction = "pca", dims = 1:15, assay = "RNA", reduction.name = "rna.TSNE")
png(paste0(OUT_DIR, "3_6_TSNE_PCA_RNA_dims30.png"), width = 2000, height = 1000, res = 150)
plot1 <- DimPlot(so_wk4, reduction = "rna.T.SNE", group.by = "seurat_clusters")
plot2 <- DimPlot(so_wk4, reduction = "rna.T.SNE", group.by = "condition")
plot1 + plot2
dev.off()
```
![3_6_TSNE_PCA_RNA_dims30](https://github.com/user-attachments/assets/6ab8ae84-5a76-4df2-8bb1-353c0641f870)


We can extract the number of cell in each cluster, each sample  
```
n_cells <- FetchData(so_wk4, 
                     vars = c("ident", "orig.ident")) %>%
        dplyr::count(ident, orig.ident) %>%
        tidyr::spread(ident, n)
View(n_cells)
     orig.ident    0    1    2    3    4   5   6   7   8   9 10
1 SeuratProject 1907 1523 1447 1084 1005 901 726 614 575 183 39

```
### 3.7 Differential expression analysis
Now we have clusters, we can extract the markers of each cluster, compared to one other cluster or all others clusters. These information could be useful if we want to annotate cell type.
```
Idents(so_wk4) <- "seurat_clusters"
rna_markers <- FindAllMarkers(so_wk4, assay = "RNA")
saveRDS(rna_markers, file = "rna_markers.rds")
```
Based on the markers that we found, we can identify the cell type corresponding each cluster 




In the case of trying this article, we can also try to find differential expression between old and young mice. 

In seurat, you can choose the Identity class labels of cells that you would like to use, for example with Seurat.cluster identification, each cell will be assigned to the cluster to whom they belong.
We can change the Identity class labels with this command : 
```
Idents(so_wk4) <- "seurat_clusters"
```
We can now try to find differential expression gene between old and young mice across clusters
```
so_wk4$celltype.condition <- paste(Idents(so_wk4), so_wk4$condition, sep="_")
so_wk4$celltype <- Idents(so_wk4)
Idents(so_wk4) <- "celltype.condition"

for (i in 0:9){ #9 or however many clusters you have
try({
ident1 <- paste0(i,"_old")
ident2 <- paste0(i,"_young")
condition.diffgenes <- FindMarkers(so_wk4, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
write.csv(condition.diffgenes, file=paste0(output_diff,i,".csv"))
})
}
```
