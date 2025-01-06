# 4. ADT integration
## 4.1 [ADT analysis : Normalization, scale, dimensional reduction](https://satijalab.org/seurat/articles/multimodal_vignette.html#identify-cell-surface-markers-for-scrna-seq-clusters-1)
After scRNA analysis, we can now attack to the part of ADT. 

We also have to normalize and scale these data and we will use all ADT features for dimensional reduction 

```
adtfeatures <- rownames(so_wk4[["ADT"]])
VariableFeatures(so_wk4) <- adtfeatures
so_wk4 <- NormalizeData(so_wk4, normalization.method = 'CLR', margin = 2) 
so_wk4 <- ScaleData(so_wk4, features = adtfeatures)
so_wk4 <- RunPCA(so_wk4, features = VariableFeatures(object = so_wk4), reduction.name ='adt.pca')
```
We can then do dimension reduction like UMAP 

```
so_wk4 <- RunUMAP(so_wk4, reduction = "adt.pca", dims = 1:25, assay = "ADT", reduction.name = "adt.umap")
png(paste0(OUT_DIR, "6_Cluster_by_ADT_weight.png"), width = 2000, height = 1000, res = 150)
plot1 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "seurat_clusters")
plot2 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "condition")
plot1 + plot2
dev.off()
```
![4_1_Cluster_by_ADT_weight](https://github.com/user-attachments/assets/50cbe942-0a27-43a1-bbcb-8500333b30ae)

## 4.2 ADT integration
Now we already have our ADT and RNA dimensional reduction, we can use the integration tools [Weighted Nearest Neighbor Analysis](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis) proposed by Seurat 

```
so_wk4 <- FindMultiModalNeighbors(
  so_wk4, 
  reduction.list = list("rna.pca", "adt.pca"), 
  dims.list = list(1:30, 1:25), 
  modality.weight.name = c("RNA.weight", "ADT.weight")
)
so_wk4 <- RunUMAP(so_wk4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
```
We can have other clusters based on this integration 
```
so_wk4 <- FindClusters(so_wk4, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
```
![4_2_UMAP_WNN](https://github.com/user-attachments/assets/595819e4-6765-4199-a788-48a31827efba)

So we can based on these integration and litterature to annotate our cells

