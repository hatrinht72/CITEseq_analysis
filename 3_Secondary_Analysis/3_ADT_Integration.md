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
We can also find the cluster based on ADT information 
```
so_wk4 <- FindNeighbors(so_wk4, assay = "ADT", reduction ="adt.pca", dims = 1:25)
so_wk4 <- FindClusters(so_wk4, assay = "ADT", redection.type = "adt.pca", dims.use = 1:25,resolution = 0.8)
so_wk4@meta.data$adt_clusters <- so_wk4@meta.data$RNA_snn_res.0.8
```
We can then do dimension reduction like UMAP 

```
so_wk4 <- RunUMAP(so_wk4, reduction = "adt.pca", dims = 1:25, assay = "ADT", reduction.name = "adt.umap")
png(paste0(OUT_DIR, "4_1_UMAP_PCA_ADT.png"), width = 3000, height = 1000, res = 150)
plot0 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "adt_clusters")
plot1 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "rna_clusters")
plot2 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "condition")
plot0 + plot1 + plot2
dev.off()
```
![4_1_UMAP_PCA_ADT](https://github.com/user-attachments/assets/a2c0f44c-9d35-472b-81b4-c0426033e3d7)

Here we can see we have some similarity in clusters based on ADT and cluster based on RNA but still remain some difference, thats why we can move into next part, integration these 2 datasets.

## 4.2 ADT integration
Now we already have our ADT and RNA dimensional reduction, we can use the integration tools [Weighted Nearest Neighbor Analysis](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis) proposed by Seurat 

```
so_wk4 <- FindMultiModalNeighbors(
  so_wk4, 
  reduction.list = list("rna.pca", "adt.pca"), 
  dims.list = list(1:30, 1:25), 
  modality.weight.name = c("RNA.weight", "ADT.weight")
)
so_wk4 <- FindClusters(so_wk4, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
so_wk4@meta.data$wnn_clusters <- so_wk4@meta.data$wsnn_res.2
so_wk4 <- RunUMAP(so_wk4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
```
We can have other clusters based on this integration 
So now we can compare all 3 reductions that we made 
```
png(paste0(OUT_DIR, "4_2_UMAP_all.png"), width = 3000, height = 1000, res = 150)
plot0 <-DimPlot(so_wk4, reduction = "rna.umap", group.by = "rna_clusters")
plot1 <-DimPlot(so_wk4, reduction = "adt.umap", group.by = "adt_clusters")
plot2 <-DimPlot(so_wk4, reduction = "wnn.umap", group.by = "wnn_clusters")
plot0 + plot1 + plot2
dev.off()
```
![4_2_UMAP_all](https://github.com/user-attachments/assets/dd4cf7b7-2b5e-41d7-b13f-ca178c140497)

only WNN reduction version 
```
png(paste0(OUT_DIR, "4_2_UMAP_wnn.png"), width = 4000, height = 2500, res = 150)
plot0 <-DimPlot(so_wk4, reduction = "wnn.umap", group.by = "rna_clusters")
plot1 <-DimPlot(so_wk4, reduction = "wnn.umap", group.by = "adt_clusters")
plot2 <-DimPlot(so_wk4, reduction = "wnn.umap", group.by = "wnn_clusters")
plot3 <-DimPlot(so_wk4, reduction = "wnn.umap", group.by = "condition")
plot0+plot1+plot2+plot3
dev.off()
```
![4_2_UMAP_wnn](https://github.com/user-attachments/assets/ba23f952-36f5-4ccc-96ff-1f208f90bee9)

Since the clusters is difference, so we can find new set of markers, a reminder : If we dont precise assay, it will automatically searching on ADT

```
Idents(so_wk4) <- "wnn_clusters"
wnn_rna_markers <- FindAllMarkers(so_wk4, assay = "RNA")

wnn_all_pos_markers = FindAllMarkers(object = so_wk4, assay= "RNA",
                                     only.pos = TRUE, # genes more expressed in the cluster compared
                                     min.pct = 0.25, # % of cell expressing the marker
                                     logfc.threshold = 0.25)

saveRDS(wnn_all_pos_markers, file = "wnn_all_pos_markers.rds")
saveRDS(wnn_rna_markers, file = "wnn_rna_markers.rds")
```
So we can based on these integration and litterature to annotate our cells

