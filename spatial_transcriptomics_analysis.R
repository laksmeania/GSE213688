library(Seurat)
library(dplyr)

# Load the M1 data
obj <- Load10X_Spatial(data.dir = "data/GSE213688/M1")
obj <- NormalizeData(obj)

# Find variable genes, then run PCA and clustering
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

# Find Marker
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Provisional labels based on marker genes
label_map <- c(
  "0" = "ECM-rich",
  "1" = "ECM-rich",
  "2" = "Ig-rich",
  "3" = "SCGB2A2-high",
  "4" = "Mixed immune",
  "5" = "Unclear",
  "6" = "ADIPOQ-high"
)

cluster_id <- as.character(Idents(obj))
labels <- unname(label_map[cluster_id])
names(labels) <- names(cluster_id)
obj$tentative_label <- labels

# Save results
saveRDS(obj, "M1_seurat_annotated.rds")
write.csv(markers, "M1_markers.csv", row.names = FALSE)

# Recreate the plot without large labels
SpatialDimPlot(obj, group.by = "tentative_label", label = FALSE)
