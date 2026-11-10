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

# Load the pathologist annotations
annotations <- read.csv(
  gzfile("data/GSE213688/GSM6592048_M1_pathologist_annotations.csv.gz")
)

# Match annotations to spots by barcode
annotations$Annotations <- trimws(as.character(annotations$Annotations))
annotations$Annotations[annotations$Annotations == ""] <- NA

obj$pathologist_annotation <- annotations$Annotations[
  match(colnames(obj), annotations$Barcode)
]

# Compare clusters with pathologist annotations
tab <- table(
  cluster = Idents(obj),
  pathology = obj$pathologist_annotation,
  useNA = "ifany"
)

# Percentages within each cluster
round(100 * prop.table(tab, margin = 1), 1)

# Percentages within each pathology category
round(100 * prop.table(tab, margin = 2), 1)

# Plot pathologist annotations on the tissue
obj$pathology_plot <- as.character(obj$pathologist_annotation)
obj$pathology_plot[is.na(obj$pathology_plot)] <- "Unannotated"

SpatialDimPlot(
  obj,
  group.by = "pathology_plot",
  label = FALSE
)
label_map <- c(
  "0" = "ECM-rich mixed stroma",
  "1" = "ECM-rich mixed stroma",
  "2" = "Ig-rich mixed stroma",
  "3" = "SCGB2A2-high, fibrosis-rich",
  "4" = "Mixed",
  "5" = "Epithelium-associated",
  "6" = "Adipose-associated"
)

cluster_id <- as.character(Idents(obj))
labels <- unname(label_map[cluster_id])
names(labels) <- names(cluster_id)

obj$tentative_label <- labels

saveRDS(obj, "M1_seurat_annotated.rds")
write.csv(markers, "M1_markers.csv", row.names = FALSE)

# Recreate the plot without large labels
SpatialDimPlot(obj, group.by = "tentative_label", label = FALSE)
