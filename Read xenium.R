1) Paths and sample labels (edit to your directories)
paths <- c(
"path/to/sample1",
"path/to/sample2",
"path/to/sample3",
"path/to/sample4",
"path/to/sample5",
"path/to/sample6",
"path/to/sample7",
"path/to/sample8",
"path/to/sample9",
"path/to/sample10",
"path/to/sample11"
)
samples <- paste0("sample", seq_along(paths))

2) Loader for Xenium output directories
LoadXenium2 <- function(data.dir, fov = "fov", assay = "Xenium") {

Read Xenium counts and spatial objects
data <- ReadXenium(
data.dir = data.dir,
type = c("centroids", "segmentations")
)
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
fov_obj <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)

Create Seurat object with gene counts and attach control assays
obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
obj[[fov]] <- fov_obj
obj
}

3) Read all runs and add sample labels
objs <- lapply(seq_along(paths), function(i) {
x <- LoadXenium2(paths[i], fov = "fov", assay = "Xenium")
x$sample <- samples[i]
x
})

4) Merge all into one Seurat object
xenium.obj <- Reduce(function(a, b) merge(a, y = b, add.cell.ids = c(a$sample[1], b$sample[1])), objs)

Optional: free memory
rm(objs); gc()

5) Basic QC (adjust thresholds as needed)
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 10 & nFeature_Xenium > 0)

6) SCTransform (SCT) + PCA + Harmony integration
xenium.obj <- SCTransform(
xenium.obj, assay = "Xenium",
conserve.memory = TRUE,
ncells = 2000, # increase if RAM allows
variable.features.n = 1000 # adjust feature count
)
xenium.obj <- RunPCA(xenium.obj, npcs = 30, reduction.name = "pca")

Harmony-based layer integration (Seurat v5 IntegrateLayers)
xenium.obj <- IntegrateLayers(
object = xenium.obj,
method = HarmonyIntegration,
assay = "SCT",
orig.reduction = "pca",
new.reduction = "integrated.cca",
verbose = TRUE
)

7) UMAP, neighbors, clustering
xenium.obj <- RunUMAP(xenium.obj, reduction = "integrated.cca", dims = 1:30, n.neighbors = 20, min.dist = 0.3)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "integrated.cca", dims = 1:30, k.param = 20)
xenium.obj <- FindClusters(xenium.obj, resolution = 1.0)

8) Quick visualizations
p_clusters <- DimPlot(xenium.obj, label = TRUE, repel = TRUE, raster = TRUE) + ggtitle("Clusters")
p_samples <- DimPlot(xenium.obj, group.by = "sample", raster = TRUE) + ggtitle("By sample")
print(p_clusters | p_samples)

9) Fix umi.assay in SCT models (needed for PrepSCTFindMarkers/FindMarkers sometimes)
if (!is.null(xenium.obj@assays$SCT@SCTModel.list)) {
for (i in seq_along(xenium.obj@assays$SCT@SCTModel.list)) {
slot(object = xenium.obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
}
}

10) Prep and compute markers
xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- PrepSCTFindMarkers(xenium.obj, assay = "SCT")
Idents(xenium.obj) <- "seurat_clusters"

markers_all <- FindAllMarkers(
xenium.obj,
assay = "SCT",
only.pos = TRUE,
logfc.threshold = 0.25,
min.pct = 0.1,
recorrect_umi = FALSE
)

Top markers per cluster (e.g., top 3 by adjusted p-value)
top_markers <- markers_all %>%
group_by(cluster) %>%
slice_min(n = 3, order_by = p_val_adj) %>%
pull(gene) %>%
unique()

11) DotPlot of top markers
p_dot <- DotPlot(xenium.obj, features = top_markers, assay = "SCT") +
scale_color_viridis(option = "magma") +
ylab("Clusters") + xlab("Marker genes") +
theme(axis.text.x = element_text(size = 8, angle = 90),
axis.text.y = element_text(size = 8))
print(p_dot)

12) Spatial plots (auto-select first FOV; guard against missing)
fov_names <- names(Filter(function(x) inherits(x, "FOV"), xenium.obj@images))
if (length(fov_names) > 0) {
DefaultFOV(xenium.obj) <- fov_names[1]
DefaultBoundary(xenium.obj[[fov_names[1]]]) <- "centroids"

Example FeaturePlot on tissue
p_feat <- ImageFeaturePlot(
xenium.obj,
features = "CXCL9", # change gene as needed
max.cutoff = "q95",
size = 1.5,
dark.background = FALSE
) + ggtitle("CXCL9 (spatial)")
print(p_feat)
}

13) Optional: cluster composition heatmap by sample or by another meta (if present)
if ("sample" %in% colnames(xenium.obj@meta.data)) {
mat_counts <- table(xenium.obj$seurat_clusters, xenium.obj$sample)
Heatmap(
mat_counts,
name = "count",
cluster_rows = TRUE, cluster_columns = TRUE,
show_row_dend = FALSE, show_column_dend = FALSE
)
}
